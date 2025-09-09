import html
from bs4 import BeautifulSoup
from jinja2 import Template
from owslib.iso import namespaces
from requests_toolbelt.multipart.encoder import MultipartEncoder
from shapely.geometry import LineString
from urllib.parse import urljoin, urlparse
from urllib.request import urlopen
import codecs
import os
import pandas as pd
import re
import requests
import shapely
import stamina
import urllib.error
import xarray as xr
import xml.etree.ElementTree as ET
import zipfile


@stamina.retry(on=urllib.error.HTTPError, attempts=3)
def _openurl_with_retry(url):
    """Thin wrapper around urlopen adding stamina."""
    return urlopen(url)


def get_ncei_accession_mapping():
    """
    Scrapes NCEI for ATN accession numbers and associated file metadata.

    Returns:
        pd.DataFrame: A DataFrame mapping accession numbers to file names,
                      download URLs, and other metadata.
    """
    print("Fetching NCEI accession mapping table...")
    # Configure namespaces for XML parsing
    namespaces.update({"gmi": "http://www.isotc211.org/2005/gmi"})
    namespaces.update({"gml": "http://www.opengis.net/gml/3.2"})
    if None in namespaces:
        del namespaces[None]

    url = "https://www.ncei.noaa.gov/access/metadata/landing-page/bin/iso?id=gov.noaa.nodc:IOOS-ATN-STP;view=xml;responseType=text/xml"
    iso = _openurl_with_retry(url)
    iso_tree = ET.parse(iso)
    root = iso_tree.getroot()

    accessions = []
    # Collect individual accession IDs
    for MD_keywords in root.iterfind(
        ".//gmd:descriptiveKeywords/gmd:MD_Keywords", namespaces
    ):
        for thesaurus_name in MD_keywords.iterfind(
            ".//gmd:thesaurusName/gmd:CI_Citation/gmd:title/gco:CharacterString",
            namespaces,
        ):
            if thesaurus_name.text == "NCEI ACCESSION NUMBER":
                for acce_no in MD_keywords.iterfind(
                    ".//gmd:keyword/gmx:Anchor", namespaces
                ):
                    accessions.append(acce_no.text)

    df_map = pd.DataFrame()

    for acc in accessions:
        try:
            ## There are a couple endpoints we can use
            #
            # https://www.ncei.noaa.gov/data/oceans/ncei/archive/metadata/approved/granule/{acc}.xml
            # https://www.ncei.noaa.gov/access/metadata/landing-page/bin/iso?id=gov.noaa.nodc:{acc};view=xml
            # http://www.ncei.noaa.gov/metadata/granule/geoportal/rest/metadata/item/IOOS-ATN-STP.{acc}/xml
            #
            url = f"https://www.ncei.noaa.gov/data/oceans/ncei/archive/metadata/approved/granule/{acc}.xml"
            iso = _openurl_with_retry(url)
            iso_tree = ET.parse(iso)
            root = iso_tree.getroot()

            # Collect terms of interest.
            title = pd.DataFrame(
                {
                    "title": [
                        root.find(".//gmd:title/gco:CharacterString", namespaces).text
                    ],
                    "accession": acc,
                    "start_date": [
                        root.find(
                            ".//gml:TimePeriod/gml:beginPosition", namespaces
                        ).text
                    ],
                    "end_date": [
                        root.find(".//gml:TimePeriod/gml:endPosition", namespaces).text
                    ],
                }
            )

            # Collect DataOne References
            for MD_AggregateInformation in root.iterfind(
                ".//gmd:MD_AggregateInformation", namespaces
            ):
                if (
                    MD_AggregateInformation.find(
                        ".//gmd:description/gco:CharacterString", namespaces
                    ).text
                    == "related dataset"
                ):
                    related_data = MD_AggregateInformation.find(
                        ".//gmd:linkage/gmd:URL", namespaces
                    ).text
                    related_data_citation = MD_AggregateInformation.find(
                        "./gmd:aggregateDataSetName/gmd:CI_Citation/gmd:title/gco:CharacterString",
                        namespaces,
                    ).text
                    #print(f"{acc} = {related_data} {related_data_citation}")
                    title["related_data_url"] = related_data
                    title["related_data_citation"] = html.escape(related_data_citation)

            for CI_OnlineResource in root.iterfind(
                ".//gmd:onLine/gmd:CI_OnlineResource", namespaces
            ):
                if (
                    CI_OnlineResource.find(
                        ".//gmd:protocol/gco:CharacterString", namespaces
                    ).text
                    == "FTP"
                ):
                    string = CI_OnlineResource.find(
                        ".//gmd:linkage/gmd:URL", namespaces
                    ).text
                    arc = re.search("(arc[0-9]{1,4})", string)
                    if arc:
                        title["arc"] = [arc.group()]
                        xml_manifest = f"https://www.ncei.noaa.gov/data/oceans/archive/{arc.group()}/{acc}/{acc}.1.1.xml"
                        title["xml"] = [xml_manifest]

                        iso_mani = _openurl_with_retry(xml_manifest)
                        iso_mani_tree = ET.parse(iso_mani)
                        root_mani = iso_mani_tree.getroot()

                        for path in root_mani.iterfind(".//path"):
                            fpath = path.text
                            fname = re.search("(data/0-data/atn_.*)", fpath)
                            if fname:
                                title["file_name"] = fname.group().split("/")[-1]

            df_map = pd.concat([df_map, title], ignore_index=True)

        except Exception as e:
            print(f"Could not process accession {acc}: {e}")

    df_map["ptt_id"] = df_map["title"].str.extract(r".*ptt ([0-9]{3,7}) .*")
    print(f"Successfully built mapping for {len(df_map)} files.")
    return df_map


def recursive_wget(url, output_dir):
    """
    Recursively downloads files from a given URL to a specified output directory,
    mirroring the directory structure of the website.

    Args:
        url (str): The URL to start downloading from.
        output_dir (str): The local directory to save files to.
    """
    print(f"Accessing: {url}")
    try:
        # --- Create the output directory if it doesn't exist ---
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
            print(f"Created directory: {output_dir}")

        # --- Send a GET request and parse the HTML ---
        response = requests.get(url)
        # Raise an exception for bad status codes (like 404 Not Found)
        response.raise_for_status()
        soup = BeautifulSoup(response.text, "html.parser")

        # --- Find all links on the page ---
        for link in soup.find_all("a"):
            href = link.get("href")

            # --- Skip invalid or parent directory links ---
            if not href or href.startswith("?") or href.startswith("/") or ".." in href:
                continue

            # --- Construct the full, absolute URL for the link ---
            absolute_url = urljoin(url, href)

            # Get the path component of the URL to create local directories/files
            path = urlparse(absolute_url).path
            # Create a valid local path from the last part of the URL path
            local_path = os.path.join(output_dir, os.path.basename(path))

            # --- If the link points to a directory, recurse into it ---
            if href.endswith("/"):
                print(f"\nEntering directory: {absolute_url}")
                # Call the function again for the new directory
                recursive_wget(absolute_url, local_path)
            # --- If the link points to a file, download it ---
            else:
                download_file(absolute_url, local_path)

    except requests.exceptions.HTTPError as e:
        print(f"HTTP Error accessing URL {url}: {e}")
    except requests.exceptions.RequestException as e:
        print(f"Error accessing URL {url}: {e}")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")


def download_file(url, local_path):
    """
    Downloads a single file from a URL and saves it to a local path.

    Args:
        url (str): The URL of the file to download.
        local_path (str): The local path where the file will be saved.
    """
    try:
        print(f"  Downloading file: {os.path.basename(local_path)}")
        # Use stream=True to efficiently download large files
        with requests.get(url, stream=True) as r:
            r.raise_for_status()
            # Open the file in binary write mode
            with open(local_path, "wb") as f:
                # Write the file in chunks
                for chunk in r.iter_content(chunk_size=8192):
                    f.write(chunk)
        # print(f"  Successfully downloaded {os.path.basename(local_path)}")
    except requests.exceptions.RequestException as e:
        print(f"  Failed to download {url}: {e}")
    except IOError as e:
        print(f"  Failed to write file {local_path}: {e}")


def create_dwc_occurrence(ds: xr.Dataset, output_csv: str, df_map: pd.DataFrame):
    """Create a Darwin Core Occurrence CSV from an xarray Dataset."""
    source_file = os.path.basename(ds.encoding.get("source"))
    # bail if we can't find the file in the mapping table.
    if source_file not in df_map["file_name"].values:
        raise KeyError(f"File {source_file} not found in NCEI Accession mapping table.")

    filename = os.path.splitext(source_file)[
        0
    ]  # "ioos_atn_{ds.ptt_id}_{start_date}_{end_date}""

    file_map_entry = df_map[df_map["file_name"] == source_file].iloc[0]

    acce_no = file_map_entry["accession"]
    related_data_url = file_map_entry["related_data_url"]

    dwc_df = pd.DataFrame()
    dwc_df["occurrenceID"] = (
        "ioos_atn_"
        + ds.ptt_id
        + "_"
        + ds["time"].dt.strftime("%Y-%m-%dT%H:%M:%SZ")
        + "_"
        + ds["z"].astype(str)
        + "_"
        + ds.animal_common_name.replace(" ", "_")
    )
    dwc_df["eventID"] = filename
    dwc_df["organismID"] = (
        ds.platform_id + "_" + ds.animal_common_name.replace(" ", "_")
    )
    dwc_df["occurrenceStatus"] = "present"
    dwc_df["basisOfRecord"] = ds["type"]
    dwc_df["eventDate"] = ds["time"].dt.strftime("%Y-%m-%dT%H:%M:%SZ")
    dwc_df["decimalLatitude"] = ds["lat"]
    dwc_df["decimalLongitude"] = ds["lon"]
    dwc_df["geodeticDatum"] = ds.crs.epsg_code
    dwc_df["scientificName"] = ds["taxon_name"].values.tolist()
    dwc_df["scientificNameID"] = ds["taxon_lsid"].values.tolist()
    dwc_df["samplingProtocol"] = "satellite telemetry"
    dwc_df["kingdom"] = ds["animal"].attrs["kingdom"]
    dwc_df["taxonRank"] = ds["animal"].attrs["rank"]
    dwc_df["lifeStage"] = ds["animal_life_stage"].values.tolist()
    dwc_df["sex"] = ds["animal_sex"].values.tolist()
    dwc_df["associatedReferences"] = (
        f"https://doi.org/10.25921/wp4e-ph20; https://www.ncei.noaa.gov/archive/accession/{acce_no}; {related_data_url}"
    )
    dwc_df["minimumDepthInMeters"] = ds["z"].values.tolist()
    dwc_df["maximumDepthInMeters"] = ds["z"].values.tolist()
    dwc_df["bibliographicCitation"] = ds.citation

    # set basisOfRecord
    dwc_df.loc[dwc_df["basisOfRecord"] == "User", "basisOfRecord"] = "HumanObservation"
    dwc_df.loc[dwc_df["basisOfRecord"] == "Argos", "basisOfRecord"] = (
        "MachineObservation"
    )
    dwc_df.loc[dwc_df["basisOfRecord"] == "FastGPS", "basisOfRecord"] = (
        "MachineObservation"
    )

    # filter to respectable locations
    # drop A, B, and Z records
    dwc_df["location_class"] = ds["location_class"].to_series()
    valid_locations = ~dwc_df["location_class"].isin(["A", "B", "Z"])
    dwc_df = dwc_df[valid_locations].copy()

    print(f"  Extracted {len(dwc_df)} occurrences with valid locations.")

    # Map location class to coordinate uncertainty
    uncertainty_map = {"nan": 0, "G": 200, "3": 250, "2": 500, "1": 1500, "0": 10000}
    dwc_df["coordinateUncertaintyInMeters"] = dwc_df["location_class"].map(
        uncertainty_map
    )

    # --- Define Occurrences: First detection per location per hour ---
    dwc_df["event_hour"] = pd.to_datetime(dwc_df["eventDate"]).dt.strftime(
        "%Y-%m-%dT%H"
    )
    dwc_df.sort_values("event_hour", inplace=True)
    duplicate_counts = dwc_df.groupby(by="event_hour").transform("size")
    dwc_df["dataGeneralizations"] = (
        "first of " + duplicate_counts.astype(str) + " records for this hour."
    )
    dwc_df.loc[
        dwc_df["dataGeneralizations"] == "first of 1 records for this hour.",
        "dataGeneralizations",
    ] = ""
    dwc_df = dwc_df.drop_duplicates(subset=["event_hour"], keep="first").copy()

    # --- Add Occurrence Remarks ---
    dwc_df["occurrenceRemarks"] = (
        f"This is a representative occurrence from a full deployment. For the complete dataset please see https://www.ncei.noaa.gov/archive/accession/{acce_no}."
    )

    print(f"  Extracted {len(dwc_df)} occurrences to first row in hour.")

    # only pick specific columns to save
    cols = [
        "eventID",
        "occurrenceID",
        "occurrenceStatus",
        "basisOfRecord",
        "organismID",
        "eventDate",
        "decimalLatitude",
        "decimalLongitude",
        "geodeticDatum",
        "scientificName",
        "scientificNameID",
        "samplingProtocol",
        "kingdom",
        "taxonRank",
        "lifeStage",
        "sex",
        "associatedReferences",
        "coordinateUncertaintyInMeters",
        "minimumDepthInMeters",
        "maximumDepthInMeters",
        "dataGeneralizations",
        "bibliographicCitation",
        "occurrenceRemarks",
    ]

    # Save the individual CSV
    dwc_df.to_csv(output_csv, columns=cols, index=False)
    print(f"  Saved data to '{output_csv}'")

    return dwc_df, cols


def create_dwc_event(ds: xr.Dataset, dwc_df: pd.DataFrame, output_csv: str):
    # create parent event that is a summary of dwc_df
    event_df = pd.DataFrame()
    event_df["eventID"] = dwc_df["eventID"].unique()
    event_df["eventDate"] = dwc_df["eventDate"].min() + "/" + dwc_df["eventDate"].max()

    ## Convex hull summary of the points
    points = list(zip(dwc_df["decimalLongitude"], dwc_df["decimalLatitude"]))
    event_df["footprintWKT"] = shapely.convex_hull(LineString(points))

    event_df["minimumDepthInMeters"] = dwc_df["minimumDepthInMeters"].min()
    event_df["maximumDepthInMeters"] = dwc_df["maximumDepthInMeters"].max()
    event_df["eventType"] = "deployment"
    event_df["countryCode"] = "US"
    event_df["samplingProtocol"] = "satellite telemetry"
    event_df["dynamicProperties"] = [str(ds.attrs)]

    event_df.to_csv(output_csv.replace("occurrence", "event"), index=False)
    print(f"  Created {len(event_df)} events.")
    print(f"  Saved data to {output_csv.replace('occurrence', 'event')}")

    return event_df


def save_eml_file(eml_metadata: dict) -> str:
    """
    Save EML dictionary in a file
    Author: Jon Pye, Angela Dini
    Maintainer: Angela Dini
    :param eml_metadata: dictionary of EML metadata
    :return: filepath of where the EML filepath will be
    """
    # Write it out to the package
    template_file = codecs.open("templates/eml.xml.j2", "r", "UTF-8").read()
    template = Template(template_file)
    result_string = template.render(eml_metadata)
    eml_file = "data/dwc/{filename}/eml.xml".format(**eml_metadata)
    fh = codecs.open(eml_file, "wb+", "UTF-8")
    fh.write(result_string)
    fh.close()
    eml_full_path = os.path.abspath(eml_file)
    print(f"  EML metadata has been written to '{eml_full_path}'.")
    return eml_full_path


def create_eml(ds: xr.Dataset, df_map: pd.DataFrame):
    eml_metadata = ds.attrs
    source_file = os.path.basename(ds.encoding.get("source"))
    file_map_entry = df_map[df_map["file_name"] == source_file].iloc[0]

    contributors = dict()
    for attr in [
        x for x in ds.attrs if re.match(r"contributor_(?!role_vocabulary\b).*", x)
    ]:
        contributors[attr] = ds.attrs[attr].split(",")

    contributors_list = [
        {key: contributors[key][i] for key in contributors}
        for i in range(len(next(iter(contributors.values()))))
    ]

    other_meta = {
        "dataset_ipt_id": None,
        "dataset_short_name": ds.encoding.get("source")
        .split("\\")[-1]
        .replace(".nc", ""),
        "data_manager_firstname": "Megan",
        "data_manager_lastname": "McKinzie",
        "data_manager_title": "Data Manager",
        "data_manager_phone": "",
        "data_manager_email": "mmckinzie@mbari.org",
        "contributors": contributors_list,
        "ncei_accession_number": file_map_entry[
            "accession"
        ],  # df_map.loc[df_map['file_name'] == ds.encoding.get('source').split("\\")[-1], 'accession'].values[0],
        "related_data_url": file_map_entry["related_data_url"],
        "related_data_citation": file_map_entry["related_data_citation"],
        "ncei_title": file_map_entry[
            "title"
        ],  # df_map.loc[df_map['file_name'] == ds.encoding.get('source').split("\\")[-1], 'title'].values[0],
        "filename": os.path.splitext(source_file)[0],
        "nc_globals": str(ds.attrs),
    }

    eml_metadata.update(other_meta)

    instrument_info = ds["instrument_tag"].attrs

    eml_metadata.update(instrument_info)

    save_eml_file(eml_metadata)

    return eml_metadata


def create_dwc_emof(ds: xr.Dataset, dwc_df: pd.DataFrame, output_csv: str):
    # --- Creates the DwC Extended Measurement or Fact (eMoF) file. ---
    vars = list(ds.keys())
    animal_vars = [x for x in vars if re.match(r"animal_(?!life_stage\b|sex\b).*", x)]
    new_rows = pd.DataFrame()

    emof_ids = {
        "animal_weight": "http://vocab.nerc.ac.uk/collection/MVB/current/MVB000019",
        "animal_length": "http://vocab.nerc.ac.uk/collection/P01/current/TL01XX01/",
    }

    emof_unit_ids = {
        "kg": "http://vocab.nerc.ac.uk/collection/P06/current/KGXX/",
        "cm": "http://vocab.nerc.ac.uk/collection/P06/current/ULCM/",
    }

    for animal_var in animal_vars:
        row = pd.DataFrame(
            {
                "measurementValue": ds[animal_var].values.tolist(),
                "measurementType": [f"{animal_var}: {ds[animal_var].long_name}"],
                "measurementTypeID": [
                    emof_ids[animal_var] if animal_var in emof_ids.keys() else ""
                ],
                "measurementMethod": ds[animal_var].attrs[animal_var],
                "measurementUnit": [
                    ds[animal_var].units if "units" in ds[animal_var].attrs else ""
                ],
                "measurementUnitID": [
                    emof_unit_ids[ds[animal_var].units]
                    if ds[animal_var].units in emof_unit_ids.keys()
                    else ""
                ],
            }
        )
        new_rows = pd.concat([new_rows, row])

    df_transmitter_serial = pd.DataFrame(
        {
            "measurementValue": [ds["instrument_tag"].attrs["serial_number"]],
            "measurementType": ["tag serial number"],
            "measurementTypeID": [
                "http://vocab.nerc.ac.uk/collection/MVB/current/MVB000189/"
            ],
            "measurementMethod": [""],
            "measurementUnit": [""],
        }
    )
    new_rows = pd.concat([new_rows, df_transmitter_serial], ignore_index=True)

    tag_manu = pd.DataFrame(
        {
            "measurementValue": [ds["instrument_tag"].attrs["manufacturer"]],
            "measurementType": ["tag manufacturer"],
            "measurementTypeID": [
                "http://vocab.nerc.ac.uk/collection/MVB/current/MVB000183/"
            ],
            "measurementMethod": [""],
            "measurementUnit": [""],
        }
    )
    new_rows = pd.concat([new_rows, tag_manu], ignore_index=True)

    tag_makemodel = pd.DataFrame(
        {
            "measurementValue": [ds["instrument_tag"].attrs["make_model"]],
            "measurementType": ["tag make and model"],
            "measurementTypeID": [
                "http://vocab.nerc.ac.uk/collection/MVB/current/MVB000185/"
            ],
            "measurementMethod": [""],
            "measurementUnit": [""],
        }
    )
    new_rows = pd.concat([new_rows, tag_makemodel], ignore_index=True)

    # sometimes we don't have attachment information
    if "attachment" in ds.attrs:
        attachment_location = pd.DataFrame(
            {
                "measurementValue": [ds.attrs["attachment"]],
                "measurementType": ["tag attachment location"],
                "measurementTypeID": [
                    "http://vocab.nerc.ac.uk/collection/MVB/current/MVB000395/"
                ],
                "measurementMethod": [""],
                "measurementUnit": [""],
            }
        )
        new_rows = pd.concat([new_rows, attachment_location], ignore_index=True)

    # --- Add eventID and occurrenceID. Assign the first occurrenceID and eventID from the dwc_df.
    new_rows["eventID"] = dwc_df["eventID"].iloc[0]
    new_rows["occurrenceID"] = dwc_df["occurrenceID"].iloc[0]

    # --- Select column order because this matters for meta.xml
    columns = [
        "eventID",
        "occurrenceID",
        "measurementValue",
        "measurementType",
        "measurementTypeID",
        "measurementMethod",
        "measurementUnit",
        "measurementUnitID",
    ]

    emof_df = new_rows[columns].copy()

    # drop any empty value rows.
    emof_df.dropna(axis=0, subset=["measurementValue"], inplace=True)

    if emof_df.empty:
        print("  no emof data found")
        return pd.DataFrame()  # Return an empty DataFrame if no observations are found
    else:
        emof_df.to_csv(output_csv.replace("occurrence", "emof"), index=False)
        print(f"  Created {len(emof_df)} emofs.")
        print(f"  Saved data to {output_csv.replace('occurrence', 'emof')}")
        return emof_df


def create_meta_xml(
    dwc_df: pd.DataFrame,
    emof_df: pd.DataFrame,
    event_df: pd.DataFrame,
    output_csv: str,
    cols: list,
):
    """
    Create meta.xml file for the Darwin Core dataset.

    Args:
        dwc_df (DataFrame): DataFrame containing Darwin Core occurrence data.
        emof_df (DataFrame): DataFrame containing eMoF data.
        event_df (DataFrame): DataFrame containing event data.
        output_csv (str): Path to the output CSV file.
        dir (str): Directory where the meta.xml will be saved.
        cols (list): List of occurrence columns to include in the meta.xml.
    """
    # Ensure the directory exists
    if not os.path.exists(output_csv):
        print(f"Missing directory: {output_csv}")

    # create and include the meta.xml and eml.xml
    # set the meta.xml paramaters by hand, using the format of the dataframes above
    meta_xml_vars = {}

    # when writing dwc occurrence file, we only save some columns
    dwc_df = dwc_df[cols].copy()

    meta_xml_vars["cols_list"] = dwc_df.columns.tolist()
    meta_xml_vars["occurrence_filename"] = output_csv
    meta_xml_vars["emof_cols_list"] = emof_df.columns.tolist()
    meta_xml_vars["emof_filename"] = output_csv.replace("occurrence", "emof")
    meta_xml_vars["event_cols_list"] = event_df.columns.tolist()
    meta_xml_vars["event_filename"] = output_csv.replace("occurrence", "event")

    # grab the template file for making meta.xml
    meta_template_file = codecs.open("templates/meta.xml.j2", "r", "UTF-8").read()
    meta_template = Template(meta_template_file)
    meta_result_string = meta_template.render(meta_xml_vars)
    dir = os.path.join(*output_csv.split("\\")[:-1])
    meta_file = f"{dir}/meta.xml"

    fh = codecs.open(meta_file, "wb+", "UTF-8")
    fh.write(meta_result_string)
    fh.close()
    meta_full_path = os.path.abspath(meta_file)
    print(f"  Meta XML has been written to '{meta_full_path}'.")


def package_dwc_zip(output_dir="data/dwc", zip_filename="data/dwc_package.zip"):
    """
    Packages all CSV and XML files in the specified output directory into a zip file.

    Args:
        output_dir (str): The directory containing the files to package.
        zip_filename (str): The path for the output zip file.
    """
    print(f"  Packaging Darwin Core files from '{output_dir}' into '{zip_filename}'...")
    with zipfile.ZipFile(zip_filename, "w", zipfile.ZIP_DEFLATED) as zipf:
        for root, _, files in os.walk(output_dir):
            for file in files:
                if file.endswith((".csv", ".xml")):
                    file_path = os.path.join(root, file)
                    arcname = os.path.relpath(file_path, output_dir)
                    zipf.write(file_path, arcname)
    print(f"  ✅ Packaged files into '{zip_filename}'")


def convert_to_dwc_individual(file_paths, output_dir="data/dwc"):
    """
    Converts a list of NetCDF files to individual Darwin Core Occurrence CSVs.

    An "occurrence" is the first detection of an animal at a specific
    location within a given hour.

    Args:
        file_paths (list): A list of paths to the .nc files.
        output_dir (str): The directory to save the individual CSV files.
    """
    print("\n--- 2. Starting Darwin Core Conversion (Individual Files) ---")
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        print(f"Created output directory: {output_dir}")

    processed_count = 0

    for nc_file in file_paths:
        base_filename = os.path.basename(nc_file)
        sub_dir = base_filename.split(".")[0]

        if not os.path.exists(f"{output_dir}/{sub_dir}"):
            os.makedirs(f"{output_dir}/{sub_dir}")
            print(f"Created output directory: {output_dir}/{sub_dir}")

        output_csv = os.path.join(
            output_dir, f"{sub_dir}/{os.path.splitext(base_filename)[0]}_occurrence.csv"
        )
        output_csv = os.path.normpath(output_csv)

        print(f"Processing {base_filename}...")

        try:
            with xr.open_dataset(nc_file, engine="netcdf4") as ds:
                df = ds.to_dataframe().reset_index()

                print(f"Found {len(df)} records.")

                # --- Data Cleaning and Preparation ---
                if "lat" not in df.columns or "lon" not in df.columns:
                    print(f"  Skipping {base_filename}: missing location data.")
                    continue

                df.dropna(subset=["lat", "lon", "time"], inplace=True)
                if df.empty:
                    print(f"  Skipping {base_filename}: no valid records.")
                    continue

                # --- Map to Darwin Core Occurrence Terms ---
                dwc_df, cols = create_dwc_occurrence(ds, output_csv, df_map)

                # Create and save eml
                create_eml(ds, df_map)

                # --- Event and eMoF (as needed) ---
                event_df = create_dwc_event(ds, dwc_df, output_csv)
                emof_df = create_dwc_emof(ds, dwc_df, output_csv)

                # --- Create meta.xml file ---
                create_meta_xml(dwc_df, emof_df, event_df, output_csv, cols)

                # --- Package into DwC-A ---
                output_dir_zip = f"data/dwc/{sub_dir}/"
                zip_filename = (
                    f"data/dwc/{sub_dir}/{base_filename.replace('.nc', '.zip')}"
                )
                package_dwc_zip(output_dir=output_dir_zip, zip_filename=zip_filename)

                processed_count += 1

        except Exception as e:
            print(f"  Could not process {base_filename}: {e}")

    print("\n--- 3. Conversion Complete ---")
    print(f"✅ Success! Processed {processed_count} files.")


def open_ipt_session(ipt_auth, ipt_url):
    """
    Begin a session with the target IPT
    Author: Jon Pye
    :param ipt_auth: Authentication details for the ipt, of the form {'email': 'email@mailserver.com', 'password':'cleartextPassword'}
    :param ipt_url: URL of the IPT we are authenticating with.
    :return: None
    """

    # relative path to IPT login form
    login_url = ipt_url + "login.do"

    s = requests.Session()  # open a session

    # retrieve the login form
    resp = s.get(login_url)

    # login forms generate a CSRF token that we have to persist in our response
    soup = BeautifulSoup(resp.text, "lxml")

    # Add it to our credentials dictionary
    ipt_auth["csrfToken"] = soup.find("input", {"name": "csrfToken"})["value"]
    login = s.post(login_url, data=ipt_auth)

    if login.status_code != 200:
        print("Login failed, status code {}".format(login.status_code))
        print(login.text)
        return None
    else:
        return s


def create_new_ipt_project(projname: str, filepath: str, ipt_url: str, ipt_session):
    """
    Create a new project on the given IPT using an existing DwC archive zip
    Author: Jon Pye
    :param projname: the project name as given by get_obis_shortname()
    :param filepath: payload resource filepath
    :param ipt_url: URL of the IPT to publish to
    :param ipt_session: authenticated requests session for the IPT
    :return: URL of the resource
    """

    path, filename = os.path.split(filepath)

    if not filename:  # if the filepath has no name in it
        print("no file specified in filepath, aborting")
        return None
    else:
        print(path, filename)
        print(filepath)

    # if there IS a file and it is not a valid DwC Archive, do we want to do anything here? The IPT runs its own checks...

    values = MultipartEncoder(
        fields={
            "create": "Create",  # hidden form fields with values
            "shortname": projname,
            "resourceType": "samplingevent",
            "__checkbox_importDwca": "true",
            "importDwca": "true",
            "file": (filename, open(filepath, "rb"), "application/x-zip-compressed"),
        }
    )
    create_dataset = ipt_session.post(
        ipt_url + "manage/create.do",
        data=values,
        headers={"Content-Type": values.content_type},
    )
    return create_dataset


def refresh_ipt_project_files(projname: str, filepath: str, ipt_url: str, ipt_session):
    """
    TODO: This should only push the .csv files, not the eml and meta

    Update data for a project on the given IPT using an existing DwC archive zip
    Author: Jon Pye
    :param projname: the project name as given by get_obis_shortname()
    :param filepath: payload resource filepath
    :param ipt_url: URL of the IPT to publish to
    :param ipt_session: authenticated requests session for the IPT
    :return: URL of the resource
    """

    path, filename = os.path.split(filepath)

    if not filename:  # if the filepath has no name in it
        print("no file specified in filepath, aborting")
        return None

    values = MultipartEncoder(
        fields={
            "add": "Add",
            "r": projname,
            "sourceType": "source-file",
            "validate": "false",
            "file": (filename, open(filepath, "rb"), "application/x-zip-compressed"),
        }
    )

    update_dataset = ipt_session.post(
        ipt_url + "manage/addsource.do",
        data=values,
        headers={"Content-Type": values.content_type},
    )
    if update_dataset.status_code == 200:
        # Handle the Are you Sure popup.
        print("Publication successful")
        return update_dataset
    else:
        print("publication error, check landing page output")
        return update_dataset


def refresh_ipt_project_metadata(
    projname: str, filepath: str, ipt_url: str, ipt_session
):
    """
    Update metadata for a project on the given IPT using an existing eml.xml file
    Author: Jon Pye
    :param projname: the project name as given by get_obis_shortname()
    :param filepath: payload resource filepath
    :param ipt_url: URL of the IPT to publish to
    :param ipt_session: authenticated requests session for the IPT
    :return: URL of the resource
    """

    path, filename = os.path.split(filepath)

    if not filename:  # if the filepath has no name in it
        print("no file specified in filepath, aborting")
        return None

    values = MultipartEncoder(
        fields={
            "emlReplace": "Replace",
            "r": projname,
            "sourceType": "source-file",
            "validateEml": "true",
            "__checkbox_validateEml": "true",
            "emlFile": (filename, open(filepath, "rb"), "application/xml"),
        }
    )

    update_metadata = ipt_session.post(
        ipt_url + "manage/replace-eml.do",
        data=values,
        headers={"Content-Type": values.content_type},
    )
    return update_metadata


def change_publishing_org_ipt_project(
    projname: str, ipt_url: str, ipt_session, new_publishing_org_name: str
):
    """
    Change the publishing organization in the given IPT project
    Author: Mathew Biddle
    :param projname: the project name as given by get_obis_shortname()
    :param ipt_url: URL of the IPT to publish to
    :param ipt_session: authenticated requests session for the IPT
    :param new_publishing_org: the new publishing organisation to set for this project. See publishingOrganizationKey in the IPT source.

    :return: URL of the resource
    """
    pub_orgs = {
        "NOAA Integrated Ocean Observing System": "1d38bb22-cbea-4845-8b0c-f62551076080",
        "No organization": "625a5522-1886-4998-be46-52c66dd566c9",
        "SCAR - AntOBIS": "104e9c96-791b-4f14-978c-f581cb214912",
        "The Marine Genome Project": "aa0b26e8-779c-4645-a569-5f39fa85d528",
        "USFWS-AK": "530fda11-7af7-4447-9649-0f9fc22e6156",
        "United States Fish and Wildlife Service": "f8dbeca7-3131-41ab-872f-bfad71041f3f",
        "United States Geological Survey": "c3ad790a-d426-4ac1-8e32-da61f81f0117",
    }

    if new_publishing_org_name not in pub_orgs:
        print(
            f"Publishing organization '{new_publishing_org_name}' not recognised as one of {pub_orgs.keys()}. Please check the name and try again."
        )
        return None

    pub_params = {
        "r": projname,  # resource = dataset name
        "publishingOrganizationKey": pub_orgs[new_publishing_org_name],
    }

    contents = ipt_session.post(
        ipt_url + "manage/resource-changePublishingOrganization.do", data=pub_params
    )
    return contents


def make_public_ipt_project(projname: str, ipt_url: str, ipt_session):
    """
    Update metadata for a project on the given IPT
    Author: Jon Pye
    :param projname: the project name as given by get_obis_shortname()
    :param ipt_url: URL of the IPT to publish to
    :param ipt_session: authenticated requests session for the IPT
    :return: URL of the resource
    """
    pub_params = {
        "r": projname,  # resource = dataset name
        "makePrivate": "Public",
    }

    contents = ipt_session.post(
        ipt_url + "manage/resource-makePublic.do", data=pub_params
    )
    return contents


def publish_ipt_project(
    projname: str, ipt_url: str, ipt_session, publishing_notes: str = ""
):
    """
    Update metadata for a project on the given IPT
    Author: Jon Pye
    :param projname: the project name as given by get_obis_shortname()
    :param ipt_url: URL of the IPT to publish to
    :param ipt_session: authenticated requests session for the IPT
    :param publishing_notes: optional message to publish this version with
    :return: URL of the resource
    """

    pub_params = {
        "r": projname,  # resource = dataset name
        "autopublish": "",
        "currPubMode": "AUTO_PUBLISH_OFF",
        "pubMode": "",
        "currPubFreq": "",
        "pubFreq": "",
        "publish": "Publish",
        "summary": publishing_notes,
    }
    contents = ipt_session.post(ipt_url + "manage/publish.do", data=pub_params)
    return contents


def register_ipt_project(projname: str, ipt_url: str, ipt_session):
    """
    Update Register the given IPT project with GBIF
    Author: Mathew Biddle
    :param projname: the project name as given by get_obis_shortname()
    :param ipt_url: URL of the IPT to publish to
    :param ipt_session: authenticated requests session for the IPT

    Need to do the following in the dialog-confirm.
       * check checkbox-confirm
       * select yes-button

    :return: URL of the resource
    """

    pub_params = {
        "r": projname,  # resource = dataset name
        "checkbox-confirm": "true",  # checkbox-confirm
        "yes-button": "Yes",
    }

    contents = ipt_session.post(
        ipt_url + "manage/resource-registerResource.do", data=pub_params
    )
    return contents


def check_if_project_exists(projname: str, ipt_url: str, ipt_session):
    """
    Test if a project exists on the IPT already
    Author: Jon Pye
    :param projname: the project name as given by get_obis_shortname()
    :param ipt_url: URL of the IPT to check for this publication
    :param ipt_session: authenticated requests session for the IPT
    :return: True if the project already exists on the IPT in question
    """

    checkUrl = "{ipt_url}ipt/resource?r={projname}".format(
        ipt_url=ipt_url, projname=projname
    )

    contents = ipt_session.post(checkUrl)

    # if it's not found, the IPT returns a 404
    if contents.status_code == 404:
        print("No existing repository by this name: '{}'".format(projname))
        return False
    elif contents.status_code == 200:
        print("Found existing project by name: '{}'".format(projname))
        return True
