import re
import shapely
import os  # replace with pathlib
from pathlib import Path
import codecs  # should not be necessary in py>=3
from jinja2 import Template
import xarray as xr
import pandas as pd


from ncei_mapping import ncei_accession_mapping as df_map


def create_dwc_occurrence(ds: xr.Dataset, output_csv: str, df_map: pd.DataFrame):
    """Create a Darwin Core Occurrence CSV from an xarray Dataset."""
    fname = ds.encoding.get("source")
    if fname is not None:
        source_file = Path(fname).name
    else:
        msg = f"Cannot find file name in {ds.encoding=}."
        raise ValueError(msg)
    file_map_entry = df_map[df_map["file_name"] == source_file].squeeze()

    # bail if we can't find the file in the mapping table.
    if file_map_entry.empty:
        raise KeyError(f"File {source_file} not found in NCEI Accession mapping table.")

    filename = Path(source_file).stem  # "ioos_atn_{ds.ptt_id}_{start_date}_{end_date}""
    dwc_df = pd.DataFrame()

    date = ds["time"].dt.strftime("%Y-%m-%dT%H:%M:%SZ")
    occurrenceID = f"ioos_atn_{ds['ptt_id']}_{date}_{ds['z']}_{ds['animal_common_name'].replace(' ', '_')}"
    organismID = f"{ds['platform_id']}_{ds['animal_common_name'].replace(' ', '_')}"
    associatedReferences = f"https://doi.org/10.25921/wp4e-ph20; https://www.ncei.noaa.gov/archive/accession/{file_map_entry['accession']}; {file_map_entry['related_data_url']}"

    dwc_df["occurrenceID"] = occurrenceID
    dwc_df["eventID"] = filename
    dwc_df["organismID"] = organismID
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
    dwc_df["associatedReferences"] = associatedReferences
    # FIXME: Should these be min/max?
    dwc_df["minimumDepthInMeters"] = ds["z"].values.tolist()
    dwc_df["maximumDepthInMeters"] = ds["z"].values.tolist()
    dwc_df["bibliographicCitation"] = ds["citation"]

    # set basisOfRecord...
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
        f"This is a representative occurrence from a full deployment. For the complete dataset please see https://www.ncei.noaa.gov/archive/accession/{file_map_entry['accession']}."
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


## Create the DwC Event file from the occurrences.
def create_dwc_event(ds: xr.Dataset, dwc_df: pd.DataFrame, output_csv: str):
    # create parent event that is a summary of dwc_df
    event_df = pd.DataFrame()
    event_df["eventID"] = dwc_df["eventID"].unique()
    event_df["eventDate"] = dwc_df["eventDate"].min() + "/" + dwc_df["eventDate"].max()

    ## Convex hull summary of the points
    points = list(zip(dwc_df["decimalLongitude"], dwc_df["decimalLatitude"]))
    event_df["footprintWKT"] = shapely.convex_hull(shapely.LineString(points))

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


## Create the DwC Extended Measurement or Fact (eMoF) file.
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
    output_csv: Path,
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
    if not output_csv.exists():
        # FIXME: Should we raise something here?
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
    # FIXME: Use pathlib and avoid Windows doulble \
    # FIXME: Do not clobber built-in `dir`.
    dir = os.path.join(*output_csv.split("\\")[:-1])
    meta_file = f"{dir}/meta.xml"

    fh = codecs.open(meta_file, "wb+", "UTF-8")
    fh.write(meta_result_string)
    fh.close()
    meta_full_path = os.path.abspath(meta_file)
    print(f"  Meta XML has been written to '{meta_full_path}'.")


def convert_to_dwc_individual(fname, output_csv: Path):
    with xr.open_dataset(fname, engine="netcdf4") as ds:
        df = ds.to_dataframe().reset_index()

        print(f"Found {len(df)} records.")

        # --- Data Cleaning and Preparation ---
        if "lat" not in df.columns or "lon" not in df.columns:
            msg = f"Skipping {fname}: missing location data."
            raise ValueError(msg)

        df.dropna(subset=["lat", "lon", "time"], inplace=True)
        if df.empty:
            msg = f"Skipping {fname}: no valid records."
            raise ValueError(msg)

        # --- Map to Darwin Core Occurrence Terms ---
        dwc_df, cols = create_dwc_occurrence(ds, output_csv, df_map)

        # Create and save eml
        create_eml(ds, df_map)

        # --- Event and eMoF (as needed) ---
        event_df = create_dwc_event(ds, dwc_df, output_csv)
        emof_df = create_dwc_emof(ds, dwc_df, output_csv)

        # --- Create meta.xml file ---
        create_meta_xml(dwc_df, emof_df, event_df, output_csv, cols)
