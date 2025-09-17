import functools
from owslib.iso import namespaces
import stamina
from urllib.request import urlopen
import urllib.error
import xml.etree.ElementTree as ET
import pandas as pd
import html
import re
from pathlib import Path


@stamina.retry(on=urllib.error.HTTPError, attempts=3)
def _openurl_with_retry(url):
    """Thin wrapper around urlopen adding stamina."""
    return urlopen(url)


@functools.lru_cache(maxsize=128)
def _get_ncei_accession_mapping():
    """
    Scrapes NCEI for ATN accession numbers and associated file metadata.

    Returns:
        pd.DataFrame: A DataFrame mapping accession numbers to file names,
                      download URLs, and other metadata.
    """
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
            # FIXME: Catching all exceptions will hide errors.
            # NB: That do we expect to fail here? `HTTPError`?, Parsing the xml?
            # The way this is we will catch import errors, other lib erros,
            # and move forward without knowing if we got it right
            print(f"Could not process accession {acc}: {e}")

    df_map["ptt_id"] = df_map["title"].str.extract(r".*ptt ([0-9]{3,7}) .*")
    return df_map



def get_ncei_accession_mapping(refresh=False):
    if refresh:
        ncei_accession_mapping = _get_ncei_accession_mapping()
    else:
        path = Path(__file__).absolute().parent
        ncei_accession_mapping = pd.read_csv(path.joinpath("ncei_accession_mapping.csv"))
    return ncei_accession_mapping
