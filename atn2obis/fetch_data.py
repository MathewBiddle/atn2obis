from siphon.catalog import TDSCatalog
from urllib.parse import urljoin


def _get_name(catalog_url):
    return catalog_url.split("catalog")[-2].strip("/")


def _opendap_urls(cat):
    return [value.access_urls.get("opendap") for value in cat.datasets.values()]


def _nested_catalogs(cat):
    # Reached end with datasets.
    if not cat.catalog_refs and cat.datasets:
        yield cat
    # Keep navigating the refs
    # FIXME: If there are datasets and refs, we will skip the dataset.
    # I'm not sure if that is common in a THREDDS server.
    if cat.catalog_refs:
        for catalog_ref in cat.catalog_refs:
            ref = urljoin(cat.catalog_url, f"{catalog_ref}/catalog.xml")
            new_cat = TDSCatalog(catalog_url=ref)
            yield from _nested_catalogs(new_cat)


def fetch_atn_catalog():
    base_catalog = TDSCatalog(
        catalog_url="https://www.ncei.noaa.gov/thredds-ocean/catalog/ioos/atn/catalog.xml"
    )
    nested_catalogs = _nested_catalogs(base_catalog)
    return {
        _get_name(nested_catalog.catalog_url): _opendap_urls(nested_catalog)
        for nested_catalog in nested_catalogs
    }
