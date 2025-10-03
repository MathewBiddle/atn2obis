import xarray as xr
import pytest
from siphon.catalog import TDSCatalog

from atn2obis.fetch_data import _get_name, _opendap_urls, _nested_catalogs


_catalog_url = "https://www.ncei.noaa.gov/thredds-ocean/catalog/ioos/atn/catalog.xml"


@pytest.fixture
def nested_catalogs():
    catalog = TDSCatalog(catalog_url=_catalog_url)
    return list(_nested_catalogs(catalog))


def test__get_name():
    assert "ioos/atn" == _get_name(_catalog_url)


def test__nested_catalogs(nested_catalogs):
    assert isinstance(nested_catalogs, list)
    assert isinstance(nested_catalogs[0], TDSCatalog)


def test__opendap_urls(nested_catalogs):
    cat = nested_catalogs[0]
    urls = _opendap_urls(cat)
    assert isinstance(urls, list)
    ds = xr.open_dataset(urls[0], drop_variables="time")
    assert isinstance(ds, xr.Dataset)
