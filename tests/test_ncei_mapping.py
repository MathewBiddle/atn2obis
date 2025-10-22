import pandas as pd
import pytest

from atn2obis.ncei_mapping import (
    _openurl_with_retry,
    get_ncei_accession_mapping,
)


_url = "https://www.ncei.noaa.gov/access/metadata/landing-page/bin/iso?id=gov.noaa.nodc:IOOS-ATN-STP;view=xml;responseType=text/xml"


def test__openurl_with_retry():
    r = _openurl_with_retry(_url)
    assert r.status == 200


@pytest.mark.parametrize("refresh", [True, False], ids=["fetch", "use_local"])
def test_get_ncei_accession_mapping(refresh):
    df = get_ncei_accession_mapping(refresh=refresh)
    assert isinstance(df, pd.DataFrame)
    assert not df.empty
    assert all(
        [
            col in df.columns
            for col in (
                "title",
                "accession",
                "start_date",
                "end_date",
                "related_data_url",
                "related_data_citation",
                "arc",
                "xml",
                "file_name",
                "ptt_id",
            )
        ]
    )
