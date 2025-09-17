def create_dwc_occurrence(
    ds: xr.Dataset, output_csv: str, ncei_accession_mapping: pd.DataFrame
):
    """Create a Darwin Core Occurrence CSV from an xarray Dataset."""
    source_file = os.path.basename(ds.encoding.get("source"))
    # bail if we can't find the file in the mapping table.
    if source_file not in ncei_accession_mapping["file_name"].values:
        raise KeyError(f"File {source_file} not found in NCEI Accession mapping table.")

    filename = os.path.splitext(source_file)[
        0
    ]  # "ioos_atn_{ds.ptt_id}_{start_date}_{end_date}""

    file_map_entry = ncei_accession_mapping[
        ncei_accession_mapping["file_name"] == source_file
    ].iloc[0]

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


def convert_to_dwc_individual(fname):
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
