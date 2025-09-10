# ATN to OBIS

## Goal: 
Create an open source package to submit summarized ATN satellite telemetry data from NCEI to OBIS-USA for inclusion in OBIS and GBIF: 
In collaboration with OBIS-USA, NCEI and IOOS, establish a pathway for ATN satellite telemetry data (sourced from NCEI) to be summarized and 
submitted to OBIS-USA Integrated Publishing Toolkit (IPT) to be shared with OBIS, GBIF, and archived at NCEI. 

The netCDF files at NCEI follow the ATN Satellite Telemetry Specification V1.0 as documented at https://ioos.github.io/ioos-atn-data/atn-sat-telem-specification-v1-0.html.

## Step-by-step:
1. Pulls source netCDF data directly from NCEI.
2. Ingests netCDF data and translates to DarwinCore format
3. Creates applicable eml metadata.
4. Creates applicable file column mapping xml file.
5. Zips everything together into a compliant DwC Archive package.
6. Automatically publish DwC Archive package to OBIS-USA IPT.
7. Provide capability to version control updates to the datasets on the OBIS-USA IPT by: updating of eml metadata, updating of data files, and republishing.

