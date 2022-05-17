# Data Analysis

## Summary

This page describes how screening data was processed in this project.

## Aim

- Identify binding interactions between P450s and compounds from UV-Visible light absorbance profiles between 220 and 800 nm from plate assay.
- Compensate for compound light absorbance in the sensitive area (390-420 nm) using control plates.
- Output a table of enzymes, compounds, a qualifitative or quantitative measure of binding and an indication of data quality.

## Procedure

During each screening experiment, data collected were:

- Plate data, which was exported to `.csv` files on the platereader host machine.
- *Echo* logs and exception reports, also `.csv` files.
- UV-Visible spectrometer readings taken when thawing, diluting and dispensing the target protein, exported to `.csv`s.

The key tasks are:

1. Use *Echo* picklists and exception reports to map compounds and their volumes to plate wells.
2. Match traces in platereader files to compounds and volumes.
3. Process and normalize traces, subtract control traces to correct for compound absorbance.
4. Detect and filter anomalies.
5. Quantify P450 absorbance profile response to compound concentration.
6. Get binding metrics $K_d$ and $V_{max}$ for each compound-P450 pair.

### Configuration and Data Processing 

Each experiment had a directory that looked like this:

```
02.0/                    
├── config.yml           # maps experiment directory for analysis
├── echo                 
│   ├── echo_logs        # surveys, exception reports etc
│   └── picklist         # 
├── img                  # output specs of 
├── nb                   # jupyter notebooks used for lab notes
├── platereader          # raw platereader data
└── uv-vis               # raw UV-Visible light spectroscopy data
```


Configuration files `config.yml` were generated in directory partly using the script `screening-fist/sxfst/scripts/config.py`, partly with human intervention.
`config.yml` contains experimental metadata that can be used to map compounds and their end volume in each well.
`yml` format was chosen because it's human friendly data structure and is easily handled with the `python` `yaml` library.
Each file looked like this:

``` yaml

echo:
  picklist:
  - ../../lab/02.0/echo/picklist/01-picklist.csv
  protocol:
  - ../../lab/02.0/echo/2022-04-23-exp02.epr
  surveys:
  - ../../lab/02.0/echo/2022-04-22/E5XX-1564_Survey_Source[4](UnknownBarCode).csv
  - ../../lab/02.0/echo/2022-04-22/E5XX-1564_Survey_Source[3](UnknownBarCode).csv
  - ../../lab/02.0/echo/2022-04-22/E5XX-1564_Survey_Source[2](UnknownBarCode).csv
  transfers:
  - ../../lab/02.0/echo/2022-04-22/E5XX-1564_Transfer_1650618552_Exceptions.csv
  - ../../lab/02.0/echo/2022-04-22/E5XX-1564_Transfer_1650618552.csv
nb:
- ../../lab/02.0/nb/design-01.ipynb
- ../../lab/02.0/nb/lab.ipynb
platereader:
  plate_1:
    control:
      date: 23/04/2022
      id1: null
      id2: null
      machine: BMG CLARIOstar
      path: ../../lab/02.0/platereader/TRno2900.CSV
      test_run_no: '2900'
      time: '16:30:14'
      user: USER
    test:
      date: 24/04/2022
      id1: null
      id2: null
      machine: BMG CLARIOstar
      path: ../../lab/02.0/platereader/24042022,0719350719.CSV
      test_run_no: '2930'
      time: 07:19:35
      user: USER
```

A script - `screening-fist/sxfst/scripts/data_proc.py` uses `config.yml` to map compounds and their volumes to absorbance traces from the platereader data.
It outputs plots of the normalized traces.
It may be more useful to have it output one big old table of compounds, traces etc.
Normalizing etc can come later

### Signal Processing and Anomaly Detection 



### Response Qualification and Quantification 
