# 01 - P450 BM3 Wild-Type vs Full Library

## Aim

Screen entire FDA library against P450 BM3 Wild-Type.

Use more than 4 concentrations for each compound. 

Not all of the data will be useable because some of the compounds interfere with measurement, but these edge cases can be caught for process development.

## Overview

```
.
├── echo
│   ├── 2022-04-21
│   │   ├── E5XX-1564_Survey_Source[2](UnknownBarCode).csv
│   │   ├── E5XX-1564_Survey_Source[3](UnknownBarCode).csv
│   │   ├── E5XX-1564_Survey_Source[4](UnknownBarCode).csv
│   │   ├── E5XX-1564_Transfer_0.csv
│   │   ├── E5XX-1564_Transfer_1650538543.csv
│   │   └── E5XX-1564_Transfer_1650538543_Exceptions.csv
│   ├── 2022-04-21-exp01-test.epr
│   └── picklist
│       └── 01-picklist.csv
├── nb
│   ├── design-01.ipynb
│   └── layouts.csv -> ../../../data/lib/layouts.csv
├── platereader
├── readme.md
└── uv-vis

6 directories, 11 files
```

## Notes

### 2022-04-21

8:30 - 10:30ish? : Filling echo source plates from FDA lib. 
11µl.
I broke source plate `[0 .. 4]` in the centrifuge - 3,500(?)(rpm/rcf?) - two stacked in each bucket, in one the two plates cracked on the flange.
Lucky one was a balance.
Transferred contents of `[0 .. 4]` to new plate.

Only did one set, will do the other tomorrow with the same source plates.

11:00 - 13:00(?) : Echo Dispensing.

- `echo/picklist/01-picklist.csv` -> `echo/2022-04-21-exp01.test.epr`
- `echo/2022-04-21` : logs

### 2022-04-22

- Echo booked morning
- Platereader booked afternoon
