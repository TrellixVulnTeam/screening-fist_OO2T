# 00 - pilot

## Aim

- Test design
- Test lab protocol
- Test data processing
- Test Dataset construction

## Requirements

- Medium-sized experiment array
    - Edge cases 
    
## Plan

Full lib vs 1 mutant.
Dump screen.
Old BM3 A82F/F87V full-length surplus.
Can clear active site with NADPH.
Metadata to database.

## Protocol

- [x] Make echo scripts
- [x] Lab (duration: 6h)

    1. Book echo, platereader
    2. Make 1 `Src` plate
    3. Dispense 2 `Dest` plates (control & test)
    4. Make buffer: `{'kpi':{'KPi-mM':100, 'pH':7.0, 'filtered':true, 'degassed':true}`
        - forgot to filter
    5. Prep `{'bm3-dm-fl':{'conc-um':10, 'buffer':'kpi'}` - specs 
        - used --- instead
    5. Dest plates: `{'control':{'compoud':{'vol-ul':1.5}, 'kpi':{'vol-ul':30}}` 
    `{'test':{'compoud':{'vol-ul':1.5}, 'kpi':{'vol-ul':30}}}` 
    6. Read `{'platereader':{'temp-c':25, 'wavelengths':'200-800'}}`
- [x] Data to bucket
- [ ] Bucket to db
- [ ] Analyze 
    - exploratory :hourglass:
    - automation 
    - end point: pytorch dataset

## Summary
- Used BM3 WT Heme & Compound racks 1-4
- data @ linode bucket `james-engleback/00-data.tar.gz`
- 6h lab time
- [Design](00-design-lab-notes-01.ipynb)
- [Lab notes](00-labNotes.ipynb)
