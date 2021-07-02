# screening library files
## csvs
- [original lib files](selleckchem-plate.xlsx)
- [layouts in plates](layouts.csv) - cleaned from original files
- [herbicide-like](herbicide-like.csv) - filtered subset of [layouts.csv](layouts.csv) based on herbicide-likeness rules

## python programs
- [readSheet.py](readSheet.py) - parses and cleans [selleckchem-plate.xlsx](selleckchem-plate.xlsx) -> [layouts.csv](layouts.csv)
- [herbicide-like.py](herbicide-like.py) filters [layouts.csv](layouts.csv) based on herbicide-likeness rules (ref) -> [herbicide-like.csv](herbicide-like.csv) - requires rdkit, outputs csv with ```;``` as sep, because the smiles strings contain commas
- [select-compounds.py](select-compounds.py) - pick n compounds from [herbicide-like.csv](herbicide-like.csv) using a maxmin picker diversity filter, outputs to [herbicide-like-selection.csv](herbicide-like-selection.csv) - requires rdkit
