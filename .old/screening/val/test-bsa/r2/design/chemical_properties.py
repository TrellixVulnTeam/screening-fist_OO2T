import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from rdkit.Chem import AllChem as Chem 
from rdkit.Chem import Draw

df = pd.read_csv('selection.csv', index_col=0).dropna()
mols = [Chem.MolFromSmiles(i) for i in df['SMILES']]
### pic
im=Draw.MolsToGridImage(mols, legends = df['Item Name'].tolist(), molsPerRow=10)
im.save('compound_selection_structures.png')
### hydrogens for property calcs
mols = [Chem.AddHs(i) for i in mols]


