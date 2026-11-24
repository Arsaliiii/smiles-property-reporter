from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors, Lipinski

mol = Chem.MolFromSmiles("CC(=O)Oc1ccccc1C(=O)O")   # aspirin