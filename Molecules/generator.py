# This module is no longer used.
from rdkit import Chem

class Basic_molecule_gen(): 
    "generate a basic (nornmal) molecule from number of carbon and functional groups"
    def __init__(self, num_carbons, functional_groups):
        self.num_carbons = num_carbons
        self.functional_groups = functional_groups
        self.smiles = self.num_carbons * "C" + "".join(self.functional_groups)
        self.mol = Chem.MolFromSmiles(self.smiles)
