"""
This module reads the molecules from text files and returns them as a dictionary.
"""
import os
import ast

class Get_smiles:
    def __init__(self):
        self.here = os.path.dirname(__file__)

class Alcohol(Get_smiles):
    def get(self, txt_path = "Alcohol.txt"):
        # Load the alcohols from the text file
        fullpath = os.path.join(self.here, txt_path)
        with open(fullpath, 'r') as f:
            txt = f.read()
            # parses the literal dict safely
            molecules = ast.literal_eval(txt)
        return molecules

class Carboxylic_acid(Get_smiles):
    def get(self, txt_path = "Carboxylic_acid.txt"):
        # Load the carboxylic acids from the text file
        fullpath = os.path.join(self.here, txt_path)
        with open(fullpath, 'r') as f:
            txt = f.read()
            # parses the literal dict safely
            molecules = ast.literal_eval(txt)
        return molecules
    
class Ester(Get_smiles):
    def get(self, txt_path = "Ester.txt"):
        # Load the esters from the text file
        fullpath = os.path.join(self.here, txt_path)
        with open(fullpath, 'r') as f:
            txt = f.read()
            # parses the literal dict safely
            molecules = ast.literal_eval(txt)
        return molecules

class AcylHalide(Get_smiles):
    def get(self, txt_path = "AcylHalide.txt"):
        # Load the acyl halides from the text file
        fullpath = os.path.join(self.here, txt_path)
        with open(fullpath, 'r') as f:
            txt = f.read()
            # parses the literal dict safely
            molecules = ast.literal_eval(txt)
        return molecules