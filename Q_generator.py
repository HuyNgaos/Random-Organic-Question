"""
This module generates questions and answers for the flash card app. With randomly selected
molecules.
"""
from Utls import *
from Molecules import ReadMolecules as RM
from Reactions import chemical_reactions
import random
from rdkit import Chem

class Question_gen:
    def read_molecules(self, functional_group):
        match functional_group: # read the molecules from the text file matching the functional group
            case "Alcohol":
                molecules = RM.Alcohol().get()
            case "Carboxylic_acid":
                molecules = RM.Carboxylic_acid().get()
            case "Ester":
                molecules = RM.Ester().get()
            case _:
                raise ValueError("Invalid functional group")
        return molecules
    
class Esterfication(Question_gen):
    def __init__(self): # Read the correct molecules from the text file and set the reaction type
        # super().__init__(question_type)
        self.Carboxylic_acid = self.read_molecules("Carboxylic_acid")
        self.Alcohol = self.read_molecules("Alcohol")
        self.reaction = "Esterification"

    def generate_question(self):
        # Generate a random carboxylic acid and alcohol
        carboxylic_acid_name, carboxylic_acid_smiles = random.choice(list(self.Carboxylic_acid.items()))
        alcohol_name, alcohol_smiles = random.choice(list(self.Alcohol.items()))
        # Generate the the reactants list and reaction
        reactants = [
            Chem.MolFromSmiles(carboxylic_acid_smiles),
            Chem.MolFromSmiles(alcohol_smiles)
        ]
        rxn = chemical_reactions.Make_reaction(self.reaction, reactants)
        # Generate the question
        question = f"What is the product of the reaction between {carboxylic_acid_name} and {alcohol_name}?"
        # Return the question and the reaction (including the product)
        return question, rxn
    
class Ester_hydrolysis(Question_gen):
    def __init__(self): # Read the correct molecules from the text file and set the reaction type
        # super().__init__(question_type)
        self.Ester = self.read_molecules("Ester")
        self.reaction = "Ester_hydrolysis"

    def generate_question(self):
        # Generate a random ester and alcohol
        ester_name, ester_smiles = random.choice(list(self.Ester.items()))
        # Generate the the reactants list and reaction
        reactants = [
            Chem.MolFromSmiles(ester_smiles),
            Chem.MolFromSmiles("O")
        ]
        rxn = chemical_reactions.Make_reaction(self.reaction, reactants)
        # Generate the question 
        question = f"What are the products of the hydrolysis {ester_name}?"
        # Return the question and the reaction (including the product)
        return question, rxn
    
# def test():
#     question_gen = Ester_hydrolysis()
#     question, answer = question_gen.generate_question()
#     print(question)
#     print(Chem.MolToSmiles(answer.products[0][0]))
#     answer.visualizing()

# test()

