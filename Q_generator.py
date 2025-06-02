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
            case "AcylHalide":
                molecules = RM.AcylHalide().get()
            case _:
                raise ValueError("Invalid functional group")
        return molecules

"""
The subclasses of Question_gen generate questions and answers for the flash card app
for each type of reaction.
"""
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
    
class AcylHalideEsterification(Question_gen):
    def __init__(self): # Read the correct molecules from the text file and set the reaction type
        # super().__init__(question_type)
        self.AcylHalide = self.read_molecules("AcylHalide")
        self.Alcohol = self.read_molecules("Alcohol")
        self.reaction = "AcylHalideEsterification"

    def generate_question(self):
        # Generate a random acyl halide and alcohol
        acyl_halide_name, acyl_halide_smiles = random.choice(list(self.AcylHalide.items()))
        alcohol_name, alcohol_smiles = random.choice(list(self.Alcohol.items()))
        # Generate the the reactants list and reaction
        reactants = [
            Chem.MolFromSmiles(acyl_halide_smiles),
            Chem.MolFromSmiles(alcohol_smiles)
        ]
        rxn = chemical_reactions.Make_reaction(self.reaction, reactants)
        # Generate the question
        question = f"What is the product of the reaction between {acyl_halide_name} and {alcohol_name}?"
        # Return the question and the reaction (including the product)
        return question, rxn
    
# class Wittig_reaction(Question_gen):
#     def __init__(self): # Read the correct molecules from the text file and set the reaction type
#         self.AlkylHalide = self.read_molecules("AlkylHalide")
#         self.Ketone = self.read_molecules("Ketone")
#         self.Aldehyde = self.read_molecules("Aldehyde")
#         self.carbonyl = 
# def test():
#     question_gen = Esterfication()
#     question, answer = question_gen.generate_question()
#     print(question)
#     print(Chem.MolToSmiles(answer.products[0][0]))
#     answer.visualizing()

# test()

