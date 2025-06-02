from rdkit import Chem
from rdkit.Chem import Draw, AllChem
from rdkit.Chem import rdChemReactions
import ast
import os

class ChemReaction():
    """Create an object that can perform chemical reactions and
      contain a dictionary of reaction SMARTS and catalysts from a text file."""
    def __init__(self, txt_path = "reactions_smarts.txt"): #text file contain the dictionary of reaction SMARTS
        # Load the reaction SMARTS from the text file
        here = os.path.dirname(__file__)
        fullpath = os.path.join(here, txt_path)
        with open(fullpath, 'r') as f:
            txt = f.read()
            # parses the literal dict safely
            self.reaction_smarts = ast.literal_eval(txt)

    def reacting_catalyst_products(self, reaction, reactants):
        # Create the reaction from the SMARTS string
        # Check if the reaction is in the dictionary
        if reaction not in self.reaction_smarts:
            print(f"Reaction '{reaction}' not found in the dictionary.")
            return None
        else:
            rsmarts = self.reaction_smarts[reaction]["smarts"] # get the reaction SMARTS from the dictionary
        rxn = AllChem.ReactionFromSmarts(rsmarts)

        # Perform the reaction
        products = rxn.RunReactants(reactants)
        # Get the catalyst from the dictionary
        catalyst = self.reaction_smarts[reaction].get("catalyst") # get the catalyst from the dictionary
        # Return the products (if any)
        if products:
            if catalyst:
                return catalyst , products  # Return the catalysts and products of the reactions
            else:
                return None, products
        else:
            print("No products found for the given reactants.")
            return None, None

class Make_reaction():
    """Create the reaction and Visualize a chemical reaction"""
    def __init__(self, reaction, reactants): #create the reaction
        self.reaction = reaction
        self.reactants = reactants
        self.catalyst, self.products = ChemReaction().reacting_catalyst_products(reaction, reactants)
        
    def visualizing(self):
        """Visualize the reaction"""
        if self.products:
            # Draw the reactants and products
            rxn = rdChemReactions.ChemicalReaction()
            for reactant in self.reactants:
                rxn.AddReactantTemplate(reactant)
            for product_set in self.products:
                for product in product_set:
                    rxn.AddProductTemplate(product)
            # Draw catalyst
            if self.catalyst:
                catalyst_mol = Chem.MolFromSmiles(self.catalyst)
                rxn.AddAgentTemplate(catalyst_mol)
            img = Draw.ReactionToImage(rxn, subImgSize=(300, 300))
            img.show()
        else:
            print("No products found for the given reactants.")

    def img_get(self):
        """Visualize the reaction"""
        if self.products:
            # Draw the reactants and products
            rxn = rdChemReactions.ChemicalReaction()
            for reactant in self.reactants:
                rxn.AddReactantTemplate(reactant)
            for product_set in self.products:
                for product in product_set:
                    rxn.AddProductTemplate(product)
            # Draw catalyst
            if self.catalyst:
                catalyst_mol = Chem.MolFromSmiles(self.catalyst)
                rxn.AddAgentTemplate(catalyst_mol)
            img = Draw.ReactionToImage(rxn, subImgSize=(200, 200))
            return img
        else:
            print("No products found for the given reactants.")
            return None


# reactants = [
#     Chem.MolFromSmiles("CC(=O)OC")
#     , Chem.MolFromSmiles("O")
# ]
# rxn = Make_reaction("Ester_hydrolysis", reactants)
# rxn.visualizing()