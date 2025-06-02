import random

def randomize_reactants(reactants):
    randomized_reactants = random.choice(list(reactants.items()))
    return randomized_reactants
