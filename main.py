import Q_generator as Qg
from rdkit import Chem

def main():
    introduction = """    Welcome to the Flash Card App!
This app is still in development.
Please check back later for updates."""
    print(introduction)
    available_reactions = ["Esterfication", "Ester_hydrolysis"]
    user_input = input("Press type the reation you wish to review:\n" + str(available_reactions) + "\nEnter: ")
    reaction_string = "Qg." + user_input + "().generate_question()"
    question, answer = eval(reaction_string)
    print("\n" + question)
    for product in answer.products:
        for p in product:
            print(Chem.MolToSmiles(p))
    answer.visualizing()
        
    

if __name__ == "__main__":
    main() 