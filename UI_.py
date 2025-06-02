import tkinter as tk
from tkinter import ttk
import Q_generator as Qg
from rdkit import Chem
import random
from PIL import ImageTk

available_reactions = ["Esterfication", "Ester_hydrolysis", "AcylHalideEsterification"]



class FlashCardApp():
    def __init__(self, master):
        self.master = master
        self.master.title("Flash Card App")

        self.intro_frame = ttk.Frame(master)
        self.intro_frame.grid(row=0, column=0)
        self.intro_label = ttk.Label(self.intro_frame, text="Welcome to the Flash Card App!\nThis app is still in development.\nPlease check back later for updates.")
        self.intro_label.grid(row=0, column=0, padx=10, pady=10)
        
        ttk.Button(self.intro_frame, text="Start", command=self.start_init).grid(row=1, column=0, sticky=tk.W)
        ttk.Button(self.intro_frame, text="Exit", command=self.exit_app).grid(row=1, column=1, sticky=tk.W)
        self.reaction = None

    def reaction_get(self, reaction_name):
            match reaction_name:
                case "Random":
                    self.reaction = random.choice(available_reactions)
                case reaction_name if reaction_name in available_reactions:
                    self.reaction = reaction_name
                case _:
                    print(f"Reaction '{reaction_name}' is not available.")

    def start(self):
        # self.intro_frame.destroy()
        self.start_frame = ttk.Frame(self.master)
        self.start_frame.grid(row=0, column=0)
        start_reaction_label = ttk.Label(self.start_frame, text="Select a reaction to review:")
        start_reaction_label.grid(row=0, column=0, padx=10, pady=10)
        # Create a combobox for selecting reactions
        reactions = ttk.Combobox(self.start_frame, values=available_reactions + ["Random"])
        # reactions.set("Random")
        reactions.bind("<<ComboboxSelected>>", lambda event: self.reaction_get(reactions.get()))
        reactions.grid(row=0, column=1, padx=10, pady=10)
        # test_button = ttk.Button(self.start_frame, text="Test", command=lambda: print(self.reaction))
        # test_button.grid(row=0, column=2, padx=10, pady=10)
        # Create a button to start the flash card display
        next_button = ttk.Button(self.start_frame, text="Next", command=self.flash_card_display)
        next_button.grid(row=1, column=0, columnspan=3, padx=10, pady=10)

    def start_init(self):
        self.intro_frame.destroy()
        self.start()

    def flash_card_data(self):
        self.question, self.answer = eval(f"Qg.{self.reaction}().generate_question()")

    def flash_card_display(self):
        self.start_frame.destroy()
        # Create the flash card display
        self.flash_card_frame_title = ttk.Frame(self.master)
        self.flash_card_frame_title.grid(row=0, column=0)
        title_label = ttk.Label(self.flash_card_frame_title, text=f"Flash Card: {self.reaction}")
        title_label.grid(row=0, column=1, padx=0, pady=10)
        self.flash_card_data()
        # Exit button
        exit_button = ttk.Button(self.flash_card_frame_title, text="Exit", command=self.exit_app)
        exit_button.grid(row=0, column = 2, sticky= "ne", pady=10)
        # Create the question frame
        self.flash_card_frame_q = ttk.Frame(self.master)
        question_label = ttk.Label(self.flash_card_frame_q, text=self.question)
        question_label.grid(row=0, column=0, padx=10, pady=10)
        self.flash_card_frame_q.grid(row=1, column=0)
        # Create the answer frame
        self.flash_card_frame_a = ttk.Frame(self.master)

        answer_string = "Answer:"
        for product in self.answer.products:
            for p in product:
                answer_string += ", " + str(Chem.MolToSmiles(p))

        answer_label = ttk.Label(self.flash_card_frame_a, text=answer_string)
        answer_label.grid(row=0, column=0, padx=10, pady=10)

        reaction_image = ImageTk.PhotoImage(self.answer.img_get())
        img_label = ttk.Label(self.flash_card_frame_a, image=reaction_image)
        img_label.image = reaction_image
        img_label.grid(row=1, column=0, padx=10, pady=10)
        #Visualize the reaction in photo
        visualize_button = ttk.Button(self.flash_card_frame_a, text="Open in default photo viewer", command=self.answer.visualizing)
        visualize_button.grid(row=2, column=0, padx=10, pady=10)
        # Create the flip button
        self.flip = False
        self.flip_button = ttk.Button(self.flash_card_frame_title, text="Flip", command=self.flip_card)
        self.flip_button.grid(row=0, column=0, padx=10, pady=10, sticky= "nw")
        # Create the reset button
        self.reset_button = ttk.Button(self.flash_card_frame_a, text="Reset", command=self.start_reset)
        self.reset_button.grid(row=3, column=0, padx=10, pady=10)

    def flip_card(self):
        self.flip = not self.flip
        if self.flip:
            self.flash_card_frame_q.grid_remove()
            self.flash_card_frame_a.grid(row=1, column=0)
        else:
            self.flash_card_frame_a.grid_remove()
            self.flash_card_frame_q.grid(row=1, column=0)

    def start_reset(self):
        """Reset the app to the initial state"""
        self.flash_card_frame_title.destroy()
        self.flash_card_frame_q.destroy()
        self.flash_card_frame_a.destroy()
        self.start()

        

    def exit_app(self):
        self.master.quit()

