import UI_ as ui
import tkinter as tk

def main():
    introduction = """    Welcome to the Flash Card App!
This app is still in development.
Please check back later for updates."""
    print(introduction)
    
    root = tk.Tk() #start the tkinter application
    app = ui.FlashCardApp(root) #create the app
    root.mainloop()
    

if __name__ == "__main__":
    main() 