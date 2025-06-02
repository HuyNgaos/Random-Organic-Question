import PyInstaller.__main__

PyInstaller.__main__.run([
    'main.py',              # Replace with your actual script
    '--onefile',            # Bundle into one executable
    '--noconsole'          # Hides the console (use this for GUI apps),
    # rememer to add , if you want to have icon
    # '--icon=Icon.ico'       # Optional: adds a custom icon
])