import tkinter as tk
from tkinter import ttk
import ui.interfaces.interface as i
def main():
    window = tk.Tk()
    window.title("Pipeline")
    # Create a notebook (tabbed interface)
    notebook = ttk.Notebook(window)
    notebook.pack(fill=tk.BOTH, expand=True)

    # Create a tab
    folder_path = "path/to/your/folder"
    tab = i.create_tab(notebook, folder_path)
    tab2 = i.create_tab(notebook, folder_path)
    tab3 = i.create_tab(notebook, folder_path)

    # Add the tab to the notebook
    notebook.add(tab, text="fastQC")
    notebook.add(tab2, text="Tab 1")
    notebook.add(tab3, text="Tab 1")

    window.mainloop()

if __name__ == "__main__":
    main()