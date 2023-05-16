import os
import tkinter as tk
from tkinter import ttk, filedialog

def get_folder_contents(folder_path):
    # Retrieve the contents of the folder
    if os.path.exists(folder_path):
        return os.listdir(folder_path)
    else:
        return []

def handle_browse_result(folder_path_entry, tree, folder_path_label):

    # Open the file dialog to select a folder
    folder_selected = filedialog.askdirectory()

    # Update the folder path
    folder_path = folder_selected

    # Update the folder path entry field
    folder_path_entry.delete(0, tk.END)
    folder_path_entry.insert(tk.END, folder_path)

    # Update the folder contents in the table
    folder_contents = get_folder_contents(folder_path)
    tree.delete(*tree.get_children())
    for index, item in enumerate(folder_contents):
        item_path = os.path.join(folder_path, item)
        item_type = "File" if os.path.isfile(item_path) else "Folder"
        tree.insert("", index, text="", values=(item, item_type))

    # Update the selected folder path label
    folder_path_label.config(text=f"Selected Directory: {folder_path}")

def on_select(event, tree, path_label, folder_path):
    # Retrieve the selected item(s)
    selection = event.widget.selection()
    if selection:
        item_text = tree.item(selection[0])["values"][0]
        item_path = os.path.join(folder_path, item_text)
        path_label.configure(text="Selected Path: " + item_path)

def process_data(reference_genome_entry, reads_folder_entry, output_directory_entry, result_text):
    # Retrieve the entered values from the input fields
    reference_genome = reference_genome_entry.get()
    reads_folder = reads_folder_entry.get()
    output_directory = output_directory_entry.get()

    # Format the result
    result = f"Reference Genome: {reference_genome}\n"
    result += f"Reads Folder: {reads_folder}\n"
    result += f"Output Directory: {output_directory}"

    # Clear previous text content
    result_text.delete('1.0', tk.END)

    # Insert the new text content
    result_text.insert(tk.END, result)

    # Disable editing of the text widget
    result_text.configure(state=tk.DISABLED)

def create_input_panel(parent):
    # Create input fields
    reference_genome_label = ttk.Label(parent, text="Reference Genome:")
    reference_genome_label.pack(anchor=tk.W)
    reference_genome_entry = ttk.Entry(parent)
    reference_genome_entry.pack(fill=tk.X, padx=5, pady=5)

    reads_folder_label = ttk.Label(parent, text="Reads Folder:")
    reads_folder_label.pack(anchor=tk.W)
    reads_folder_entry = ttk.Entry(parent)
    reads_folder_entry.pack(fill=tk.X, padx=5, pady=5)

    output_directory_label = ttk.Label(parent, text="Output Directory:")
    output_directory_label.pack(anchor=tk.W)
    output_directory_entry = ttk.Entry(parent)
    output_directory_entry.pack(fill=tk.X, padx=5, pady=5)

    # Create a button to process the entered data
    process_button = ttk.Button(parent, text="Process",
                                command=lambda: process_data(reference_genome_entry,
                                                             reads_folder_entry, output_directory_entry,
                                                             result_text))
    process_button.pack(pady=10)

    # Create a label for displaying the result
    result_text = tk.Text(parent, width=200, height=50, relief="sunken")
    result_text.pack()

    # Return the result_text widget
    return result_text

def create_tab(window, folder_path):
    # Create a frame for the tab
    tab = ttk.Frame(window)

    # Create an entry field to display the selected folder path
    folder_path_entry = ttk.Entry(tab, width=50)
    folder_path_entry.pack(anchor=tk.NW, padx=10, pady=10)

    # Create a side panel for the table
    side_panel = ttk.Frame(tab)
    side_panel.pack(side=tk.LEFT, fill=tk.Y)

    # Create a Treeview widget for the table
    tree = ttk.Treeview(side_panel, selectmode="extended")

    # Define the columns
    tree["columns"] = ("Name", "Type")

    # Format the columns
    tree.column("#0", width=0, stretch=tk.NO)
    tree.column("Name", width=200, anchor=tk.W)
    tree.column("Type", width=100, anchor=tk.W)

    # Create a button to browse for folder
    browse_button = ttk.Button(tab, text="Browse", command=lambda: handle_browse_result(folder_path_entry, tree, folder_path_label))
    browse_button.pack(anchor=tk.NW, padx=10, pady=10)

    # Create a label to display the selected directory
    folder_path_label = ttk.Label(tab, text="Selected Directory: ")
    folder_path_label.pack(anchor=tk.SE)

    # Populate the table with folder contents
    contents = get_folder_contents(folder_path)
    for index, item in enumerate(contents):
        item_path = os.path.join(folder_path, item)
        item_type = "File" if os.path.isfile(item_path) else "Folder"
        tree.insert("", index, text="", values=(item, item_type))

    # Create a frame for the path display box
    path_frame = ttk.Frame(side_panel)
    path_frame.pack(pady=10)

    # Create a label for displaying the path
    path_label = ttk.Label(path_frame, text="Selected Path: ")
    path_label.pack(anchor=tk.E)

    # Bind the selection event
    tree.bind("<<TreeviewSelect>>", lambda event: on_select(event, tree, path_label, folder_path))

    # Pack the Treeview widget
    tree.pack(expand=True, fill=tk.BOTH)

    # Create the main panel for input fields
    main_panel = ttk.Frame(tab)
    main_panel.pack(side=tk.RIGHT, fill=tk.BOTH, expand=True, padx=10, pady=10)

    create_input_panel(main_panel)

    return tab


def main():
    window = tk.Tk()
    window.title("Pipeline")
    # Create a notebook (tabbed interface)
    notebook = ttk.Notebook(window)
    notebook.pack(fill=tk.BOTH, expand=True)

    # Create a tab
    folder_path = "path/to/your/folder"
    tab = create_tab(notebook, folder_path)
    tab2 = create_tab(notebook, folder_path)
    tab3 = create_tab(notebook, folder_path)

    # Add the tab to the notebook
    notebook.add(tab, text="fastQC")
    notebook.add(tab2, text="Tab 1")
    notebook.add(tab3, text="Tab 1")

    window.mainloop()

if __name__ == "__main__":
    main()

