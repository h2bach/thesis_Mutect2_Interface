import os
import tkinter as tk
from tkinter import ttk, filedialog


def get_folder_contents(folder_path):
    # Retrieve the contents of the folder
    if os.path.exists(folder_path):
        return os.listdir(folder_path)
    else:
        return []


class Interface:
    def __init__(self, parent, folder_path="/home/hhbach/gui/"):
        self.parent = parent
        # Create a notebook (tabbed interface)
        self.notebook = ttk.Notebook(parent)
        self.notebook.pack(side=tk.RIGHT,fill=tk.BOTH, expand=True)
        self.folder_path = folder_path

    def create_input_panel(self, tab):
        # Create a frame for the input panel
        input_panel = ttk.Frame(tab)
        input_panel.pack(side=tk.RIGHT, fill=tk.BOTH, padx=10, pady=10)

        # Configure column weights
        input_panel.grid_columnconfigure(1, weight=1)

        # Create input fields
        reference_genome_label = ttk.Label(input_panel, text="Reference Genome:")
        reference_genome_label.grid(row=0, column=0, sticky=tk.W)
        reference_genome_entry = ttk.Entry(input_panel)
        reference_genome_entry.grid(row=0, column=1, padx=5, pady=5, sticky=tk.W+tk.E)

        reads_folder_label = ttk.Label(input_panel, text="Reads Folder:")
        reads_folder_label.grid(row=1, column=0, sticky=tk.W)
        reads_folder_entry = ttk.Entry(input_panel)
        reads_folder_entry.grid(row=1, column=1, padx=5, pady=5, sticky=tk.W+tk.E)

        output_directory_label = ttk.Label(input_panel, text="Output Directory:")
        output_directory_label.grid(row=2, column=0, sticky=tk.W)
        output_directory_entry = ttk.Entry(input_panel)
        output_directory_entry.grid(row=2, column=1, padx=5, pady=5, sticky=tk.W+tk.E)

        # Create a button to process the entered data
        process_button = ttk.Button(input_panel, text="Process", command=lambda: self.process_data(reference_genome_entry, reads_folder_entry, output_directory_entry, result_text))
        process_button.grid(row=3, column=0, columnspan=2, pady=10)

        # Create a label for displaying the result
        result_text = tk.Text(input_panel, width=200, height=50, relief="sunken")
        result_text.grid(row=4, column=0, columnspan=2, sticky=tk.W+tk.E+tk.N+tk.S)

        # Create a frame for the export button
        export_frame = ttk.Frame(input_panel)
        export_frame.grid(row=5, column=0, columnspan=2, pady=10, sticky=tk.E)

        # Create a button to export the result
        export_button = ttk.Button(export_frame, text="Export", command=lambda: self.export_data(result_text))
        export_button.pack(padx=5, pady=5)

        # Return the result_text widget
        return result_text



    def export_data(self, result_text):
        # Get the content of the result text widget
        content = result_text.get("1.0", tk.END)

        # Open a file dialog to select the export file
        file_path = filedialog.asksaveasfilename(defaultextension=".sh", filetypes=[("Bash Files", "*.sh")])

        # Write the content to the selected file
        with open(file_path, "w") as file:
            file.write(content)
    
    def process_data(self, reference_genome_entry, reads_folder_entry, output_directory_entry, result_text):
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

    def create_tab(self, tab_name, folder_path):
        # Create a tab
        tab = ttk.Frame(self.parent)

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
        browse_button = ttk.Button(tab, text="Browse", command=lambda: self.handle_browse_result(folder_path_entry, folder_path_label, tree))
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
        tree.bind("<<TreeviewSelect>>", lambda event: self.on_select(event, tree, path_label, folder_path))

        # Pack the Treeview widget
        tree.pack(expand=True, fill=tk.BOTH)
        self.create_input_panel(tab)

        # Add the tab to the notebook
        self.notebook.add(tab, text=tab_name)

    def handle_browse_result(self, folder_path_entry, folder_path_label, tree):
        # Open the file dialog to select a folder
        folder_selected = filedialog.askdirectory()

        # Update the folder path
        folder_path = folder_selected

        # Update the folder path entry field
        folder_path_entry.delete(0, tk.END)
        folder_path_entry.insert(tk.END, folder_path)

        # Update the selected directory label
        folder_path_label.configure(text="Selected Directory: " + folder_path)

        # Update the folder contents in the table
        folder_contents = get_folder_contents(folder_path)
        tree.delete(*tree.get_children())
        for index, item in enumerate(folder_contents):
            item_path = os.path.join(folder_path, item)
            item_type = "File" if os.path.isfile(item_path) else "Folder"
            tree.insert("", index, text="", values=(item, item_type))
    
    def on_select(self, event, tree, path_label, folder_path):
    # Retrieve the selected item(s)
        selection = event.widget.selection()
        if selection:
            item_text = tree.item(selection[0])["values"][0]
            item_path = os.path.join(folder_path, item_text)
            path_label.configure(text="Selected Path: " + item_path)

def main():
    root = tk.Tk()
    root.title("Pipeline")
    interface = Interface(root)
    interface.create_tab("Tab 1", "/path/to/folder1")
    interface.create_tab("Tab 2", "/path/to/folder2")
    root.mainloop()


if __name__ == "__main__":
    main()

