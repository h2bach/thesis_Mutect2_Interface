import tkinter as tk
from tkinter import ttk
from tab_interfaces.general_interface import Interface, get_folder_contents

class fastQC(Interface):
    def __init__(self, notebook, parent, folder_path="/home/hhbach/gui/"):
        self.notebook = notebook
        super().__init__(parent, folder_path)

    def create_input_panel(self, tab):
        # Create a frame for the input panel
        input_panel = ttk.Frame(tab)
        input_panel.pack(side=tk.RIGHT, fill=tk.BOTH, padx=10, pady=10)

        # Configure column weights
        input_panel.grid_columnconfigure(1, weight=1)

        # Create input fields
        input_directory_label = ttk.Label(input_panel, text="Input Directory:")
        input_directory_label.grid(row=0, column=0, sticky=tk.W)
        input_directory_entry = ttk.Entry(input_panel)
        input_directory_entry.grid(row=0, column=1, padx=5, pady=5, sticky=tk.W+tk.E)

        output_directory_label = ttk.Label(input_panel, text="Output Directory:")
        output_directory_label.grid(row=2, column=0, sticky=tk.W)
        output_directory_entry = ttk.Entry(input_panel)
        output_directory_entry.grid(row=2, column=1, padx=5, pady=5, sticky=tk.W+tk.E)

        # Create a button to process the entered data
        process_button = ttk.Button(input_panel, text="Process", command=lambda: self.process_data(input_directory_entry, output_directory_entry, result_text))
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
 
    
    def process_data(self, input_directory_entry, output_directory_entry, result_text):
        input_directory = input_directory_entry.get()
        output_directory = output_directory_entry.get()

        read_list = [read for read in get_folder_contents(input_directory) if '.fastq' in read]
        print(read_list)
        
        # Format the result
        result = f"fastqc"
        for read in read_list:
            result += f" {read}"
        
        result += f" -o {output_directory}"

        # Clear previous text content
        result_text.delete('1.0', tk.END)

        # Insert the new text content
        result_text.insert(tk.END, result)

        # Disable editing of the text widget
        result_text.configure(state=tk.DISABLED)
    
    def create_tab(self, notebook, tab_name, folder_path):
        self.notebook = notebook
        return super().create_tab(tab_name, folder_path)

def main():
    root = tk.Tk()
    root.title("Pipeline")
    interface = fastQC(root)
    interface.create_tab("fastQC", "/home/hhbach/gui/data")
    interface.create_tab("Tab 2", "/path/to/folder2")
    root.mainloop()

if __name__ == "__main__":
    main()