import tkinter as tk
from tkinter import ttk
from tab_interfaces.general_interface import Interface, get_folder_contents
# from general_interface import Interface, get_folder_contents

class BWA(Interface):
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
        num_threads_label = ttk.Label(input_panel, text="Number of Threads:")
        num_threads_label.grid(row=0, column=0, sticky=tk.W)
        num_threads_entry = ttk.Entry(input_panel)
        num_threads_entry.grid(row=0, column=1, padx=5, pady=5, sticky=tk.W+tk.E)
        
        sra_id_label = ttk.Label(input_panel, text="SRA-ID:")
        sra_id_label.grid(row=1, column=0, sticky=tk.W)
        sra_id_entry = ttk.Entry(input_panel)
        sra_id_entry.grid(row=1, column=1, padx=5, pady=5, sticky=tk.W+tk.E)

        platform_label = ttk.Label(input_panel, text="Platform:")
        platform_label.grid(row=2, column=0, sticky=tk.W)
        platform_entry = ttk.Entry(input_panel)
        platform_entry.grid(row=2, column=1, padx=5, pady=5, sticky=tk.W+tk.E)

        ref_genome_label = ttk.Label(input_panel, text="Reference genome:")
        ref_genome_label.grid(row=3, column=0, sticky=tk.W)
        ref_genome_entry = ttk.Entry(input_panel)
        ref_genome_entry.grid(row=3, column=1, padx=5, pady=5, sticky=tk.W+tk.E)
        
        fwd_label = ttk.Label(input_panel, text="Forward read:")
        fwd_label.grid(row=4, column=0, sticky=tk.W)
        fwd_entry = ttk.Entry(input_panel)
        fwd_entry.grid(row=4, column=1, padx=5, pady=5, sticky=tk.W+tk.E)

        bwd_label = ttk.Label(input_panel, text="Backward reads:")
        bwd_label.grid(row=5, column=0, sticky=tk.W)
        bwd_entry = ttk.Entry(input_panel)
        bwd_entry.grid(row=5, column=1, padx=5, pady=5, sticky=tk.W+tk.E)

        output_label = ttk.Label(input_panel, text="Output Directory:")
        output_label.grid(row=6, column=0, sticky=tk.W)
        output_entry = ttk.Entry(input_panel)
        output_entry.grid(row=6, column=1, padx=5, pady=5, sticky=tk.W+tk.E)

        # Create a button to process the entered data
        process_button = ttk.Button(input_panel, text="Process", 
                            command=lambda: self.process_data(num_threads_entry, sra_id_entry, platform_entry, ref_genome_entry, fwd_entry, bwd_entry, output_entry, result_text))
        process_button.grid(row=7, column=0, columnspan=2, pady=10)

        # Create a label for displaying the result
        result_text = tk.Text(input_panel, width=200, height=50, relief="sunken")
        result_text.grid(row=8, column=0, columnspan=2, sticky=tk.W+tk.E+tk.N+tk.S)

        # Create a frame for the export button
        export_frame = ttk.Frame(input_panel)
        export_frame.grid(row=9, column=0, columnspan=2, pady=10, sticky=tk.E)

        # Create a button to export the result
        export_button = ttk.Button(export_frame, text="Export", command=lambda: self.export_data(result_text))
        export_button.pack(padx=5, pady=5)

        # Return the result_text widget
        return result_text
 
    
    def process_data(self, num_threads_entry, sra_id_entry, platform_entry, ref_genome_entry, fwd_entry, bwd_entry, output_entry, result_text):
        num_threads = num_threads_entry.get()
        sra_id =  sra_id_entry.get()
        platform = platform_entry.get()
        ref_genome = ref_genome_entry.get()
        fwd = fwd_entry.get()
        bwd = bwd_entry.get()
        output = output_entry.get()
        
        rg_info = '\"@RG\\tID:'+ sra_id + '\\tPL:' + platform + '\\tSM:'+ sra_id +'\"'
        
        result = f"bwa mem -t {num_threads} {rg_info} {ref_genome} {fwd} {bwd} {output}"

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
    interface = BWA(root)
    interface.create_tab("fastQC", "/home/hhbach/gui/data")
    interface.create_tab("BWA-MEM", "/path/to/folder2")
    root.mainloop()

if __name__ == "__main__":
    main()