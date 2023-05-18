import tkinter as tk
from tkinter import ttk
from tab_interfaces.general_interface import Interface, get_folder_contents
# from general_interface import Interface, get_folder_contents

class ANNOVAR(Interface):
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
        input_filtered_variant_label = ttk.Label(input_panel, text="Input unfiltered variant calls:")
        input_filtered_variant_label.grid(row=0, column=0, sticky=tk.W)
        input_filtered_variant_entry = ttk.Entry(input_panel)
        input_filtered_variant_entry.grid(row=0, column=1, padx=5, pady=5, sticky=tk.W+tk.E)

        ref_genome_label = ttk.Label(input_panel, text="Input aligned tumor BAM/SAM:")
        ref_genome_label.grid(row=1, column=0, sticky=tk.W)
        ref_genome_entry = ttk.Entry(input_panel)
        ref_genome_entry.grid(row=1, column=1, padx=5, pady=5, sticky=tk.W+tk.E)

        interval_list_label = ttk.Label(input_panel, text="Interval list:")
        interval_list_label.grid(row=2, column=0, sticky=tk.W)
        interval_list_entry = ttk.Entry(input_panel)
        interval_list_entry.grid(row=2, column=1, padx=5, pady=5, sticky=tk.W+tk.E)
    
        funcotator_data_sources_label = ttk.Label(input_panel, text="Popular variant sites:")
        funcotator_data_sources_label.grid(row=3, column=0, sticky=tk.W)
        funcotator_data_sources_entry = ttk.Entry(input_panel)
        funcotator_data_sources_entry.grid(row=3, column=1, padx=5, pady=5, sticky=tk.W+tk.E)

        funcotator_filt_ann_output_label = ttk.Label(input_panel, text="F1R2:")
        funcotator_filt_ann_output_label.grid(row=4, column=0, sticky=tk.W)
        funcotator_filt_ann_output_entry = ttk.Entry(input_panel)
        funcotator_filt_ann_output_entry.grid(row=4, column=1, padx=5, pady=5, sticky=tk.W+tk.E)

        out_format_label = ttk.Label(input_panel, text="Output file name:")
        out_format_label.grid(row=5, column=0, sticky=tk.W)
        out_format_entry = ttk.Entry(input_panel)
        out_format_entry.grid(row=5, column=1, padx=5, pady=5, sticky=tk.W+tk.E)

        # Create a button to process the entered data
        process_button = ttk.Button(input_panel, text="Process", 
                            command=lambda: self.process_data(input_filtered_variant_entry, ref_genome_entry, funcotator_data_sources_entry, funcotator_filt_ann_output_entry, out_format_entry, result_text))
        process_button.grid(row=6, column=0, columnspan=2, pady=10)

        # Create a label for displaying the result
        result_text = tk.Text(input_panel, width=200, height=50, relief="sunken")
        result_text.grid(row=7, column=0, columnspan=2, sticky=tk.W+tk.E+tk.N+tk.S)

        # Create a frame for the export button
        export_frame = ttk.Frame(input_panel)
        export_frame.grid(row=8, column=0, columnspan=2, pady=10, sticky=tk.E)

        # Create a button to export the result
        export_button = ttk.Button(export_frame, text="Export", command=lambda: self.export_data(result_text))
        export_button.pack(padx=5, pady=5)

        # Return the result_text widget
        return result_text
 
    def process_data(self, input_filtered_variant_entry, ref_genome_entry, funcotator_data_sources_entry, funcotator_filt_ann_output_entry, out_format_entry, result_text):
        input_filtered_variant = input_filtered_variant_entry.get()
        ref_genome = ref_genome_entry.get()
        funcotator_data_sources = funcotator_data_sources_entry.get()
        funcotator_filt_ann_output = funcotator_filt_ann_output_entry.get()
        out_format = out_format_entry.get()
        
        result = f"gatk FuncotatorDataSourceDownloader --somatic --validate_integrity --extract-after-download\n"
        
        result += f"gatk Funcotator --variant {input_filtered_variant} --reference {ref_genome} --data-sources-path {funcotator_data_sources} --output {funcotator_filt_ann_output} --output-file-format {out_format}"

        # http://www.openbioinformatics.org/annovar/download/0wgxR2rIVP/annovar.latest.tar.gz       
        # {curr_dir}/
        # gatk Funcotator --variant {input_filtered_variant}
        # --reference {ref_genome} \
        # --ref-version {ref_ver} \
        # --data-sources-path {funcotator_data_sources} \
        # --output {funcotator_filt_ann_output}
        # --output-file-format {out_format}

        # Clear previous text content
        result_text.delete('1.0', tk.END)

        # Insert the new text content
        result_text.insert(tk.END, result)

        # Disable editing of the text widget
        result_text.configure(state=tk.DISABLED)

    def create_tab(self, notebook, tab_name, folder_path):
        self.notebook = notebook
        return super().create_tab(tab_name, folder_path)