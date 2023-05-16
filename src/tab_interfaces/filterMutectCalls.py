import tkinter as tk
from tkinter import ttk
from tab_interfaces.general_interface import Interface, get_folder_contents
# from general_interface import Interface, get_folder_contents

class FilterMutectCalls(Interface):
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
        input_tumor_label = ttk.Label(input_panel, text="Input aligned tumor BAM/SAM:")
        input_tumor_label.grid(row=0, column=0, sticky=tk.W)
        input_tumor_entry = ttk.Entry(input_panel)
        input_tumor_entry.grid(row=0, column=1, padx=5, pady=5, sticky=tk.W+tk.E)

        input_normal_label = ttk.Label(input_panel, text="Input aligned normal BAM/SAM:")
        input_normal_label.grid(row=1, column=0, sticky=tk.W)
        input_normal_entry = ttk.Entry(input_panel)
        input_normal_entry.grid(row=1, column=1, padx=5, pady=5, sticky=tk.W+tk.E)

        ref_genome_label = ttk.Label(input_panel, text="Reference genome:")
        ref_genome_label.grid(row=2, column=0, sticky=tk.W)
        ref_genome_entry = ttk.Entry(input_panel)
        ref_genome_entry.grid(row=2, column=1, padx=5, pady=5, sticky=tk.W+tk.E)

        interval_list_label = ttk.Label(input_panel, text="Interval list:")
        interval_list_label.grid(row=3, column=0, sticky=tk.W)
        interval_list_entry = ttk.Entry(input_panel)
        interval_list_entry.grid(row=3, column=1, padx=5, pady=5, sticky=tk.W+tk.E)
    
        germline_resource_label = ttk.Label(input_panel, text="Germline resources:")
        germline_resource_label.grid(row=4, column=0, sticky=tk.W)
        germline_resource_entry = ttk.Entry(input_panel)
        germline_resource_entry.grid(row=4, column=1, padx=5, pady=5, sticky=tk.W+tk.E)

        pon_label = ttk.Label(input_panel, text="Panel of Normals:")
        pon_label.grid(row=5, column=0, sticky=tk.W)
        pon_entry = ttk.Entry(input_panel)
        pon_entry.grid(row=5, column=1, padx=5, pady=5, sticky=tk.W+tk.E)

        f1r2_tar_gz_label = ttk.Label(input_panel, text="F1R2:")
        f1r2_tar_gz_label.grid(row=6, column=0, sticky=tk.W)
        f1r2_tar_gz_entry = ttk.Entry(input_panel)
        f1r2_tar_gz_entry.grid(row=6, column=1, padx=5, pady=5, sticky=tk.W+tk.E)

        output_label = ttk.Label(input_panel, text="Output Directory:")
        output_label.grid(row=7, column=0, sticky=tk.W)
        output_entry = ttk.Entry(input_panel)
        output_entry.grid(row=7, column=1, padx=5, pady=5, sticky=tk.W+tk.E)

        output_bam_label = ttk.Label(input_panel, text="BAM-output Directory:")
        output_bam_label.grid(row=8, column=0, sticky=tk.W)
        output_bam_entry = ttk.Entry(input_panel)
        output_bam_entry.grid(row=8, column=1, padx=5, pady=5, sticky=tk.W+tk.E)

        # Create a button to process the entered data
        process_button = ttk.Button(input_panel, text="Process", 
                            command=lambda: self.process_data(input_tumor_entry, input_normal_entry, interval_list_entry, ref_genome_entry, germline_resource_entry, pon_entry, f1r2_tar_gz_entry, output_entry, output_bam_entry, result_text))
        process_button.grid(row=9, column=0, columnspan=2, pady=10)

        # Create a label for displaying the result
        result_text = tk.Text(input_panel, width=200, height=50, relief="sunken")
        result_text.grid(row=10, column=0, columnspan=2, sticky=tk.W+tk.E+tk.N+tk.S)

        # Create a frame for the export button
        export_frame = ttk.Frame(input_panel)
        export_frame.grid(row=11, column=0, columnspan=2, pady=10, sticky=tk.E)

        # Create a button to export the result
        export_button = ttk.Button(export_frame, text="Export", command=lambda: self.export_data(result_text))
        export_button.pack(padx=5, pady=5)

        # Return the result_text widget
        return result_text
 
    
    def process_data(self, input_tumor_entry, input_normal_entry, interval_list_entry, ref_genome_entry, germline_resource_entry, pon_entry, f1r2_tar_gz_entry, output_entry, output_bam_entry,result_text):
        input_tumor = input_tumor_entry.get()
        input_normal = input_normal_entry.get()
        interval_list = interval_list_entry.get()
        ref_genome = ref_genome_entry.get()
        germline_resource = germline_resource_entry.get()
        pon = pon_entry.get()
        f1r2_tar_gz =f1r2_tar_gz_entry.get()
        output = output_entry.get()
        output_bam = output_bam_entry.get()
        result = f"gatk Mutect2 -R {ref_genome} -L {interval_list} -I {input_tumor} -I {input_normal} --germline-resource {germline_resource} -pon {pon} --f1r2-tar-gz {f1r2_tar_gz} -O {output} -bamout {output_bam}"

        # Clear previous text content
        result_text.delete('1.0', tk.END)

        # Insert the new text content
        result_text.insert(tk.END, result)

        # Disable editing of the text widget
        result_text.configure(state=tk.DISABLED)

    def create_tab(self, notebook, tab_name, folder_path):
        self.notebook = notebook
        return super().create_tab(tab_name, folder_path)