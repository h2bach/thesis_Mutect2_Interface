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
        unfiltered_vcf_label = ttk.Label(input_panel, text="Input unfiltered variant calls:")
        unfiltered_vcf_label.grid(row=0, column=0, sticky=tk.W)
        unfiltered_vcf_entry = ttk.Entry(input_panel)
        unfiltered_vcf_entry.grid(row=0, column=1, padx=5, pady=5, sticky=tk.W+tk.E)

        tumor_bam_label = ttk.Label(input_panel, text="Input aligned tumor BAM/SAM:")
        tumor_bam_label.grid(row=1, column=0, sticky=tk.W)
        tumor_bam_entry = ttk.Entry(input_panel)
        tumor_bam_entry.grid(row=1, column=1, padx=5, pady=5, sticky=tk.W+tk.E)

        interval_list_label = ttk.Label(input_panel, text="Interval list:")
        interval_list_label.grid(row=2, column=0, sticky=tk.W)
        interval_list_entry = ttk.Entry(input_panel)
        interval_list_entry.grid(row=2, column=1, padx=5, pady=5, sticky=tk.W+tk.E)
    
        pop_variant_sites_label = ttk.Label(input_panel, text="Popular variant sites:")
        pop_variant_sites_label.grid(row=3, column=0, sticky=tk.W)
        pop_variant_sites_entry = ttk.Entry(input_panel)
        pop_variant_sites_entry.grid(row=3, column=1, padx=5, pady=5, sticky=tk.W+tk.E)

        f1r2_tar_gz_label = ttk.Label(input_panel, text="F1R2:")
        f1r2_tar_gz_label.grid(row=4, column=0, sticky=tk.W)
        f1r2_tar_gz_entry = ttk.Entry(input_panel)
        f1r2_tar_gz_entry.grid(row=4, column=1, padx=5, pady=5, sticky=tk.W+tk.E)

        output_vcf_label = ttk.Label(input_panel, text="Output file name:")
        output_vcf_label.grid(row=5, column=0, sticky=tk.W)
        output_vcf_entry = ttk.Entry(input_panel)
        output_vcf_entry.grid(row=5, column=1, padx=5, pady=5, sticky=tk.W+tk.E)

        # Create a button to process the entered data
        process_button = ttk.Button(input_panel, text="Process", 
                            command=lambda: self.process_data(unfiltered_vcf_entry, tumor_bam_entry, interval_list_entry, pop_variant_sites_entry, f1r2_tar_gz_entry, output_vcf_entry,result_text))
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
 
    
    def process_data(self, unfiltered_vcf_entry, tumor_bam_entry, interval_list_entry, pop_variant_sites_entry, f1r2_tar_gz_entry, output_vcf_entry,result_text):
        unfiltered_vcf = unfiltered_vcf_entry.get()
        tumor_bam = tumor_bam_entry.get()
        interval_list = interval_list_entry.get()
        pop_variant_sites = pop_variant_sites_entry.get()
        f1r2_tar_gz =f1r2_tar_gz_entry.get()
        output_vcf = output_vcf_entry.get()
        result = f"gatk LearnOrientationModel -I {f1r2_tar_gz} -O read-orientation-model.tar.gz\n"
        result += f"gatk GetPileupSummaries -I {tumor_bam} -V {pop_variant_sites} -L {interval_list} -O pileupsummaries.table\n"
        result += f"gatk CalculateContamination -I pileupsummaries.table -tumor-segmentation segments.table -O contamination.table\n"
        result += f"gatk FilterMutectCalls -V {unfiltered_vcf} --tumor-segmentation segments.table --contamination-table contamination.table --ob-priors read-orientation-model.tar.gz -O {output_vcf}"

        # gatk FilterMutectCalls -V unfiltered.vcf \
        # --tumor-segmentation segments.table \
        # --contamination-table contamination.table] \
        # --ob-priors read-orientation-model.tar.gz \
        # -O filtered.vcf

        # Clear previous text content
        result_text.delete('1.0', tk.END)

        # Insert the new text content
        result_text.insert(tk.END, result)

        # Disable editing of the text widget
        result_text.configure(state=tk.DISABLED)

    def create_tab(self, notebook, tab_name, folder_path):
        self.notebook = notebook
        return super().create_tab(tab_name, folder_path)