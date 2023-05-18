import tkinter as tk
import tab_interfaces.general_interface as Interface
import tab_interfaces.fastqc as FastQC
import tab_interfaces.bwa as BWA
import tab_interfaces.markdup_sort as MarkDuplicatesSpark
import tab_interfaces.bqsr as BQSR
import tab_interfaces.mutect2 as Mutect2
import tab_interfaces.filterMutectCalls as FMC
import tab_interfaces.annotation_funcotator as Funcotator
import tab_interfaces.bash_run as BashRunner

def main():
    root = tk.Tk()
    root.title("Pipeline")
    notebook = tk.ttk.Notebook(root)
    notebook.pack(fill=tk.BOTH, expand=True)

    fastqc = FastQC.fastQC(notebook,root)
    fastqc.create_tab(notebook, "fastQC", fastqc.folder_path)
    bwa = BWA.BWA(notebook,root)
    bwa.create_tab(notebook, "BWA-MEM",bwa.folder_path)
    markDuplicatesSpark = MarkDuplicatesSpark.MarkDuplicatesSpark(notebook,root)
    markDuplicatesSpark.create_tab(notebook, "MarkDuplicatesSpark",bwa.folder_path)
    bqsr = BQSR.BQSR(notebook, root)
    bqsr.create_tab(notebook,"BQSR",bqsr.folder_path)
    mutect2 = Mutect2.Mutect2(notebook, root)
    mutect2.create_tab(notebook,"Mutect2",mutect2.folder_path)
    fmc = FMC.FilterMutectCalls(notebook, root)
    fmc.create_tab(notebook,"FilterMutectCalls", fmc.folder_path)
    func_ann = Funcotator.Funcotator(notebook, root)
    func_ann.create_tab(notebook,"Funcotator", fmc.folder_path)
    bash_run = BashRunner.BashRunner(notebook, root)
    bash_run.create_tab(notebook,"BashRunner", fmc.folder_path)
    root.mainloop()

if __name__ == "__main__":
    main()