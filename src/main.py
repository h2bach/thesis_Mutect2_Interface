import tkinter as tk
import tab_interfaces.general_interface as Interface
import tab_interfaces.fastqc as FastQC
import tab_interfaces.bwa as BWA

def main():
    root = tk.Tk()
    root.title("Pipeline")
    notebook = tk.ttk.Notebook(root)
    notebook.pack(fill=tk.BOTH, expand=True)

    fastqc = FastQC.fastQC(notebook,root)
    fastqc.create_tab(notebook, "fastQC", fastqc.folder_path)
    bwa = BWA.BWA(notebook,root)
    bwa.create_tab(notebook, "BWA-MEM",bwa.folder_path)
    root.mainloop()

if __name__ == "__main__":
    main()