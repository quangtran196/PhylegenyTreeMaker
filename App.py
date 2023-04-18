import tkinter as tk
from tkinter import filedialog
from io import StringIO
from Bio import SeqIO, AlignIO, Phylo
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
import matplotlib.pyplot as plt


class App(tk.Tk):
    def __init__(self):
        super().__init__()

        self.title("Sequence Alignment and Phylogenetic Tree Construction")
        self.geometry("500x500")

        self.text = tk.Text(self, width=50, height=10)
        self.text.pack(pady=10)

        self.upload_button = tk.Button(self, text="Upload File", command=self.upload_file)
        self.upload_button.pack()

        self.submit_button = tk.Button(self, text="Submit", command=self.submit)
        self.submit_button.pack(pady=10)

    def upload_file(self):
        self.file_path = filedialog.askopenfilename()
        file = open(self.file_path)
        content = file.read()
        self.text.insert(tk.END, content)
        file.close()

    def submit(self):
        content = self.text.get("1.0", tk.END)
        alignment = align_sequences(content)
        draw_tree(alignment)


def align_sequences(content):
    # Parse the input sequences
    records = list(SeqIO.parse(StringIO(content), "fasta"))

    # Find the length of the longest sequence
    max_length = max([len(record.seq) for record in records])

    # Pad the sequences to make them the same length
    for record in records:
        seq = record.seq
        if len(seq) < max_length:
            seq += "-" * (max_length - len(seq))
        record.seq = seq

    # Create a MultipleSeqAlignment object
    alignment = MultipleSeqAlignment(records)

    return alignment


def draw_tree(alignment):
    # Calculate pairwise distances
    calculator = DistanceCalculator('identity')
    dm = calculator.get_distance(alignment)

    # Construct a neighbor-joining tree
    constructor = DistanceTreeConstructor(calculator)
    tree = constructor.nj(dm)

    # Plot tree using matplotlib
    fig, ax = plt.subplots(figsize=(8, 8))
    Phylo.draw(tree, axes=ax, show_confidence=False)
    plt.show()


if __name__ == '__main__':
    app = App()
    app.mainloop()
