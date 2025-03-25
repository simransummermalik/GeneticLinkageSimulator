import tkinter as tk
from tkinter import filedialog, messagebox
from chromosome import Chromosome
from recombination import simulate_crossover
from chi_square import chi_square_linkage
from visualize import plot_linkage_map
from ml_model import predict_snp_recombination

class GeneticSimulatorGUI:
    def __init__(self, root):
        self.root = root
        self.root.title("Genetic Linkage Simulator")
        self.root.geometry("600x500")

        tk.Label(root, text="Genetic Linkage Simulator", font=("Arial", 14, "bold")).pack(pady=10)

        self.snp_file_button = tk.Button(root, text="Upload SNP Data", command=self.upload_snp_file)
        self.snp_file_button.pack(pady=5)

        tk.Label(root, text="Enter SNP Position for Prediction:").pack()
        self.snp_entry = tk.Entry(root, width=40)
        self.snp_entry.pack()
        self.snp_predict_button = tk.Button(root, text="Predict SNP Recombination", command=self.predict_snp_recombination)
        self.snp_predict_button.pack(pady=5)

    def upload_snp_file(self):
        """Allows the user to upload a VCF file."""
        file_path = filedialog.askopenfilename(filetypes=[("VCF Files", "*.vcf")])
        if file_path:
            messagebox.showinfo("SNP Data", f"SNP data loaded from {file_path}")

    def predict_snp_recombination(self):
        """Predicts recombination likelihood based on user input SNP position."""
        try:
            position = int(self.snp_entry.get())
            result = predict_snp_recombination(position)
            messagebox.showinfo("Recombination Prediction", f"Predicted Recombination Rate: {result[0]}")
        except:
            messagebox.showerror("Error", "Invalid SNP position. Enter a number.")

if __name__ == "__main__":
    root = tk.Tk()
    app = GeneticSimulatorGUI(root)
    root.mainloop()
