import random
import numpy as np
import csv
import tkinter as tk
from tkinter import filedialog, messagebox
import matplotlib.pyplot as plt
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
from scipy.stats import chi2_contingency
class Chromosome:
    def __init__(self, genes, positions):
        self.genes = genes
        self.positions = positions

    def display(self): 
        print("\nChromosome Gene Map:")
        for gene, pos in zip(self.genes, self.positions):
            print(f"Gene {gene}: {pos} cM")
def simulate_crossover(chromosome):
    crossovers = []
    for i in range(len(chromosome.genes) - 1):
        dist = chromosome.positions[i+1] - chromosome.positions[i]
        prob = 1 - np.exp(-2 * dist / 100)
        if random.random() < prob:
            crossovers.append(chromosome.genes[i])
    return crossovers
def chi_square_linkage(observed_counts):
    total = sum(observed_counts)
    expected = [total / len(observed_counts)] * len(observed_counts)
    chi2, p_value = chi2_contingency([observed_counts, expected])[:2]
    return p_value < 0.05
def plot_linkage_map(chromosome):
    positions = chromosome.positions
    genes = chromosome.genes
    plt.figure(figsize=(8, 2))
    plt.scatter(positions, [1] * len(positions), marker="|", color="red")
    for i, gene in enumerate(genes):
        plt.text(positions[i], 1.02, gene, ha="center")
    plt.title("Genetic Linkage Map")
    plt.xlabel("Genetic Distance (cM)")
    plt.yticks([])
    plt.show()
def read_snp_csv(file_path):
    snps = []
    try:
        with open(file_path, 'r') as file:
            reader = csv.reader(file)
            next(reader)  
            for row in reader:
                if len(row) < 4:
                    print(f"Skipping invalid row: {row}")
                    continue
                snps.append({
                    "chrom": row[0],
                    "position": int(row[1]),
                    "ref": row[2],
                    "alt": row[3]
                })
    except FileNotFoundError:
        print("Error: File not found.")
    return snps
def train_ml_with_snps(csv_file):
    snps = read_snp_csv(csv_file)
    if not snps:
        print("Error: No SNP data available.")
        return None
    data = np.array([[snp["position"]] for snp in snps])  
    labels = np.random.choice([0, 1], size=len(snps))
    X_train, X_test, y_train, y_test = train_test_split(data, labels, test_size=0.2)
    model = RandomForestClassifier()
    model.fit(X_train, y_train)
    return model

ml_model = train_ml_with_snps("example.csv")

def predict_snp_recombination(position):
    return ml_model.predict([[position]]) if ml_model else None
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
        file_path = filedialog.askopenfilename(filetypes=[("CSV Files", "*.csv")])
        if file_path:
            global ml_model
            ml_model = train_ml_with_snps(file_path)
            messagebox.showinfo("SNP Data", f"SNP data loaded from {file_path}")

    def predict_snp_recombination(self):
        try:
            position = int(self.snp_entry.get())
            result = predict_snp_recombination(position)
            if result is not None:
                messagebox.showinfo("Recombination Prediction", f"Predicted Recombination Rate: {result[0]}")
            else:
                messagebox.showerror("Error", "No trained model available.")
        except ValueError:
            messagebox.showerror("Error", "Invalid SNP position. Enter a number.")

if __name__ == "__main__":
    root = tk.Tk()
    app = GeneticSimulatorGUI(root)
    root.mainloop()