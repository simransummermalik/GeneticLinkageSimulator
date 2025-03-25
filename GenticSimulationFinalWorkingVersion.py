import random
import numpy as np
import csv
import tkinter as tk
from tkinter import filedialog, messagebox
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
from scipy.stats import chi2_contingency
from scipy.stats import pearsonr

# Chromosome Class
class Chromosome:
    def __init__(self, genes, positions):
        self.genes = genes
        self.positions = positions

    def display(self): 
        print("\nChromosome Gene Map:")
        for gene, pos in zip(self.genes, self.positions):
            print(f"Gene {gene}: {pos} cM")

# Simulate Crossover
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

# Plot Linkage Map with Color-coded Chromosomes and Tooltips
def plot_linkage_map(snps):
    chromosomes = {}
    for snp in snps:
        chrom = snp['chrom']
        if chrom not in chromosomes:
            chromosomes[chrom] = []
        chromosomes[chrom].append(snp['position'])
    
    plt.figure(figsize=(10, 6))
    colors = plt.cm.get_cmap('tab10', len(chromosomes))
    for i, (chrom, positions) in enumerate(chromosomes.items()):
        plt.scatter(positions, [i] * len(positions), marker="|", color=colors(i), label=f"Chromosome {chrom}")
    
    plt.title("Genetic Linkage Map")
    plt.xlabel("Genetic Distance (cM)")
    plt.ylabel("Chromosomes")
    plt.legend()
    plt.show()

# Heatmap of Recombination Frequencies
def plot_recombination_heatmap(snps):
    recomb_matrix = np.zeros((len(snps), len(snps)))
    for i in range(len(snps)):
        for j in range(len(snps)):
            recomb_matrix[i, j] = abs(snps[i]['position'] - snps[j]['position'])
    
    plt.figure(figsize=(8, 6))
    sns.heatmap(recomb_matrix, cmap="coolwarm", annot=False)
    plt.title("Recombination Frequency Heatmap")
    plt.show()

# Read SNP Data from CSV
def read_snp_csv(file_path):
    snps = []
    with open(file_path, 'r') as file:
        reader = csv.reader(file)
        next(reader)
        for row in reader:
            if not row or len(row) < 4:
                continue
            snps.append({
                "chrom": row[0],
                "position": int(row[1]),
                "ref": row[2],
                "alt": row[3]
            })
    return snps

# Correlation Analysis
def correlation_analysis(snps):
    positions = [snp['position'] for snp in snps]
    labels = np.random.choice([0, 1], size=len(snps))
    correlation, _ = pearsonr(positions, labels)
    return correlation

# Train ML Model with SNP Data
def train_ml_with_snps(csv_file):
    snps = read_snp_csv(csv_file)
    data = np.array([[snp["position"]] for snp in snps])  
    labels = np.random.choice([0, 1], size=len(snps))
    X_train, X_test, y_train, y_test = train_test_split(data, labels, test_size=0.2)
    model = RandomForestClassifier()
    model.fit(X_train, y_train)
    return model

ml_model = train_ml_with_snps("example.csv")

def predict_snp_recombination(position):
    return ml_model.predict([[position]])
def haplotype_block_detection(snps):
    haplotype_blocks = {}
    for snp in snps:
        chrom = snp['chrom']
        if chrom not in haplotype_blocks:
            haplotype_blocks[chrom] = []
        haplotype_blocks[chrom].append(snp['position'])
    return haplotype_blocks
def perform_haplotype_block_detection(self):
        blocks = haplotype_block_detection(self.snps)
        messagebox.showinfo("Haplotype Blocks", f"Detected Haplotype Blocks: {blocks}")
# Heatmap of Recombination Frequencies
def plot_recombination_heatmap(snps):
    recomb_matrix = np.zeros((len(snps), len(snps)))
    for i in range(len(snps)):
        for j in range(len(snps)):
            recomb_matrix[i, j] = abs(snps[i]['position'] - snps[j]['position'])
    
    plt.figure(figsize=(8, 6))
    sns.heatmap(recomb_matrix, cmap="coolwarm", annot=False)
    plt.title("Recombination Frequency Heatmap")
    plt.show()
# GUI
class GeneticSimulatorGUI:
    def __init__(self, root):
        self.root = root
        self.root.title("Genetic Linkage Simulator")
        self.root.geometry("900x800")
        self.snp_file_path = ""
        self.snps = []
        
        tk.Label(root, text="Genetic Linkage Simulator", font=("Arial", 16, "bold")).pack(pady=10)
        self.info_label = tk.Label(root, text="This tool helps visualize genetic linkage and SNP recombination.", font=("Arial", 10))
        self.info_label.pack()
        
        self.snp_file_button = tk.Button(root, text="Upload SNP Data", command=self.upload_snp_file)
        self.snp_file_button.pack(pady=5)
        
        self.snp_listbox = tk.Listbox(root, height=10, width=80)
        self.snp_listbox.pack()
        
        tk.Label(root, text="Enter SNP Position for Prediction:").pack()
        self.snp_entry = tk.Entry(root, width=40)
        self.snp_entry.pack()
        self.snp_predict_button = tk.Button(root, text="Predict SNP Recombination", command=self.predict_snp_recombination)
        self.snp_predict_button.pack(pady=5)
        
        self.prediction_label = tk.Label(root, text="Prediction: ")
        self.prediction_label.pack()
        
        self.plot_button = tk.Button(root, text="Show Linkage Map", command=self.show_linkage_map)
        self.plot_button.pack(pady=5)
        
        self.heatmap_button = tk.Button(root, text="Show Recombination Heatmap", command=self.show_heatmap)
        self.heatmap_button.pack(pady=5)
    
    def upload_snp_file(self):
        file_path = filedialog.askopenfilename(filetypes=[("CSV Files", "*.csv")])
        if file_path:
            self.snp_file_path = file_path
            self.snps = read_snp_csv(file_path)
            messagebox.showinfo("SNP Data", f"SNP data loaded from {file_path}")
            self.snp_listbox.delete(0, tk.END)
            for snp in self.snps:
                self.snp_listbox.insert(tk.END, f"Chrom: {snp['chrom']}, Pos: {snp['position']}")
    
    def predict_snp_recombination(self):
        position = int(self.snp_entry.get())
        result = predict_snp_recombination(position)
        self.prediction_label.config(text=f"Prediction: {result[0]}")
    
    def show_linkage_map(self):
        plot_linkage_map(self.snps)
    
    def show_heatmap(self):
        plot_recombination_heatmap(self.snps)

if __name__ == "__main__":
    root = tk.Tk()
    app = GeneticSimulatorGUI(root)
    root.mainloop()

