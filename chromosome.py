class Chromosome:
    def __init__(self, genes, positions):

        self.genes = genes
        self.positions = positions  

    def display(self):

        print("\nChromosome Gene Map:")
        for gene, pos in zip(self.genes, self.positions):
            print(f"Gene {gene}: {pos} cM")