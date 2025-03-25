GeneticLinkageSimulator
Bioinformatics / Machine Learning / Data Visualization
Technologies: Python, Tkinter, NumPy, Matplotlib, Seaborn, scikit-learn, SciPy
Author: Summer (Simran) Malik

About:
The Genetic Linkage Simulator is Python application that allows users to visualize SNP data to explore genetic linkage and where recombination may occur - these are essential concepts in genetic and genomic research.

This uses:
Statistical Analysis
Machine Learning (Random Forest)
Interactive Data Visualization
Tkinter (GUI)

Background (Genetics):
What are SNPS? Why do they Matter?

If you are familiar with programming but not genetics, here is a clear way to break it down:
DNA is like a source code of living organisms, it contains the instructions for everything your body does. These instructions are made of four characters A, T, C, G (kind of like quanternary). An SNP (single nucleotide polymorphism) is a single-letter difference in DNA sequences such as -
Normal DNA: ---AAGTC---
DNA sequence with SNPS: --AAG**G*C--

SNPs can affect how a gene works, or how traits get passed on, and are often used to:
Understand inherited diseases, Trace ancestry, Analyze genetic diversity

So now imagine you have a huge instruction manual that tells your body how to grow, function, and look... that’s your DNA. Inside this manual are chapters called chromosomes, and each chapter has lots of small "spelling variations" called SNPs (pronounced "snips"), or Single Nucleotide Polymorphisms. Some of these SNPs help determine things like, whether you have blue or brown eyes, how your body responds to medicine, what diseases you might be more likely to get...

But here’s the tricky part; when parents pass down their DNA to children, the "manual" isn’t copied perfectly.
When cells divide to create sperm or egg cells, they don’t copy DNA 1:1 — instead, parts of the DNA swap places. This is called genetic recombination (or crossover).

If two SNPs are physically close together on the chromosome, they are likely to be inherited together... this is called genetic linkage.

This project simulates and analyzes those relationships between SNPs, helping to:
Predict where crossovers happen
Visualize SNP inheritance patterns
Identify SNPs that tend to be passed down together

Techincal Overview:

Language Chosen: Python. Python was the ideal choice for this project because of its scientific libraries like NumPy, SciPy, scikit-learn, its data visualization capabillities with Matplotlib and Seaborn, easy GUI usage (Tkinter), and Bioinformatics Libraries like Pysam which was initially considered.

SNP Data Input (CSV parsing): The user uploads a csv data file of any size with the following columns
Chrom: chromosome number or label
Position: genomic position (int)
ref: reference allele
alt: alternate allele

The file is read using Python’s built-in csv module, and SNPs are stored as dictionaries, structured for easy access and grouped by chromosome. This method was used because of its simple and univeral format as well as making it easier for students or researchers to generate and test mock data. This also avoided complex genomic formats due to compatibility issues (see below)

Initial Attempt:
Before using CSV files, I intended to use PySAM, a Python wrapper for reading real genomic data (like BAM/VCF files). However, Miniconda produced persistent dependency conflicts. Even after installing Windows Subsystem for Linux (WSL), errors persisted due to PySAM's reliance on native C libraries and system-level dependencies. These barriers made the project less accessible to non-experts, so I chose to simulate SNP input using CSVs and simplified trait labeling, enabling focus on the core logic and visualization.

Statistical Analysis:
The Chi-Square Test (from scipy.stats): Used to test linkage independence. This compares observed crossover counts with expected counts under the assumption of independence. If SNPs are linked, their recombination patterns will deviate from randomness.

Pearson Correlation Coefficient: Simulates a relationship between SNP positions and binary phenotypes. Random labels are generated for demonstration. Calculates correlation using pearsonr from scipy.stats.

Haplotype Block Detection: Groups SNPs by chromosome and spatial proximity. A basic detection method that mimics how certain SNPs tend to be inherited together. Output is grouped SNPs per chromosome.

Genetic Crossover Simulation:
Simulates crossover likelihood between gene pairs using the Haldane mapping function:
P(recombination) = 1 - exp(-2 * distance / 100)
Where distance is in centiMorgans (cM). Higher distances -----> more likely recombination. Used in biological linkage studies to estimate recombination frequency from physical distance.

Machine Learning Model (Random Forest Classifier):
Library: scikit-learn
Input Features: SNP positions
Labels: Randomly assigned (simulated binary trait)
Split: train_test_split (80/20)
Model: RandomForestClassifier()

Although this doesn't use real trait-SNP associations (due to lack of real phenotype data), it demonstrates how machine learning can be trained on genetic data and the ability to predict recombination events based on position.

Why Random Forest?
Handles non-linearity well
Fast training and interpretable output
Robust to noise (ideal for simulated data)

Visualization:
Linkage Map - Matplotlib. Chromosomes displayed as color-coded horizontal lines and each SNP as a tick mark positioned by its position value.
Recombination Map - Seaborn. Distance based between all SNP Pairs.

Code Structure Summary:

pgsql
Copy
Edit
genetic_simulator.py  
├── Chromosome class: Handles gene-position mapping  
├── CSV parsing: Reads and structures SNP data  
├── Statistical analysis functions  
├── Crossover simulation  
├── Random Forest training + prediction  
├── Linkage map and heatmap visualization  
└── GUI class (Tkinter-based)
Programming Language: Python — Rich libraries, prototyping speed, accessibility
Data Format: CSV — Universally readable, avoids complex parsing
ML Algorithm: Random Forest — Non-linear, robust, easy to interpret
GUI Toolkit: Tkinter — Built-in, lightweight, sufficient for functionality
Visualization: Matplotlib & Seaborn — Standard for scientific plots
Stats Libraries: SciPy — Reliable, proven, and widely used in biology
PySAM: Rejected due to install issues — Maintains accessibility and portability







