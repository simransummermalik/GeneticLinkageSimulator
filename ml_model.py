import numpy as np
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
from snp_data import read_vcf

def train_ml_with_snps(vcf_file):
 
    snps = read_vcf(vcf_file, "1")
    data = np.array([[snp["position"]] for snp in snps])  
    labels = np.random.choice([0, 1], size=len(snps))  # Random binary labels (low/high recombination)

    X_train, X_test, y_train, y_test = train_test_split(data, labels, test_size=0.2)
    model = RandomForestClassifier()
    model.fit(X_train, y_train)
    
    return model
ml_model = train_ml_with_snps("example.csv")

def predict_snp_recombination(position):
    return ml_model.predict([[position]])
