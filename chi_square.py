from scipy.stats import chi2_contingency

def chi_square_linkage(observed_counts):
    """Performs Chi-Square test for genetic linkage."""
    total = sum(observed_counts) 
    expected = [total / len(observed_counts)] * len(observed_counts)
    chi2, p_value = chi2_contingency([observed_counts, expected])[:2]
    
    return p_value < 0.05 