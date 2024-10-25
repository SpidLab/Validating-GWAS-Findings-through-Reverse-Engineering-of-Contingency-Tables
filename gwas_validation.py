import numpy as np
from scipy.stats import chi2_contingency

def convert_to_contingency_tables(data):
    tables = []
    
    # Iterate over all SNPs (columns)
    for snp_id in data.columns:
        count = np.bincount(data[snp_id].values, minlength=3)
        
        # Skip SNPs where counts are all zero
        if np.sum(count) == 0:
            continue
        
        tables.append(count)
    
    return tables

def filter_zero_columns(target_count, ref_count):
    combined_count = np.array([target_count, ref_count])
    non_zero_columns = np.any(combined_count != 0, axis=0)
    filtered_target_count = np.array(target_count)[non_zero_columns]
    filtered_ref_count = np.array(ref_count)[non_zero_columns]
    return filtered_target_count, filtered_ref_count

def calc_chi_pvalue(target_count, ref_count):
    filtered_target_count, filtered_ref_count = filter_zero_columns(target_count, ref_count)
    if len(filtered_target_count) < 2 or len(filtered_ref_count) < 2:
        return 1.0
    table = np.array([filtered_target_count, filtered_ref_count])
    return chi2_contingency(table)[1]

def calculate_chi_square_p_values(target_tables, reference_tables):
    p_values = []
    
    # Iterate over all contingency tables
    for target_count, ref_count in zip(target_tables, reference_tables):
        
        p_value = calc_chi_pvalue(target_count, ref_count)
        p_values.append(p_value)
    
    return np.array(p_values)

def get_significant_snps(target_tables, reference_tables, alpha=0.05):
    # Calculate p-values for all SNPs
    p_values = calculate_chi_square_p_values(target_tables, reference_tables)
    
    significant_snps = {}
    
    # Iterate over p-values and filter significant SNPs
    for snp_idx, p_value in enumerate(p_values):
        if p_value < alpha:
            significant_snps[snp_idx] = p_value
    
    return significant_snps

def find_best_case_distribution(n_samples, control_row, target_p_value=0.05):
    best_diff = np.inf
    best_distribution = None
    best_p_value = None
    
    # Iterate through all possible value distributions for the case group
    for s0 in range(n_samples + 1):
        for s1 in range(n_samples - s0 + 1):
            s2 = n_samples - s0 - s1
            case_row = [s0, s1, s2]
            table = np.array([case_row, control_row])
            p_value = calc_chi_pvalue(case_row, control_row)
            p_value_diff = abs(p_value - target_p_value)
            
            if p_value_diff < best_diff:
                best_diff = p_value_diff
                best_distribution = case_row
                best_p_value = p_value
                
    return best_distribution

def calculate_maf_from_contingency_table(table):
    total_alleles = 2 * np.sum(table)
    if total_alleles == 0:
        return 0
    
    # Calculate minor allele frequency
    maf = (table[1] + 2 * table[2]) / total_alleles
    return min(maf, 1 - maf)

def calculate_average_maf_difference(target_mafs, p_values_reporting, control_tables):
    avg_diff_list = []
    
    # Iterate over each reported SNP
    for snp_idx, p_value in p_values_reporting.items():
        control_row = control_tables[snp_idx]
        n_samples = sum(control_row) 
        
        # Find the best matching contingency table for the reported SNP
        best_table = find_best_case_distribution(n_samples, tuple(control_row), target_p_value=p_value)
        
        # Calculate MAF for the best match table
        case_maf = calculate_maf_from_contingency_table(best_table)
        
        # Get the target MAF for the SNP
        target_maf = target_mafs[snp_idx]
        
        # Calculate the difference in MAFs
        maf_diff = abs(case_maf - target_maf)
        avg_diff_list.append(maf_diff)
    
    # Calculate the average difference
    avg_diff = np.mean(avg_diff_list)
    return avg_diff

def evaluate_gwas_reliability(average_maf_difference, threshold=0.1):
    return average_maf_difference <= threshold



## Researcher conducts GWAS
# target_data: The genomic dataset used for GWAS (columns represent SNPs, rows represent individuals)
# reference_data: Public reference dataset with the same SNPs as in target_data

target_data = load_genomic_data()  # Load target genomic dataset
reference_data = load_reference_data()  # Load public reference dataset

# Convert genomic datasets into contingency tables for each SNP
target_contingency_tables = convert_to_contingency_tables(target_data)
reference_contingency_tables = convert_to_contingency_tables(reference_data)

# Identify significant SNPs with p-value < 0.05 in GWAS outcomes
significant_snps = get_significant_snps(target_contingency_tables, reference_contingency_tables, alpha=0.05)




## Verifier validates GWAS outcomes
# target_mafs: MAFs for SNPs included in the GWAS report from public phenotype-specific datasets
target_mafs = get_public_mafs_for_snps(significant_snps)  # Retrieve MAFs for significant SNPs

# Calculate average MAF difference (Hamming distance) between target and reference MAFs
average_maf_difference = calculate_average_maf_difference(target_mafs, significant_snps, reference_contingency_tables)

# Verifier determines reliability of GWAS outcomes based on average Hamming distance
is_reliable = evaluate_gwas_reliability(average_maf_difference, threshold=0.1)  # Adjust threshold as needed
