import csv
import numpy as np
import pandas as pd

# Scripts for Breed Conformity Assessment Using Microsatellite Data


# Calculates genetic similarity between two individuals based on their 
# microsatellite genotypes using the proportion of shared alleles (PSA) method
def psa_similarity(a, b):
    shared = 0
    loci = len(a) // 2    
    for i in range(loci):
        if a[2*i] in b[2*i:2*i+2] or a[2*i+1] in b[2*i:2*i+2]:
            shared += 1            
    return shared / loci


def jaccard_similarity(a, b):
    intersection = len(set(a).intersection(set(b)))
    union = len(set(a).union(set(b)))
    return intersection / union

dataset_file_name = 'dataset_tazy_microsat_223.csv'

with open(dataset_file_name, mode='r', encoding='windows-1252') as csv_file:
    data = [line.strip() for line in csv_file]

genotypes = []
individuals = []

for row in data:
    row_values = row.split(";")
    individuals.append(row_values[0])
    genotypes.append(list(map(int, row_values[1:])))
individuals = np.array(individuals)
genotypes = np.array(genotypes)
    

def calc_similarity_row(individual, genotypes, method):
    n = len(genotypes)
    similarity_row = np.zeros(n)
    for i in range(n):
        if method:
            similarity_row[i] = jaccard_similarity(individual, genotypes[i])
        else:
            similarity_row[i] = psa_similarity(individual, genotypes[i])
    return similarity_row
              
    
def calc_similarity_matrix(genotypes, method):
    n = len(genotypes)
    similarity_matrix = np.full((n, n), 1.0)
    for i in range(n):
        for j in range(i + 1, n):
            if method:
                similarity = jaccard_similarity(genotypes[i], genotypes[j])
            else:
                similarity = psa_similarity(genotypes[i], genotypes[j])
            similarity_matrix[i, j] = similarity
            similarity_matrix[j, i] = similarity
    return similarity_matrix


# similarity_method: 0 - PSA; 1 = Jaccard
def breed_microsatellite(individual, genotypes, similarity_matrix, similarity_method=0):
    max_median_similarity = np.max(np.median(similarity_matrix, axis=1))
    individual_similarity_row = calc_similarity_row(individual, genotypes, similarity_method)
    individual_median_similarity = np.median(individual_similarity_row)
    correspondence_percentage = min(individual_median_similarity/max_median_similarity, 1.0)*100
    return correspondence_percentage


def one_hot_encode(genotype):
    # Define the possible genotypes including the 'missing' category
    possible_genotypes = ['0/0', '0/1', 'missing']
    # Check if the genotype is missing or not in the possible genotypes
    if genotype not in possible_genotypes:
        genotype = 'missing'
    # Create a dictionary for one-hot encoding
    encoding = {g: 1 if g == genotype else 0 for g in possible_genotypes}    
    # Return the encoding as a list
    return list(encoding.values())


# Load the SNP dataset
file_path = 'dataset_tazy_snp_39.xlsm'
data = pd.read_excel(file_path)

# Replace missing genotype value with 'missing'
data.replace(44927, 'missing', inplace=True)
data.replace('./.', 'missing', inplace=True) 

# Initialize a list to collect the encoded vectors for each dog
encoded_vectors_list = []
individuals_list = []

# Iterate over each dog column
for column in data.columns[3:]:  # Skip the first three columns (ID, REF, ALT)
    # Apply one-hot encoding to each genotype
    encoded = data[column].apply(one_hot_encode)
    # Convert the list of dictionaries into a DataFrame and then to a numpy array
    encoded_df = pd.DataFrame(list(encoded))
    encoded_array = encoded_df.to_numpy()
    # Flatten the array to create a single vector for each dog
    encoded_vectors_list.append(encoded_array.flatten())    
    individuals_list.append(column)

# Convert the list of vectors into a 2D NumPy matrix
# vectors is a 2D NumPy matrix where each row represents a dog
# and columns represent the one-hot encoded genotypes across all loci.
individuals = np.array(individuals_list)
vectors = np.array(encoded_vectors_list)


def individual_index(individual_codes, code):
    return np.where(individual_codes == code)[0][0]

def cosine_similarity_matrix(vectors):
    # Normalize each feature vector to unit length
    norm_vectors = vectors / np.linalg.norm(vectors, axis=1, keepdims=True)
    # Calculate the cosine similarity matrix
    return np.dot(norm_vectors, norm_vectors.T)


def breed_snp(individual_id, similarity_matrix):
    max_median_similarity = np.max(np.median(similarity_matrix, axis=1))
    individual_similarity_row = similarity_matrix[individual_id]
    individual_median_similarity = np.median(individual_similarity_row)
    correspondence_percentage = min(individual_median_similarity/max_median_similarity, 1.0)*100
    return correspondence_percentage


def assess_external_individual_microsatellite(alleles, genotypes, similarity_matrix, similarity_method=0):
    """
    Calculate breed correspondence percentage for an external individual 
    using microsatellite genotype data.

    Parameters:
        alleles (list): List of integers representing microsatellite alleles (e.g., [89, 89, 118, 118, ...])
        genotypes (np.ndarray): Genotype matrix from the dataset (n_samples x 2N)
        similarity_matrix (np.ndarray): Precomputed similarity matrix of dataset individuals
        similarity_method (int): 0 = PSA, 1 = Jaccard

    Returns:
        float: Breed correspondence percentage
    """
    individual_similarity_row = calc_similarity_row(alleles, genotypes, similarity_method)
    max_median_similarity = np.max(np.median(similarity_matrix, axis=1))
    individual_median_similarity = np.median(individual_similarity_row)

    if max_median_similarity == 0:
        return 0.0

    correspondence_percentage = min(individual_median_similarity / max_median_similarity, 1.0) * 100
    return correspondence_percentage



def assess_external_individual_snp(genotype_list, vectors, similarity_matrix):
    """
    Calculate breed correspondence percentage for an external individual 
    using SNP genotype data.

    Parameters:
        genotype_list (list): List of genotype strings (e.g., ['0/0', '0/1', 'missing', ...])
        vectors (np.ndarray): Existing one-hot encoded dataset matrix (n_samples x features)
        similarity_matrix (np.ndarray): Precomputed cosine similarity matrix for the dataset

    Returns:
        float: Breed correspondence percentage
    """
    # Reuse one-hot encoding
    encoded = [one_hot_encode(gt) for gt in genotype_list]
    flat_vector = np.array(encoded).flatten()

    # Normalize external individual vector
    norm = np.linalg.norm(flat_vector)
    if norm > 0:
        flat_vector = flat_vector / norm

    # Normalize dataset vectors if not already normalized
    norm_vectors = vectors / np.linalg.norm(vectors, axis=1, keepdims=True)

    # Compute similarity to all dataset individuals
    similarity_row = np.dot(norm_vectors, flat_vector)

    max_median_similarity = np.max(np.median(similarity_matrix, axis=1))
    individual_median_similarity = np.median(similarity_row)

    if max_median_similarity == 0:
        return 0.0

    correspondence_percentage = min(individual_median_similarity / max_median_similarity, 1.0) * 100
    return correspondence_percentage






# -----------------------------------------
# EXAMPLE 1: Microsatellite External Sample
# -----------------------------------------

# This individual is not in the dataset
external_microsat = [
    87, 87, 120, 122, 162, 166, 212, 212, 226, 232,
    99, 101, 131, 131, 208, 216, 250, 250, 284, 286,
    123, 126, 150, 152, 234, 238, 98, 104, 156, 172,
    192, 208, 229, 235, 268, 276
]

# Calculate similarity matrix from the original dataset
sim_matrix_psa = calc_similarity_matrix(genotypes, method=0)  # PSA
sim_matrix_jaccard = calc_similarity_matrix(genotypes, method=1)  # Jaccard

# Compute correspondence for external individual
conformity_psa = assess_external_individual_microsatellite(external_microsat, genotypes, sim_matrix_psa, similarity_method=0)
conformity_jaccard = assess_external_individual_microsatellite(external_microsat, genotypes, sim_matrix_jaccard, similarity_method=1)

print(f"[Microsatellite] PSA Correspondence: {conformity_psa:.2f}%")
print(f"[Microsatellite] Jaccard Correspondence: {conformity_jaccard:.2f}%")

# -----------------------------------------
# EXAMPLE 2: SNP External Sample (Fixed)
# -----------------------------------------

# Determine the number of SNPs in the dataset
num_snps = vectors.shape[1] // 3

# Create a valid SNP vector matching the dataset
# For testing, we'll alternate some genotypes
external_snp = ['0/0'] * (num_snps // 3) + ['0/1'] * (num_snps // 3) + ['missing'] * (num_snps - 2 * (num_snps // 3))

# Check length match
assert len(external_snp) == num_snps, "Mismatch in SNP vector length"

# Compute correspondence
sim_matrix_cosine = cosine_similarity_matrix(vectors)
conformity_snp = assess_external_individual_snp(external_snp, vectors, sim_matrix_cosine)

print(f"[SNP] Cosine Correspondence: {conformity_snp:.2f}%")

