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


def individual_index(individual_codes, code):
    return np.where(individual_codes == code)[0][0]


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
    

            
# An example of calculating the Breed Conformity Assessment for an individual dog using Microsatellite data
similarity_method = 0 # 0 - PSA; 1 = Jaccard
similarity_matrix = calc_similarity_matrix(genotypes, similarity_method)
individual = np.copy(genotypes[0])
correspondence_percentage = breed_microsatellite(individual, genotypes, similarity_matrix, similarity_method)
print(f"Correspondence percentage is {correspondence_percentage:.2f}%")


# Scripts for Breed Conformity Assessment Using SNP Data

# Function to perform one-hot encoding for a given genotype
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
    

similarity_matrix = cosine_similarity_matrix(vectors)
print('Individual codes of all dogs in the dataset:', individuals)


# An example of calculating the Breed Conformity Assessment for an individual dog using SNP data
individual_code = 'T26'
individual_id = individual_index(individuals, individual_code)
correspondence_percentage = breed_snp(individual_id, similarity_matrix)
print(f"Correspondence percentage is {correspondence_percentage:.2f}%")


# Calculating the Breed Conformity Assessment for all dogs in the dataset using SNP data
breed_snp_cosine = []
for individual_code in individuals:
    individual_id = individual_index(individuals, individual_code)
    correspondence_percentage = breed_snp(individual_id, similarity_matrix)
    breed_snp_cosine.append(correspondence_percentage)




# Compilation and saving of results into a single file

snp_individuals = ['T17','T18','T26','T50','T52','T57','T74','T59','T77','T98','T103','T105','T107','T114','T118','T130','T142','T164','T173','T175','T176','T177','T64','T70','T81','T84','T85','T86','T90','T109','T112','T120','T138','T150','T153','T155','T167','T168','T178']
snp_individual_indices = []
for individual_code in snp_individuals:
    snp_individual_indices.append(individual_index(individuals, individual_code))
snp_individual_indices = np.array(snp_individual_indices)

genotypes223 = np.copy(genotypes)
genotypes39 = np.copy(genotypes[snp_individual_indices])

sim_mat_psa_223 = calc_similarity_matrix(genotypes223, 0)
sim_mat_psa_39 = calc_similarity_matrix(genotypes39, 0)
sim_mat_jaccard_223 = calc_similarity_matrix(genotypes223, 1)
sim_mat_jaccard_39 = calc_similarity_matrix(genotypes39, 1)

microsatellite_psa_223 = []
microsatellite_psa_39 = []
microsatellite_jaccard_223 = []
microsatellite_jaccard_39 = []

for individual_code in snp_individuals:
    individual_id = individual_index(individuals, individual_code)
    individual = np.copy(genotypes[individual_id])
    
    microsatellite_psa_223.append(breed_microsatellite(individual, genotypes223, sim_mat_psa_223, 0))
    microsatellite_jaccard_223.append(breed_microsatellite(individual, genotypes223, sim_mat_jaccard_223, 1))
    
    microsatellite_psa_39.append(breed_microsatellite(individual, genotypes39, sim_mat_psa_39, 0))
    microsatellite_jaccard_39.append(breed_microsatellite(individual, genotypes39, sim_mat_jaccard_39, 1))
    
    
rows = zip(snp_individuals, microsatellite_psa_39, microsatellite_psa_223, microsatellite_jaccard_39, microsatellite_jaccard_223, breed_snp_cosine)
csv_filename = 'breed_conformity_assessment_results.csv'
header = ['Individuals', 'Microsatellite PSA (39)', 'Microsatellite PSA (223)', 'Microsatellite Jaccard (39)', 'Microsatellite Jaccard (223)', 'SNP Cosine']
with open(csv_filename, 'w', newline='') as csvfile:
    csv_writer = csv.writer(csvfile, delimiter=';')
    csv_writer.writerow(header)    
    for row in rows:
        csv_writer.writerow(row)
print(f'The data has been saved to {csv_filename}')