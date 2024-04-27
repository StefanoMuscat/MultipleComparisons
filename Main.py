import glob
from utils import extract_dataset

####################### INPUTS #################
# List of type of similarity matrices
n_arys = ['AC']# 'BUB' , 'CT1', 'CT2', 'CT3', 'CT4', 'Fai',
        #   'Gle', 'Ja', 'Ja0', 'JT', 'RT', 'RR', 'SM', 'SS1', 'SS2']

# Dissimilarity thresholds
c_thresholds = ['dissimilar']

# Number of molecules to select
max_n = 5

# Seed selection:
# medoid = start from medoid
# random = select random initial seed
# out = start from outlier
start = 'medoid'
input_files = ["ECS_MeDiv/sample_calculations/normalized_CYP_SH2_data/CYP_complex_MD_10_final_formatted_uniform.npy"]

####################### MAIN #################
# Extract dataset
extract_dataset(input_files, n_arys, c_thresholds, max_n, start, output_format='pkl')