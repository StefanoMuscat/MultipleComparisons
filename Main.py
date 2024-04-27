import glob
from utils import extract_dataset

####################### INPUTS #################
# Here you can specify the types of similarity matrices to be processed.
# Uncomment the desired options or add new ones as needed.
n_arys = [
    'AC', # Uncomment for Angular Coefficient-based similarity
    #'BUB', # Uncomment for Bubble-based similarity
    #'CT1', # Add as many as required
    #'CT2',
    #'CT3',
    #'CT4',
    #'Fai',
    #'Gle',
    #'Ja',
    #'Ja0',
    #'JT',
    #'RT',
    #'RR',
    #'SM',
    #'SS1',
    #'SS2'
]

# Set the threshold types for dissimilarity. Extend this list as required.
c_thresholds = [
    'dissimilar'  # Currently only processing 'dissimilar' threshold. Add more as required.
]

# Specify the maximum number of molecules to be selected from the dataset.
max_n = 5

# Seed selection strategy:
# - 'medoid': Starts selection from the medoid of the dataset.
# - 'random': Begins selection at a random molecule.
# - 'out': Commences selection from an outlier.
# Choose one according to the specific requirements of the dataset.
start = 'medoid'

# Specify the path to the input files. Ensure the path is correct.
# Add multiple paths if necessary by extending the list.
input_files = [
    "ECS_MeDiv/sample_calculations/normalized_CYP_SH2_data/CYP_complex_MD_10_final_formatted_uniform.npy"
]

####################### MAIN #################
# Extracts and processes the dataset based on the inputs specified above.
# The function 'extract_dataset' is defined in 'utils.py'.
# It processes files based on the specified parameters and saves the result in the specified output format.
# Supported formats include 'pkl' for Pickle files. Add other formats as required.
extract_dataset(input_files, n_arys, c_thresholds, max_n, start, output_format='pkl')
