# MultipleComparisons

# About
MultipleComparisons allows to calculate and process extended (e.g., n-ary) similarity indices.
The basic theory is detailed in: "Extended similarity indices: the benefits of comparing more than two objects simultaneously. Part 1: Theory and characteristics", R. A. Miranda-Quintana, D. Bajusz, A. Rácz, K. Héberger; J. Cheminformatics https://jcheminf.biomedcentral.com/articles/10.1186/s13321-021-00505-3

Some applications are presented in: "Extended similarity indices: the benefits of comparing more than two objects simultaneously. Part 2: speed, consistency, diversity selection", R. A. Miranda-Quintana, A. Rácz, D. Bajusz, K. Héberger; J. Cheminformatics https://jcheminf.biomedcentral.com/articles/10.1186/s13321-021-00504-4

# License
MultipleComparisons is distributed under GPL License version 3 (GPLv3).

# Dependencies
Python >= 3.3;  http://www.python.org/

Numpy >= 1.9.1;  http://www.numpy.org/

SciPy >= 0.11.0;  http://www.scipy.org/

Matplotlib >= 1.0;  http://matplotlib.org/

# Utils Module Documentation

The `utils` module contains a suite of functions designed to facilitate complex data processing tasks, specifically geared towards handling, calculating, and saving datasets related to similarity matrices, selection strategies, and data serialization.

## Module Functions

### 1. **`calculate_counters(data_sets, c_threshold=None, w_factor="fraction")`**
Calculates various similarity and dissimilarity counters based on provided datasets.

#### Parameters:
- **`data_sets`**: An array of arrays, each containing matrix column sums and the number of fingerprints.
- **`c_threshold`**: Defines the coincidence threshold (None, 'dissimilar', or an integer).
- **`w_factor`**: Defines the weight function used ('fraction', 'power_n', or other).

#### Returns:
- A dictionary with weighted and non-weighted counters.

#### Usage:
Ideal for calculating detailed similarity metrics which can be adjusted by thresholds and weighting factors.

### 2. **`calculate_medoid(total_data, n_ary='RR', weight='nw')`**
Calculates the medoid of a dataset based on a specific similarity type and weight.

#### Parameters:
- **`total_data`**: The dataset from which to calculate the medoid.
- **`n_ary`**: Type of n-ary index to use.
- **`weight`**: Weight type to use.

#### Usage:
Useful for identifying a central or most typical element in a dataset.

### 3. **`gen_sim_dict(data_sets, c_threshold=None, w_factor="fraction")`**
Generates a dictionary of similarity indices based on the specified parameters.

#### Usage:
This function is a higher-level abstraction used to compile various similarity indices into a single dictionary, facilitating easy access and manipulation of these indices.

### 4. **`calculate_outlier(total_data, n_ary='RR', weight='nw')`**
Identifies the outlier of a dataset using specified similarity measures.

#### Usage:
Effective in detecting anomalies or outliers in a dataset which can be pivotal for certain analyses.

### 5. **`get_single_index(...)`**
Applies a binary tie-breaker selection criterion to select a unique index from a dataset.

### 6. **`get_new_index_n(...)`**
Selects a diverse object using the ECS_MeDiv algorithm based on current dataset conditions.

### 7. **`load_data(file_path)`**
Loads data from a specified file path, handling different file formats (`.parquet`, `.pkl`).

### 8. **`save_data(data, n_ary, c_threshold, max_n, base_name, format='pkl', directory='output')`**
Saves data to a file in the specified format. Handles the creation of new files or appending to existing files based on the data size and specified maximum entries (`max_n`).

#### Parameters:
- **`data`**: Data to be saved.
- **`n_ary`, `c_threshold`, `max_n`, `base_name`**: Parameters that define the file naming convention.
- **`format`**: The format to save the data in (currently supports 'pkl' for pickle and 'parquet').
- **`directory`**: The directory to save the file in.

#### Usage:
Enables flexible data storage with an intelligent file management system that prevents data overwrites and maintains version control through indexing.

### 9. **`extract_dataset(files, n_arys, c_thresholds, max_n, start, output_format='pkl')`**
Main function to process and save datasets based on user-defined specifications.

#### Usage:
Central to executing the dataset extraction and processing workflow, allowing for customization across various parameters and easy scalability.

## Getting Started
To use the `utils` module, ensure that you have all dependencies installed, including `numpy`, `pandas`, and `pyarrow` for handling Parquet files. Import the functions you need from the module into your main script and follow the usage examples as guidelines for implementing functionality in your projects.

## Contribution
Contributions to enhance the module are welcome. Please ensure to maintain the existing coding style, add appropriate comments, and update the README.md file as necessary.
"""

# Other Functionality
1. ECS_MeDiv: linearly-scaling extended similarity-based diversity selection algorithm with binary similarity tie breaker criterion. 

# Reference
Please, cite both the associated manuscripts:

"Extended similarity indices: the benefits of comparing more than two objects simultaneously. Part 1: Theory and characteristics", R. A. Miranda-Quintana, D. Bajusz, A. Rácz, K. Héberger; J. Cheminformatics 13 (32), 2021, https://doi.org/10.1186/s13321-021-00505-3

"Extended similarity indices: the benefits of comparing more than two objects simultaneously. Part 2: speed, consistency, diversity selection", R. A. Miranda-Quintana, A. Rácz, D. Bajusz, K. Héberger; J. Cheminformatics 13 (33), 2021, https://doi.org/10.1186/s13321-021-00504-4


And this repository:

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3997606.svg)](https://doi.org/10.5281/zenodo.3997606)

