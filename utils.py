import numpy as np
from math import ceil
from math import log
import pickle
import os
import pandas as pd
from collections import defaultdict

def calculate_counters(data_sets, c_threshold=None, w_factor="fraction"):
    """Calculate 1-similarity, 0-similarity, and dissimilarity counters

    Arguments
    ---------
    data_sets : np.ndarray
        Array of arrays. Each sub-array contains m + 1 elements,
        with m being the length of the fingerprints. The first
        m elements are the column sums of the matrix of fingerprints.
        The last element is the number of fingerprints.

    c_threshold : {None, 'dissimilar', int}
        Coincidence threshold.
        None : Default, c_threshold = n_fingerprints % 2
        'dissimilar' : c_threshold = ceil(n_fingerprints / 2)
        int : Integer number < n_fingerprints

    w_factor : {"fraction", "power_n"}
        Type of weight function that will be used.
        'fraction' : similarity = d[k]/n
                     dissimilarity = 1 - (d[k] - n_fingerprints % 2)/n_fingerprints
        'power_n' : similarity = n**-(n_fingerprints - d[k])
                    dissimilarity = n**-(d[k] - n_fingerprints % 2)
        other values : similarity = dissimilarity = 1

    Returns
    -------
    counters : dict
        Dictionary with the weighted and non-weighted counters.

    Notes
    -----
    Please, cite the original papers on the n-ary indices:
    https://jcheminf.biomedcentral.com/articles/10.1186/s13321-021-00505-3
    https://jcheminf.biomedcentral.com/articles/10.1186/s13321-021-00504-4
    """
    # Setting matches
    total_data = np.sum(data_sets, axis=0)
    n_fingerprints = int(total_data[-1])
    c_total = total_data[:-1]
    
    # Assign c_threshold
    if not c_threshold:
        c_threshold = n_fingerprints % 2
    if isinstance(c_threshold, str):
        if c_threshold != 'dissimilar':
            raise TypeError("c_threshold must be None, 'dissimilar', or an integer.")
        else:
            c_threshold = ceil(n_fingerprints / 2)
    if isinstance(c_threshold, int):
        if c_threshold >= n_fingerprints:
            raise ValueError("c_threshold cannot be equal or greater than n_fingerprints.")
        c_threshold = c_threshold
    
    # Set w_factor
    if w_factor:
        if "power" in w_factor:
            power = int(w_factor.split("_")[-1])
            def f_s(d):
                return power**-float(n_fingerprints - d)
    
            def f_d(d):
                return power**-float(d - n_fingerprints % 2)
        elif w_factor == "fraction":
            def f_s(d):
                return d/n_fingerprints
    
            def f_d(d):
                return 1 - (d - n_fingerprints % 2)/n_fingerprints
        else:
            def f_s(d):
                return 1
    
            def f_d(d):
                return 1
    else:
        def f_s(d):
            return 1
    
        def f_d(d):
            return 1
    
    # Calculate a, d, b + c
    a = 0
    w_a = 0
    d = 0
    w_d = 0
    total_dis = 0
    total_w_dis = 0
    for s in c_total:
        if 2 * s - n_fingerprints > c_threshold:
            a += 1
            w_a += f_s(2 * s - n_fingerprints)
        elif n_fingerprints - 2 * s > c_threshold:
            d += 1
            w_d += f_s(abs(2 * s - n_fingerprints))
        else:
            total_dis += 1
            total_w_dis += f_d(abs(2 * s - n_fingerprints))
    total_sim = a + d
    total_w_sim = w_a + w_d
    p = total_sim + total_dis
    w_p = total_w_sim + total_w_dis
    
    counters = {"a": a, "w_a": w_a, "d": d, "w_d": w_d,
                "total_sim": total_sim, "total_w_sim": total_w_sim,
                "total_dis": total_dis, "total_w_dis": total_w_dis,
                "p": p, "w_p": w_p}
    
    return counters

def calculate_medoid(total_data, n_ary = 'RR', weight = 'nw'):
    """Calculate the medoid of a set"""
    index = len(total_data[0]) + 1
    min_sim = 3.08
    total_sum = np.sum(total_data, axis = 0)
    for i, pixel in enumerate(total_data):
        i_sum = total_sum - total_data[i]
        data_sets = [np.append(i_sum, len(total_data) - 1)]
        Indices = gen_sim_dict(data_sets)
        sim_index = Indices[weight][n_ary]
        if sim_index < min_sim:
            min_sim = sim_index
            index = i
        else:
            pass
    return index

def gen_sim_dict(data_sets, c_threshold=None, w_factor="fraction"):
    counters = calculate_counters(data_sets, c_threshold=c_threshold, w_factor=w_factor)
    # Indices
    # AC: Austin-Colwell, BUB: Baroni-Urbani-Buser, CTn: Consoni-Todschini n
    # Fai: Faith, Gle: Gleason, Ja: Jaccard, Ja0: Jaccard 0-variant
    # JT: Jaccard-Tanimoto, RT: Rogers-Tanimoto, RR: Russel-Rao
    # SM: Sokal-Michener, SSn: Sokal-Sneath n
    epsilon = 1e-10

    # Weighted Indices
    ac_w = (2/np.pi) * np.arcsin(np.sqrt(counters['total_w_sim']/
                                         counters['w_p']))
    bub_w = ((counters['w_a'] * counters['w_d'])**0.5 + counters['w_a'])/\
            ((counters['w_a'] * counters['w_d'])**0.5 + counters['w_a'] + counters['total_w_dis'])
    ct1_w = (log(1 + counters['w_a'] + counters['w_d']))/\
            (log(1 + counters['w_p']))
    ct2_w = (log(1 + counters['w_p']) - log(1 + counters['total_w_dis']))/\
            (log(1 + counters['w_p']))
    ct3_w = (log(1 + counters['w_a']))/\
            (log(1 + counters['w_p']))
    ct4_w = (log(1 + counters['w_a']))/\
            (log(1 + counters['w_a'] + counters['total_w_dis'])+ epsilon)
    fai_w = (counters['w_a'] + 0.5 * counters['w_d'])/\
            (counters['w_p'])
    gle_w = (2 * counters['w_a'])/\
            (2 * counters['w_a'] + counters['total_w_dis'] + epsilon)
    ja_w = (3 * counters['w_a'])/\
           (3 * counters['w_a'] + counters['total_w_dis'] + epsilon)
    ja0_w = (3 * counters['total_w_sim'])/\
            (3 * counters['total_w_sim'] + counters['total_w_dis'] + epsilon)
    jt_w = (counters['w_a'])/\
           (counters['w_a'] + counters['total_w_dis'] + epsilon)
    rt_w = (counters['total_w_sim'])/\
           (counters['w_p'] + counters['total_w_dis'] + epsilon)
    rr_w = (counters['w_a'])/\
           (counters['w_p'])
    sm_w =(counters['total_w_sim'])/\
          (counters['w_p'])
    ss1_w = (counters['w_a'])/\
            (counters['w_a'] + 2 * counters['total_w_dis'] + epsilon)
    ss2_w = (2 * counters['total_w_sim'])/\
            (counters['w_p'] + counters['total_w_sim'] + epsilon)


    ## Non-Weighted Indices
    ac_nw = (2/np.pi) * np.arcsin(np.sqrt(counters['total_w_sim']/
                                          counters['p']))
    bub_nw = ((counters['w_a'] * counters['w_d'])**0.5 + counters['w_a'])/\
             ((counters['a'] * counters['d'])**0.5 + counters['a'] + counters['total_dis'])
    ct1_nw = (log(1 + counters['w_a'] + counters['w_d']))/\
             (log(1 + counters['p']))
    ct2_nw = (log(1 + counters['w_p']) - log(1 + counters['total_w_dis']))/\
             (log(1 + counters['p']))
    ct3_nw = (log(1 + counters['w_a']))/\
             (log(1 + counters['p']))
    ct4_nw = (log(1 + counters['w_a']))/\
             (log(1 + counters['a'] + counters['total_dis']) + epsilon)
    fai_nw = (counters['w_a'] + 0.5 * counters['w_d'])/\
             (counters['p'])
    gle_nw = (2 * counters['w_a'])/\
             (2 * counters['a'] + counters['total_dis'] + epsilon)
    ja_nw = (3 * counters['w_a'])/\
            (3 * counters['a'] + counters['total_dis'] + epsilon)
    ja0_nw = (3 * counters['total_w_sim'])/\
             (3 * counters['total_sim'] + counters['total_dis'] + epsilon)
    jt_nw = (counters['w_a'])/\
            (counters['a'] + counters['total_dis'] + epsilon)
    rt_nw = (counters['total_w_sim'])/\
            (counters['p'] + counters['total_dis'] + epsilon)
    rr_nw = (counters['w_a'])/\
            (counters['p'])
    sm_nw =(counters['total_w_sim'])/\
           (counters['p'])
    ss1_nw = (counters['w_a'])/\
             (counters['a'] + 2 * counters['total_dis'] + epsilon)
    ss2_nw = (2 * counters['total_w_sim'])/\
             (counters['p'] + counters['total_sim'] + epsilon)

    # Dictionary with all the results
    Indices = {'nw': {'AC':ac_nw,
                      'BUB':bub_nw,
                      'CT1':ct1_nw,
                      'CT2':ct2_nw,
                      'CT3':ct3_nw,
                      'CT4':ct4_nw,
                      'Fai':fai_nw,
                      'Gle':gle_nw,
                      'Ja0':ja0_nw,
                      'Ja':ja_nw,
                      'JT':jt_nw,
                      'RT':rt_nw,
                      'RR':rr_nw,
                      'SM':sm_nw,
                      'SS1':ss1_nw,
                      'SS2':ss2_nw},
                'w': {'AC':ac_w,
                      'BUB':bub_w,
                      'CT1':ct1_w,
                      'CT2':ct2_w,
                      'CT3':ct3_w,
                      'CT4':ct4_w,
                      'Fai':fai_w,
                      'Gle':gle_w,
                      'Ja0':ja0_w,
                      'Ja':ja_w,
                      'JT':jt_w,
                      'RT':rt_w,
                      'RR':rr_w,
                      'SM':sm_w,
                      'SS1':ss1_w,
                      'SS2':ss2_w}}
    return Indices

def calculate_outlier(total_data, n_ary = 'RR', weight = 'nw'):
    """Calculate the outlier of a set"""
    index = len(total_data[0]) + 1
    max_sim = -3.08
    total_sum = np.sum(total_data, axis = 0)
    for i, pixel in enumerate(total_data):
        i_sum = total_sum - total_data[i]
        data_sets = [np.append(i_sum, len(total_data) - 1)]
        Indices = gen_sim_dict(data_sets)
        sim_index = Indices[weight][n_ary]
        if sim_index > max_sim:
            max_sim = sim_index
            index = i
        else:
            pass
    return index

def get_single_index(total_data, indices, selected_n, c_threshold=None, n_ary='RR', weight='nw'):
    """Binary tie-breaker selection criterion"""
    index = -1  # inizializzato a -1 per indicare "non trovato"
    min_value = 3.08
    for i in indices:
        v = 0
        for j in selected_n:
            c_total = total_data[j] + total_data[i]
            data_sets = [np.append(c_total, len(selected_n) + 1)]  # +1 per includere i nell'analisi
            Indices = gen_sim_dict(data_sets, c_threshold=c_threshold)
            sim_index = Indices[weight][n_ary]
            v += sim_index
        av_v = v / (len(selected_n) + 1)
        if av_v < min_value:
            index = i
            min_value = av_v

    if index == -1:
        print("Attenzione: nessun indice valido trovato; potrebbe essere necessario rivedere min_value o i dati di input.")
    return index

def get_new_index_n(total_data, selected_condensed, n, select_from_n, selected_n, c_threshold=None,
                    n_ary = 'RR', weight = 'nw'):
    """Select a diverse object using the ECS_MeDiv algorithm"""
    n_total = n + 1
    # min value that is guaranteed to be higher than all the comparisons
    min_value = 3.08
    
    # placeholder index
    indices = [len(total_data[0]) + 1]
    
    # for all indices that have not been selected
    for i in select_from_n:
        # column sum
        c_total = selected_condensed + total_data[i]
        # calculating similarity
        data_sets = [np.append(c_total, n_total)]
        Indices = gen_sim_dict(data_sets, c_threshold=c_threshold)
        sim_index = Indices[weight][n_ary]
        # if the sim of the set is less than the similarity of the previous diverse set, update min_value and index
        if sim_index < min_value:
            indices = [i]
            min_value = sim_index
        elif sim_index == min_value:
            indices.append(i)
    if len(indices) == 1:
        index = indices[0]
    else:
        # Use average of binary similarities as tie-breaker
        index = get_single_index(total_data, indices, selected_n, c_threshold=None, n_ary = n_ary, weight = 'nw')
    return index

def load_data(file_path):
    # Controlla l'estensione del file per determinare come caricarlo
    extension = os.path.splitext(file_path)[1]
    if extension == '.parquet':
        return pd.read_parquet(file_path)
    elif extension == '.pkl':
        with open(file_path, 'rb') as f:
            return pickle.load(f)
    else:
        raise ValueError(f"Unsupported file format: {extension}")

def save_data(data, n_ary, c_threshold, max_n, base_name, format='pkl', directory='output'):
    if not os.path.exists(directory):
        os.makedirs(directory)

    base_file_name = f"{n_ary}_{c_threshold}_{max_n}_{base_name}"
    file_path = os.path.join(directory, base_file_name + f'.{format}')

    if os.path.exists(file_path):
        # Carica i dati esistenti utilizzando la funzione load_data
        existing_data = load_data(file_path)
        if len(existing_data) >= max_n:
            counter = 1
            while True:
                new_file_path = os.path.join(directory, f"{base_file_name}_{counter}.{format}")
                if not os.path.exists(new_file_path):
                    file_path = new_file_path
                    break
                counter += 1
        else:
            # Combina i dati esistenti con i nuovi dati, se possibile
            if isinstance(data, list) and isinstance(existing_data, list):
                data = existing_data + data
            elif isinstance(data, pd.DataFrame) and isinstance(existing_data, pd.DataFrame):
                data = pd.concat([existing_data, data], ignore_index=True)

    # Salva i dati nel formato appropriato
    if format == 'parquet':
        if not isinstance(data, pd.DataFrame):
            data = pd.DataFrame(data)
        data.to_parquet(file_path, index=False)
    else:
        with open(file_path, 'wb') as f:
            pickle.dump(data, f)
    
    print(f"Data saved to {file_path} in {format} format.")

def extract_dataset(files, n_arys, c_thresholds, max_n, start, output_format='pkl'):
    # Main loop
    DiversityDict = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))
    for n_ary in n_arys:
        for c_threshold in c_thresholds:
            for file in files:
                if os.path.exists(file):
                    base_name = os.path.splitext(os.path.basename(file))[0]
                    extension = os.path.splitext(os.path.basename(file))[1]
                    total_data = np.load(file)
                else:
                    print(f"File not found: {file}")
                    break
                if 'rank' in file:
                    pass
                else:
                    # Gestisci il file in base alla sua estensione
                    if extension == '.npy':
                        # Carica un file numpy
                        data = np.load(file)
                        print("Caricato file .npy:", data)
                        # Implementa ulteriori operazioni specifiche per i file .npy

                    elif extension == '.csv':
                        # Carica un file CSV
                        data = pd.read_csv(file)
                        print("Caricato file .csv:", data)
                        # Implementa ulteriori operazioni specifiche per i file .csv

                    elif extension == '.pkl':
                        # Carica un file pickle
                        with open(file, 'rb') as f:
                            data = pickle.load(f)
                        print("Caricato file .pkl:", data)
                        # Implementa ulteriori operazioni specifiche per i file .pkl

                    else:
                        print("Estensione del file non supportata:", extension)
                        break
                    # total number of fingerprints
                    fp_total = len(total_data)
                    
                    # indices of all the fingerprints
                    total_indices = np.array(range(fp_total))
                    
                    # starting point
                    if start =='medoid':
                        seed = calculate_medoid(total_data, n_ary = n_ary)
                    elif start == 'random':
                        seed = random.randint(0, fp_total - 1)
                    elif start == 'out':
                        seed = calculate_outlier(total_data, n_ary = n_ary)
                    else:
                        print('Select a correct starting point')
                    selected_n = [seed]
                    
                    # vector with the column sums of all the selected fingerprints
                    selected_condensed = total_data[seed]
                    
                    # number of fingerprints selected
                    n = 1
                    while len(selected_n) < max_n:
                        # indices from which to select the new fingerprints
                        select_from_n = np.delete(total_indices, selected_n)
                        
                        # new index selected
                        new_index_n = get_new_index_n(total_data, selected_condensed, n, select_from_n, selected_n,
                                                    c_threshold=c_threshold, n_ary = n_ary)
                        
                        if new_index_n == -1:
                            print("No more valid indices found, stopping.")
                            break  
                        # Exit the while loop if no valid index is found
    
                        # updating column sum vector
                        selected_condensed += total_data[new_index_n]
                        
                        # updating selected indices
                        selected_n.append(new_index_n)
                        
                        # updating n
                        n = len(selected_n)
                        print(selected_n)
                    DiversityDict[n_ary][c_threshold][base_name] = selected_n
                    # file_name = f"{n_ary}_{c_threshold}_{max_n}_{base_name}.pkl"
                    # #Create a pickle file with the selected molecules
                    # with open(file_name, 'wb') as f:
                    #     pickle.dump(selected_n, f)
                    save_data(selected_n, n_ary, c_threshold, max_n, base_name, format=output_format, directory='output')