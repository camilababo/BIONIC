import numpy as np
from itertools import combinations
from multiprocessing import Pool

from scipy.stats import pearsonr


def compute_edge_weight(args):
    src, trg, df = args
    src_values = df.loc[src].values
    trg_values = df.loc[trg].values
    weight = 1 - abs(pearsonr(src_values, trg_values)[0])
    return f"{src} {trg} {weight}\n"


if __name__ == '__main__':
    # Load data into a numpy array
    data = np.loadtxt('data.txt', skiprows=1)
    genes = data[:, 0]
    values = data[:, 1:]

    # Compute correlations in parallel
    with Pool() as pool:
        args = [(src, trg, values) for src, trg in combinations(genes, 2)]
        results = pool.map(compute_edge_weight, args)

    # Write output to file
    with open('corrected_edgelist.txt', 'w') as f:
        f.writelines(results)
