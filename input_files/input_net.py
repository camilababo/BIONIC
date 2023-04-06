import itertools

import pandas as pd
from itertools import combinations
from scipy.stats import pearsonr
from tqdm import tqdm


def data_to_edgelist(expression_data_file, edgelist_file):
    """
    This function takes the expression data and creates an edgelist file for the network.

    :param expression_data_file: The expression data file
    :param edgelist_file: The edgelist file
    :return: an edgelist file
    """

    # Read the expression data
    df = pd.read_csv(expression_data_file, sep=',', index_col=0, header=0)

    # Get a list of all the genes
    all_genes = list(df.T.columns)

    # Remove genes of chloroplast and mitochondria organs (ATC and ATM)
    not_organ_genes = [x for x in all_genes if not x.startswith('ATC') and not x.startswith('ATM')]
    expression_data_df = df.loc[not_organ_genes]
    filtered_genes = list(expression_data_df.T.columns)

    # Create a list of edges and their weights
    edges = []
    num_edges = len(list(combinations(filtered_genes, 2)))
    with tqdm(total=num_edges) as pbar:
        for i, (src, trg) in enumerate(combinations(filtered_genes, 2)):
            src_values = expression_data_df.loc[src].values
            trg_values = expression_data_df.loc[trg].values

            weight = 1 - abs(pearsonr(src_values, trg_values)[0])

            # Add the edge and weight to the list
            edges.append((src, trg, weight))

            # Update the progress bar every 1000 edges
            if i % 1000 == 0:
                pbar.update(1000)

        # Write the list of edges to the output file
        with open(edgelist_file, 'w') as f:
            for src, trg, weight in edges:
                f.write(f"{src} {trg} {weight}\n")

    # Update the progress bar to 100% after the loop has finished
    pbar.update(num_edges - pbar.n)


if __name__ == '__main__':
    data_to_edgelist('../outputs/reduced_output_corrected.csv', 'edgelist_corrected.txt')
