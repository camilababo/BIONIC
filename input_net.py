import itertools

import pandas as pd
from itertools import combinations
from scipy.stats import pearsonr


def data_to_edgelist(expression_data):
    """
    This function takes the expression data and creates an edgelist file for the network.

    :param expression_data: The expression data file
    :return: an edgelist file
    """

    # Read the expression data
    df = pd.read_csv(expression_data, sep=',', index_col=0, header=0)

    # Get a list of all the genes
    genes = list(df.T.columns)

    # Write an edge file where edges are added to the file one at a time to avoid memory errors
    with open('edgelist.txt', 'w') as f:
        for src, trg in combinations(genes, 2):
            src_values = df.loc[src].values
            trg_values = df.loc[trg].values

            weight, p_value = pearsonr(src_values, trg_values)

            # Write the edge and weight to the output file
            f.write(f"{src} {trg} {weight}\n")


if __name__ == '__main__':
    data_to_edgelist('reduced_output.csv')
