from urllib.error import HTTPError

import pandas as pd
from Bio import Entrez
from tqdm import tqdm
import argparse

Entrez.email = "camila.ribeirodebabo@wur.nl"


def collect_gene_symbols(id_mapping_file):
    """
    Collects the gene symbols from the NCBI database using the gene IDs.
    :param id_mapping_file: File containing the mapping between the locus tag and the gene ID.
    :return: A txt file with the gene symbols and locus tags.
    """
    print("Function collect_gene_symbols() is running...")

    # Implement a dictionary to store the locus tag and corresponding gene ID for NCBI
    gene_id_dict = {}
    with open(id_mapping_file) as f:
        gene_info = [line.strip().split() for line in f]
        gene_id_dict = {line[0]: line[1] for line in gene_info} # {locus_tag: gene_id}

    # Get the gene symbol from the NCBI database using the gene ID dictionary through the Entrez API
    gene_symbol = []
    batch_size = 1000
    with tqdm(total=len(gene_id_dict.keys())) as pbar:
        for i in range(0, len(gene_id_dict.keys()), batch_size):
            batch_ids = list(gene_id_dict.keys())[i:i + batch_size]
            try:
                handle = Entrez.efetch(db='gene', id=batch_ids, retmode='xml')
                records = Entrez.read(handle)
                for record in records:
                    try:
                        gene_symbol.append(record['Entrezgene_gene']['Gene-ref']['Gene-ref_locus'])
                        # print(gene_symbol)
                    except:
                        print(f"Skipping record {record.get('Id', 'unknown')} - no 'Gene-ref_locus' key found.")

            except Exception as e:
                print(f"Error fetching records: {e}")
                continue

            except HTTPError as e:
                # If there is an HTTPError, print the error and skip to the next batch
                print(f"Error: {e}")
                continue

            # Update the progress bar for each batch
            pbar.update(batch_size)

    # Save the gene symbol and locus tag in a file
    with open('gene_symbol_mapping.txt', 'w') as f:
        for i, (locus_tag, gene_id) in enumerate(gene_info):
            try:
                f.write(f"{locus_tag} {gene_symbol[i]} \n")
            except:
                print(f"Skipping {locus_tag} - no gene symbol found.")

    print("Function collect_gene_symbols() is complete.")


def feature_df_gene_symbol(data_file, gene_symbol_file):
    """
    Replace the locus tag in the dataframe with the gene symbol.
    :param df: The dataframe containing the features.
    :param gene_symbol_file: The file containing the mapping between the locus tag and the gene symbol.
    :return: The dataframe with the gene symbol instead of the locus tag.
    """
    print("Function feature_df_gene_symbol() is running...")

    gene_symbol_dict = {}
    with open(gene_symbol_file) as f:
        gene_info = [line.strip().split() for line in f]
        gene_symbol_dict = {line[0]: line[1] for line in gene_info}

    df = pd.read_csv(data_file, sep=',', header=0, index_col=0)
    df['probe_name'] = df['locus_tag'].map(gene_symbol_dict)

    print("Function feature_df_gene_symbol() is complete.")

    return df


def edgelist_gene_symbol(edgelist_file, gene_symbol_file):
    """
    Replace the locus tag in the edgelist with the gene symbol.
    :param edgelist_file: The edgelist file.
    :param gene_symbol_file: The file containing the mapping between the locus tag and the gene symbol.
    :return: The edgelist with the gene symbol instead of the locus tag.
    """
    print("Function edgelist_gene_symbol() is running...")

    gene_symbol_dict = {}
    with open(gene_symbol_file) as f:
        gene_info = [line.strip().split() for line in f]
        gene_symbol_dict = {line[0]: line[1] for line in gene_info}

    with open(edgelist_file) as f:
        edgelist = [line.strip().split() for line in f]

    for i, j in enumerate(edgelist):
        try:
            edgelist[i][0] = gene_symbol_dict[j[0]]
            edgelist[i][1] = gene_symbol_dict[j[1]]
        except:
            print(f"Skipping {j} - no gene symbol found.")

    with open('edgelist_gene_symbol.txt', 'w') as f:
        for i in edgelist:
            f.write(f"{i[0]} {i[1]} \n")

    print("Function edgelist_gene_symbol() is complete.")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Get the gene symbol from the locus tag through the NCBI database '
                                                 'using the gene IDs.')
    parser.add_argument('-f', '--data_file', type=str, help='The dataframe containing the features.')
    parser.add_argument('-m', '--mapping_file', type=str, help='The file containing the mapping between the locus tag '
                                                               'and the gene ID.')
    parser.add_argument('-e', '--edgelist_file', type=str, help='The edgelist file.')
    parser.add_argument('function', type=str, help='The function to run.', choices=['collect_gene_symbols',
                                                                                    'feature_df_gene_symbol',
                                                                                    'edgelist_gene_symbol'])

    args = parser.parse_args()

    if args.function == 'collect_gene_symbols':
        collect_gene_symbols(args.mapping_file)
    elif args.function == 'feature_df_gene_symbol':
        feature_df_gene_symbol(args.edgelist_file, args.mapping_file)
    elif args.function == 'edgelist_gene_symbol':
        edgelist_gene_symbol(args.edgelist_file, args.mapping_file)

    # collect_gene_symbols('ncbi_id_mapping.txt')
    # get_gene_symbol_from_edgelist('../outputs/edgelist.txt', 'ncbi_id_mapping.txt')
