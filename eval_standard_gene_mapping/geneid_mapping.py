from urllib.error import HTTPError

import pandas as pd
from Bio import Entrez
from tqdm import tqdm
import argparse
import csv

Entrez.email = "camila.ribeirodebabo@wur.nl"


def collect_gene_symbols(ncbi_id_mapping_file):
    """
    Collects the gene symbols from the NCBI database using the gene IDs.
    :param ncbi_id_mapping_file: File containing the mapping between the locus tag and the gene ID.
    :return: A txt file with the gene symbols and locus tags.
    """
    print("Function collect_gene_symbols() is running...")

    # Implement a dictionary to store the locus tag and corresponding gene ID for NCBI
    gene_id_dict = {}
    with open(ncbi_id_mapping_file) as f:
        gene_info = [line.strip().split() for line in f]
        gene_id_dict = {line[0]: line[1] for line in gene_info}  # {gene_id: locus_tag}
        print(gene_id_dict)

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
                        print(gene_symbol)
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
    with open('gene_symbol_test.txt', 'w') as f:
        for i, (locus_tag, gene_id) in enumerate(gene_info):
            try:
                f.write(f"{locus_tag} {gene_symbol[i]} \n")
            except:
                print(f"Skipping {locus_tag} - no gene symbol found.")

    print("Function collect_gene_symbols() is complete.")


def locus_tag_preprocessing(locus_tag_file):
    """
    Preprocess the locus tag file to remove the probes from organs (ATC and ATM) and to separate the probes and aliases
    by ';'.
    :param locus_tag_file: The file containing the gene expression data.
    :return: A new file with the processed data.
    """
    print("Function locus_tag_preprocessing() is running...")

    # Read the locus tag file and separates the probes and aliases by ';'
    with open(locus_tag_file, 'r') as file:
        reader = csv.reader(file)
        header = next(reader)
        new_rows = []

        for row in reader:
            probe_names = [row[0][i:i + 9] for i in range(0, len(row[0]), 9)]

            # Filter out probes from organs (ATC and ATM)
            probe_names = [probe for probe in probe_names if
                           not probe.startswith('ATC') and not probe.startswith('ATM')]

            # Remove ATM and ATC from the list of probe names if they are aliases
            for i in range(len(probe_names)):
                if probe_names[i] == 'ATC' or probe_names[i] == 'ATM':
                    probe_names.pop(i)

            new_row = [';'.join(probe_names)] + row[1:]
            new_rows.append(new_row)

    # Open the CSV file for writing
    with open('probes_alias_sep.csv', 'w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(header)
        writer.writerows(new_rows)

    print("Function locus_tag_preprocessing() is complete.")


def map_locus_tag_to_gene_symbol(input_file, mapping_file):
    """
    Map the locus tag to the gene symbol.
    :param input_file: The file containing the gene expression data with the locus tag.
    :param mapping_file: The file containing the mapping between the locus tag and the gene symbol.
    :return: A new file with the gene expression data with the gene symbol.
    """

    print("Function map_locus_tag_to_gene_symbol() is running...")

    probe_to_gene = {}
    with open(mapping_file, 'r') as map_file:
        for line in map_file:
            probe, gene = line.strip().split()
            probe_to_gene[probe] = gene

    input_df = pd.read_csv(input_file)

    gene_symbols = []
    changed_genes = 0
    for probe_name in input_df['probe_name']:
        if probe_name in probe_to_gene:
            gene_symbols.append(probe_to_gene[probe_name])
            changed_genes += 1
        else:
            gene_symbols.append(probe_name)

    # Replace the "probe_name" column with "gene_symbol"
    input_df['gene_symbol'] = gene_symbols
    input_df = input_df.drop(columns=['probe_name'])
    input_df = input_df.set_index('gene_symbol')

    input_df.to_csv('symbol_test_num.csv')

    # Print number of altered genes
    print(f"Number of altered genes: {changed_genes}")
    print("Function map_locus_tag_to_gene_symbol() is complete.")


def check_flowering_genes(data_file, flowering_genes_file):
    """
    This functions checks how many flowering genes are in the dataset.
    :param data_file: A csv file containing the gene expression data.
    :param flowering_genes_file: The file containing the flowering genes
    :return: The number of flowering genes in the edgelist file
    """

    print("Function check_flowering_genes() is running...")


    flowering_genes = set()
    with open(flowering_genes_file, 'r') as file:
        for line in file:
            flowering_genes.add(line.strip())

    # Read the data_file and count the number of flowering genes
    data_file = pd.read_csv(data_file, sep=',', header=0, index_col=0)
    separated_genes = set()
    for gene in data_file.index:
        genes = gene.split(';')
        for separated_gene in genes:
            separated_genes.add(separated_gene.strip())

    matching_genes = flowering_genes.intersection(separated_genes)
    num_matching_genes = len(matching_genes)
    percentage_matching_genes = num_matching_genes / len(data_file.index) * 100

    print("Function check_flowering_genes() is complete.")
    print(f"Number of flowering genes and its percentage is {num_matching_genes} and {percentage_matching_genes}, with the present genes being {matching_genes}")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Get the gene symbol from the locus tag through the NCBI database '
                                                 'using the gene IDs.')
    parser.add_argument('-d', '--data_file', type=str, help='The file containing the gene expression data.')
    parser.add_argument('-m', '--mapping_file', type=str, help='The file containing the mapping between the '
                                                               'locus tag and the NCBI gene ID.')
    parser.add_argument('function', type=str, help='The function to run.', choices=['collect_gene_symbols',
                                                                                    'locus_tag_preprocessing',
                                                                                    'map_locus_tag_to_gene_symbol',
                                                                                    'check_flowering_genes'])

    args = parser.parse_args()

    if args.function == 'collect_gene_symbols':
        collect_gene_symbols(args.mapping_file)
    elif args.function == 'locus_tag_preprocessing':
        locus_tag_preprocessing(args.data_file)
    elif args.function == 'map_locus_tag_to_gene_symbol':
        map_locus_tag_to_gene_symbol(args.data_file, args.mapping_file)
    elif args.function == 'check_flowering_genes':
        check_flowering_genes(args.data_file, args.mapping_file)

    # collect_gene_symbols('ncbi_id_mapping.txt')
    # get_gene_symbol_from_edgelist('../outputs/edgelist.txt', 'ncbi_id_mapping.txt')
    # check_flowering_genes('../outputs/genesym.csv', 'flowering_genes.txt')
