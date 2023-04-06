from Bio import Entrez

Entrez.email = "camila.ribeirodebabo@wur.nl"


def get_gene_symbol(locus_tag, feature_dataframe, mapping_file):
    """
    Gets the gene symbol from the locus tag through the NCBI database using the gene IDs.
    :param mapping_file: The file containing the mapping between the locus tag and the gene ID.
    :param locus_tag: The locus tag of the gene.
    :param feature_dataframe: The dataframe containing the locus tags.
    :return: A dataframe with the gene symbol of the genes.
    """

    # Implement a dictionary to store the locus tag and corresponding gene ID for NCBI
    gene_id_dict = {}
    gene_info = open(mapping_file, 'r')

    for line in gene_info:
        line = line.split()
        gene_id_dict[line[0]] = line[1]

    gene_info.close()

    # Transform the locus tag column into a gene ID column
    feature_dataframe['gene_id'] = feature_dataframe[locus_tag].map(gene_id_dict)

    # Get the gene symbol from the NCBI database using the gene IDs through the Entrez API
    gene_symbol = []
    for gene_id in feature_dataframe['gene_id']:
        handle = Entrez.efetch(db='gene', id=gene_id, retmode='xml')
        records = Entrez.read(handle)
        gene_symbol.append(records[0]['Entrezgene_gene']['Gene-ref']['Gene-ref_locus'])

    # Add the gene symbol column to the dataframe
    feature_dataframe['gene_symbol'] = gene_symbol

    return feature_dataframe
