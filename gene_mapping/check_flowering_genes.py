def check_flowering_genes(edgelist_file, flowering_genes_file):
    """
    This functions checks how many flowering genes are in the edgelist file after filtering.
    :param edgelist_file: The edgelist file
    :param flowering_genes_file: The file containing the flowering genes
    :return: The number of flowering genes in the edgelist file
    """

    # Read the flowering genes file
    with open(flowering_genes_file, 'r') as f:
        flowering_genes = {}
        for line in f:
            genes = line.strip().split(',')
            flowering_genes[genes[0]] = genes[1:]

    # Check if the genes or their aliases are present in the edgelist file
    flowering_genes_in_edgelist = set()
    with open(edgelist_file, 'r') as f:
        for line in f:
            genes = line.strip().split()
            gene1, gene2 = genes[0], genes[1]
            if gene1 in flowering_genes:
                flowering_genes_in_edgelist.add(gene1)
            elif gene2 in flowering_genes:
                flowering_genes_in_edgelist.add(gene2)
            else:
                for alias in flowering_genes.get(gene1, []):
                    if alias in genes:
                        flowering_genes_in_edgelist.add(gene1)
                        break
                for alias in flowering_genes.get(gene2, []):
                    if alias in genes:
                        flowering_genes_in_edgelist.add(gene2)
                        break

    # Return the number of flowering genes in the edgelist file
    return len(flowering_genes_in_edgelist)

