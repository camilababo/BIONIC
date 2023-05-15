import json


def string_file_to_eval_stand(input_file):
    """
    This function takes a file containing hierarchical STRING clusters and their proteins and returns them in a standard
    format for evaluation of the clustering algorithms of the BIONIC method.
    :param input_file: A txt file containing hierarchical STRING clusters and their proteins.
    :return: A json file containing a dictionary of all clusters and their proteins.
    """
    cluster_dict = {}
    with open(input_file, 'r') as f:
        for line in f:
            if line.startswith("#"):
                continue
            tax_id, cluster_id, protein_id = line.strip().split("\t")
            cluster_id = cluster_id.replace("CL:", "CL")
            protein_id = protein_id.split(".")[1]
            if cluster_id not in cluster_dict:
                cluster_dict[cluster_id] = [protein_id]
                print(cluster_dict)
            else:
                cluster_dict[cluster_id].append(protein_id)

    return json.dumps(cluster_dict, indent=4)


if __name__ == "__main__":
    string_file_to_eval_stand("3702.clusters.proteins.v11.5.txt")

