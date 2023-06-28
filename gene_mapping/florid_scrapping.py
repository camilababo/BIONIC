import pickle
import requests
import json


def florid_interaction_modules(florid_pickle, output_file):
    """
    Retrieves interaction complexes regarding flowering time from the JSON container of the FLORID database. It
    accesses each JSON container from a list of the locus tag genes and saves the locus tag of each gene in a
    dictionary as a key and its downstream and upstream interactions as a values."
    :param florid_pickle: Pickle file with the FLORID database
    :param output_file: Name of the output file
    :return: Dictionary with the downstream and upstram interactions of each gene
    """

    # Read the pickle file
    with open(florid_pickle, 'rb') as file:
        data = pickle.load(file)

    result = {}

    for entry in data['locustag']:
        url = f'http://www.phytosystems.ulg.ac.be/florid/details/?gene={entry}&type=json'
        response = requests.get(url)

        # Check if the request was successful
        if response.status_code == 200:
            json_data = response.json()

            if 'interactors' not in json_data['identifiers']:
                print(f'No interactors for {entry}.')
                continue

            upstream_ids = [item['id_tair'] for item in json_data['identifiers']['interactors'].get('upstream', [])]
            downstream_ids = [item['id_tair'] for item in json_data['identifiers']['interactors'].get('downstream', [])]

            if not upstream_ids and not downstream_ids:
                continue

            # Joint list of upstream and downstream interactions
            interaction_ids = upstream_ids + downstream_ids

            # Separated list of upstream and downstream interactions
            # result[entry] = {
            #     'upstream_ids': upstream_ids,
            #     'downstream_ids': downstream_ids
            # }

            result[entry] = interaction_ids

        else:
            print(f'Request failed for {entry}.')

    with open(output_file, 'w') as file:
        json.dump(result, file, indent=4)


if __name__ == '__main__':
    florid_interaction_modules('flor_id_flowering_genes.pkl', 'florid_interactions.json')