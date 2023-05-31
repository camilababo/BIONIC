import pandas as pd
import numpy as np
import json
from sklearn.metrics import confusion_matrix


def compare_predictions_with_true_labels(tsv_file, json_file):
    tsv_data = pd.read_csv(tsv_file, sep='\t', index_col=0)

    gene_names = tsv_data.index
    cluster_scores = tsv_data.values

    with open(json_file, 'r') as f:
        gene_labels_data = json.load(f)

    for gene_name, cluster_score in zip(gene_names, cluster_scores):
        true_cluster_label = gene_labels_data.get(gene_name, 'Label not found')

        print(f"Gene: {gene_name}")

        for cluster, score in enumerate(cluster_score, start=1):
            print(f"Cluster {cluster} Score: {score}")

        print(f"True Cluster Label: {true_cluster_label}")
        print("------")


def compare_predictions_confusion_matrix(tsv_file, json_file):
    # Read the TSV file
    tsv_data = pd.read_csv(tsv_file, sep='\t', index_col=0)

    # Extract gene names and cluster scores
    gene_names = tsv_data.index
    cluster_scores = tsv_data.values

    # Read the JSON file
    with open(json_file, 'r') as f:
        gene_labels_data = json.load(f)

    # Create lists to store true labels and predicted labels
    true_labels = []
    predicted_labels = []

    # Iterate over gene names, cluster scores, and true labels
    for gene_name, cluster_score in zip(gene_names, cluster_scores):
        true_label = gene_labels_data.get(gene_name, 'Label not found')
        true_labels.append(true_label)

        predicted_label = np.argmax(cluster_score) + 1  # Get the index of the maximum score as the predicted label
        predicted_labels.append(predicted_label)

        # print(f"Gene: {gene_name}")
        # print(f"True Cluster Label: {true_label}")
        # print(f"Predicted Cluster Label: {predicted_label}")
        # print("------")

    # Check if gene names are matching between TSV file and JSON file
    print(f"Number of matching gene names: {len(true_labels)}")

    # Convert true labels and predicted labels to numpy arrays
    true_labels = np.array(true_labels)
    predicted_labels = np.array(predicted_labels)

    # Compute confusion matrix
    confusion_mat = confusion_matrix(true_labels, predicted_labels)

    # Print confusion matrix
    print("Confusion Matrix:")
    print(confusion_mat)


if __name__ == '__main__':
    # compare_predictions_with_true_labels('../labelled_run/2labeled_features_features.tsv', 'labeled_output.json')
    compare_predictions_confusion_matrix('../labelled_run/2labeled_features_features.tsv', 'labeled_output.json')
