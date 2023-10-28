import json

# Load the file
with open("output.json", "r") as file:
    data = json.load(file)

# Create the transformed dictionary
transformed_dict = {}
for complex_name, probes in data.items():
    for probe in probes:
        if probe not in transformed_dict:
            transformed_dict[probe] = []
        transformed_dict[probe].append(complex_name)

with open("labeled_output.json", "w") as file:
    json.dump(transformed_dict, file, indent=4)
