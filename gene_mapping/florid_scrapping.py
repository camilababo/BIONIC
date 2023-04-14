import requests
from bs4 import BeautifulSoup

url = 'http://www.phytosystems.ulg.ac.be/florid/databases/gene_list/flowering'

# Send a GET request to the URL and get the response object
response = requests.get(url)

# Parse the HTML content using BeautifulSoup
soup = BeautifulSoup(response.content, 'html.parser')

# Find the table containing the flowering genes
table = soup.find('table', class_='list')

# Extract the gene names from the second column of each row and save them to a text file
with open('flowering_genes.txt', 'w') as f:
    for row in table.find_all('tr'):
        cells = row.find_all('td')
        if len(cells) > 1:
            gene_name = cells[1].get_text().strip()
            f.write(f"{gene_name}\n")
