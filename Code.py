import requests
from Bio import SeqIO
import networkx as nx
import matplotlib.pyplot as plt

# Obtener secuencia proteica de UniProt
def get_protein_sequence(uniprot_id):
    url = f"https://www.uniprot.org/uniprot/{uniprot_id}.fasta"
    response = requests.get(url)
    sequence = SeqIO.read(response.text.splitlines(), "fasta")
    return sequence.seq

# Obtener interacciones proteína-proteína de BioGRID
def get_protein_interactions(protein_id):
    url = f"https://webservice.thebiogrid.org/interactions/?searchNames=true&geneList={protein_id}&format=json"
    response = requests.get(url)
    # Procesar la respuesta y obtener las interacciones
    interactions = process_response(response.json())
    return interactions

# Procesar la respuesta de BioGRID y obtener las interacciones relevantes
def process_response(json_data):
    interactions = []
    for interaction in json_data['results']:
        interactor_a = interaction['interactorA']
        interactor_b = interaction['interactorB']
        interactions.append((interactor_a['symbol'], interactor_b['symbol']))
    return interactions

# Obtener secuencias proteicas
protein1_sequence = get_protein_sequence("P12345")
protein2_sequence = get_protein_sequence("P67890")

# Obtener interacciones proteína-proteína
protein1_interactions = get_protein_interactions("P12345")
protein2_interactions = get_protein_interactions("P67890")

# Construir la red de interacciones proteína-proteína
G = nx.Graph()
G.add_node("Protein 1")
G.add_node("Protein 2")
for interaction in protein1_interactions:
    G.add_edge("Protein 1", interaction[1])
for interaction in protein2_interactions:
    G.add_edge("Protein 2", interaction[1])

# Dibujar la red de interacciones proteína-proteína
pos = nx.spring_layout(G)
nx.draw_networkx(G, pos=pos, with_labels=True)
plt.show()