# directory structure for the project
import os
import subprocess
if not os.path.exists('./data'): 
    os.makedirs('./data')
if not os.path.exists('./data/fastas'): 
    os.makedirs('./data/fastas')
if not os.path.exists('./data/pdbs'): 
    os.makedirs('./data/pdbs')
if not os.path.exists('./data/pqrs'):
    os.makedirs('./data/pqrs')

# required files for the project
urls = {"https://drive.google.com/uc?export=download&id=1N1Z-oy5obzkASoe2L94sWxZs5x-SEp_Y": os.path.abspath("./data/gbr_model1.joblib"),
        "https://drive.google.com/uc?export=download&id=1sy-kBCufGMH0kEcSyFEcPHkkWoscMy1u": os.path.abspath("./data/pca1.joblib"),
        "https://drive.google.com/uc?export=download&id=16nsXGQsxao55S9xFoCvYjF-jxBbM5KRf": os.path.abspath("./data/prokaryotes_348columns.csv"),
        "https://drive.google.com/uc?export=download&id=1MmKTMPqI8ZiWHanYxFKip58jzqyAMLlJ": os.path.abspath("./data/prokaryotes1.joblib"),
        "https://drive.google.com/uc?export=download&id=15-5v0Ajhokl53_2om08yGevDegcYgGN9": os.path.abspath("./data/scaler1.joblib"),
        "https://drive.google.com/uc?export=download&id=1S7TvKjm-dqm2E6LyZGZJw2lIi3uJutU5": os.path.abspath("./data/essential_proteins_prediction_merged.csv"),
        "https://drive.google.com/uc?export=download&id=1qLkD_ui7bSS2kiHQPS0dwVeT8tnp9XiA": os.path.abspath("./data/essential_proteins.csv"),
        "https://drive.google.com/uc?export=download&id=1McFlKIsZQD7DkvVYNSqHIU401BLlM7q9": os.path.abspath("./data/prokaryotes_322columns.csv"),
        "https://drive.google.com/uc?export=download&id=1NaX3aK3goCVUfYTR4hR75VIULGSMOA8-": os.path.abspath("./data/prokaryotes_all.csv"),}
for url, output_path in urls.items():
    if not os.path.exists(output_path):
        print(f"Downloading {output_path}")
        command = [
        "curl", #download from url
        "-L",   #follow redirects
        "-o", output_path,
        url
        ]
        subprocess.run(command, check=True)
    else:
        print(f"{output_path} already exists")
        
# download of 2 example PDB files
pdburls = {'https://drive.google.com/uc?export=download&id=1Vt4tsez5u2PheQXRh8NFUYX3evqkj8KP':os.path.abspath('./data/pdbs/AF-P94501-F1.pdb'),
            'https://drive.google.com/uc?export=download&id=14Y7dra6qzjAQS7KGN0rdncpcNYspRQdh': os.path.abspath('./data/pdbs/AF-P09030-F2.pdb'),}
for url, output_path in pdburls.items():
    if not os.path.exists(output_path):
        print(f"Downloading {output_path}")
        command = [
        "curl", #download from url
        "-L",   #follow redirects
        "-o", output_path,
        url
        ]
        subprocess.run(command, check=True)
    else:
        print(f"{output_path} already exists")
        
from .proteinclass import Protein
from .SPARC import SPARC
from .ThERMOS import ThERMOS, ThERMless
__all__ = ['Protein', 'SPARC', 'ThERMOS', 'ThERMless']