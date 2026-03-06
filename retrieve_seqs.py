from Bio import Entrez

Entrez.email = "kiked38026@docsfy.com" 

def retrieve_ref_sequence():
    handle = Entrez.efetch(db="nucleotide", id="NC_045512.2", rettype="fasta", retmode="text")
    data = handle.read()
    handle.close()
    
    with open("data/ref.fasta", "w") as f:
        f.write(data)

def retrieve_sequences(qtd: int, term: str = '"Severe acute respiratory syndrome coronavirus 2"[Organism] AND "complete genome"[Title]',):
    search_handle = Entrez.esearch(db="nucleotide", term=term, retmax=qtd, idtype="acc")
    record = Entrez.read(search_handle)
    search_handle.close()

    id_list = record["IdList"]
    
    if id_list:
        print(f"Downloading sequences...")
        fetch_handle = Entrez.efetch(db="nucleotide", id=id_list, rettype="fasta", retmode="text")
        data = fetch_handle.read()
        fetch_handle.close()
        
        with open("data/sequences.fasta", "w") as f:
            f.write(data)
    else:
        print("None sequence found.")   
    
if __name__ == "__main__":
    retrieve_ref_sequence()
    retrieve_sequences(100)