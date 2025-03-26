from Bio import Entrez
from config import email
from utils import subtypes

def run():
    # Provide user email for NCBI account.
    Entrez.email = email

    # Uses each H5 subtype as the search term in NCBI's nucleotide database
    # and writes sequences to fasta file.
    for H5 in subtypes:
        acc_list, webenv, query_key = get_stream_data(H5)
        get_sequences(acc_list, webenv, query_key, H5)

# Return NCBI stream data. Ideal when dealing with many sequences.
def get_stream_data(H5):
    print("Fetching steam data for {0}.".format(H5))
    stream = Entrez.esearch(db="nucleotide", term=H5, usehistory="y", idtype="acc", retmax=1000000)
    stream_data = Entrez.read(stream)
    stream.close()
    return stream_data["IdList"], stream_data["WebEnv"], stream_data["QueryKey"]

# Use stream data to access get sequences and write to file.
def get_sequences(acc_list, webenv, query_key, H5):
    batch_size = 1000
    count = len(acc_list)
    path = "data/subtype_sequences/{0}/00_{0}_raw.fasta".format(H5)
    output = open(path, "w")
    print("Fetching sequences for {0}.".format(H5))
    for start in range(0, count, batch_size):
        end = min(count, start + batch_size)
        print("Downloading record %i to %i" % (start + 1, end))

        stream = Entrez.efetch(
            db="nucleotide",
            rettype="fasta",
            retmode="text",
            retstart=start,
            retmax=batch_size,
            webenv=webenv,
            query_key=query_key,
            idtype="acc",
        )   
        
        data = stream.read()
        stream.close()
        output.write(data)

    output.close()
    print("Created 00_{0}_raw.fasta.".format(H5))

if __name__ == "__main__":
    run()
