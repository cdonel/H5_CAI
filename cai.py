import importlib
import utils
import os

check_directories = importlib.import_module('00_check_directories')
retrieve_codon_usage_index = importlib.import_module('10_retrieve_codon_usage_index') 
parse_codon_usage_index = importlib.import_module('11_parse_codon_usage_index') 
retrieve_genbank = importlib.import_module('12_retrieve_genbank') 
genbank_to_fasta = importlib.import_module('13_genbank_to_fasta') 
retrieve_sequences = importlib.import_module('20_retrieve_sequences') 
clean_sequences = importlib.import_module('21_clean_sequences') 
split_sequences = importlib.import_module('22_split_sequences') 
score_cai = importlib.import_module('30_score_cai')
figures = importlib.import_module('31_figures')

def cai():

    check_directories.run()

    hosts()

    H5Nx()

    score_cai.run()

    figures.run()


def H5Nx():
    
    # Retrieve H5Nx sequences
    for H5 in utils.subtypes:
        if not os.path.exists('data/subtype_sequences/{0}/00_{0}_raw.fasta'.format(H5)):
            retrieve_sequences.run(H5)

    for H5 in utils.subtypes:
        if not os.path.exists('data/subtype_sequences/{0}/01_{0}.fasta'.format(H5)):
            clean_sequences.run(H5)

    clean_sequences.run()

    split_sequences.run()


def hosts():

    # Retrieve codon usage index.
    for host in utils.host_names:
        if not os.path.exists('data/codon_usage_database/{}.txt'.format(host)):
            retrieve_codon_usage_index.run(host)

    # Parse codon usage indexes
    for host in utils.host_names:
        if not os.path.exists('data/genbank/host_info/{}.csv'.format(host)):
            parse_codon_usage_index.run(host)

    # Retrieve host genbank records.
    for host in utils.host_names:
        if not os.path.exists('data/genbank/{}.gb'.format(host)):
            retrieve_genbank.run(host)
        
    # Convert host genbank to fasta.
    for host in utils.host_names:
        if not os.path.exists('data/host_sequences/{}.fasta'.format(host)):
            genbank_to_fasta.run(host)



if __name__ == '__main__':
    cai()