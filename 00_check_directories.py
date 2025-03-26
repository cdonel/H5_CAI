import os

def run():
    for dir in dirs:
        if not os.path.exists(dir):
            os.makedirs(dir)
            print('{} created'.format(dir))
    
dirs = [
    'data',
    'data/cai_results',
    'data/codon_usage_database',
    'data/genbank',
    'data/genbank/host_info',
    'data/host_sequences',
    'data/plots',
    'data/plots/chicken_strain_to_chicken_host',
    'data/plots/duck_strain_to_duck_host',
    'data/plots/H5N1',
    'data/plots/H5N6',
    'data/plots/human_strain_to_human_host',
    'data/plots/swine_strain_to_swine_host',
    'data/subtype_sequences',
    'data/subtype_sequences/H5N1',
    'data/subtype_sequences/H5N2',
    'data/subtype_sequences/H5N3',
    'data/subtype_sequences/H5N4',
    'data/subtype_sequences/H5N5',
    'data/subtype_sequences/H5N6',
    'data/subtype_sequences/H5N7',
    'data/subtype_sequences/H5N8',
    'data/subtype_sequences/meta',
    'data/subtype_sequences/phylogenetic_sequences'
    ]