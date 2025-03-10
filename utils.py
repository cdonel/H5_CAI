from Bio import SeqIO
from Bio import Entrez
import pandas as pd
import re
import config

# Dictionary of host names and URLs that point to codon usage index text files.
# Codon usage database: http://www.kazusa.or.jp/codon/
codon_usage_index_urls = {
                        "chicken" : "http://www.kazusa.or.jp/codon/current/species/208526",
                        "bovine" : "http://www.kazusa.or.jp/codon/current/species/9913",
                        "swine" : "http://www.kazusa.or.jp/codon/current/species/9823",
                        "duck" : "http://www.kazusa.or.jp/codon/current/species/8839",
                        "human" : "http://www.kazusa.or.jp/codon/current/species/9606"
}

# List of all subtypes for H5 influenza A.
subtypes = [
                        "H5N1",
                        "H5N2",
                        "H5N3",
                        "H5N4",
                        "H5N5",
                        "H5N6",
                        "H5N7",
                        "H5N8"
        ]

# List of hosts that will be used in CAI analysis
host_names = [
                        'duck',
                        'chicken',
                        'swine',
                        'bovine',
                        'human'
    ]

# Kewords for duck host names
duck_keywords = [
                        'duck', 
                        'mallard'
            ]

# List of host names for duck.
duck_host_names = [
                        'waterfowl', 
                        'peafowl', 
                        'water fowl', 
                        'wildfowl', 
                        'anas platyrhynchos',
                        'merganser',
                        'mergus'
            ]

# Keywords for chicken host names
chicken_keywords = [
                        'chicken',
                        'gallus'
            ]

# List of host names for chicken.
chicken_host_names = [
                        'poultry'
                        'guinea fowl', 
                        'guineafowl',
            ]

# List of keywords for bovine.
bovine_keywords = [
                        'bovine', 
                        'cow'
            ]

bovine_host_names = [
                        'dairy cattle',
                        'cattle',
                        'cattle milk product'
]

# List of host names for swine.
swine_host_names = [
                        'pig', 
                        'swine'
            ]

# List of host names for human.
human_host_names = [
                        'human',
                        'antofagasta',
                        'anhui',
                        'bangladesh',
                        'cambodia',
                        'china',
                        'egypt',
                        'england',
                        'guangdong',
                        'guangzhou',
                        'hongkong',
                        'hong kong',
                        'indonesia',
                        'iraq',
                        'jiangsu',
                        'jianxi',
                        'kienGiang',
                        'nanjing',
                        'puerto rico',
                        'prachinburi',
                        'shenzhen',
                        'sichuan',
                        'thailand',
                        'vietnam',
                        'viet nam',
                        'washington'
                        'xinjiang',
                        'yangzhou',
                        'yunnan',
                        'zhejiang'
]

# Read fasta sequences into list. Each element is a biopython seqrecord.
def read_fasta(path):
    seq_records = []
    for record in SeqIO.parse(path, "fasta"):
        seq_records.append(record)
    return seq_records

# Write fasta file by passing list containing seqrecords.
def write_fasta(seq_records, path):
    with open(path, "w") as output_handle:
        for i in seq_records:
            SeqIO.write(i, output_handle, "fasta")

# Retrieves genbank files and returns handle.
def get_genbank_handle(accession_ids):
    Entrez.email = config.email
    handle = Entrez.efetch(db='nucleotide', id=accession_ids, rettype='gb', retmode='text')
    return handle

# Writes genbank file using handle from get_genbank_handle
def write_genbank(handle, path):
    with open(path, 'w') as f:
        f.write(handle.read())

# Read genbank files to list.
def read_genbank(path):
    gb_records = []
    for record in SeqIO.parse(path, "genbank"):
        gb_records.append(record)
    return gb_records

# Write dataframe to CSV.
def df_to_csv(df, path):
    df.to_csv(path, index=False)

# Read CSV to datafrmae
def csv_to_df(path):
    return pd.read_csv(path)

# Scores to dataframe
def cai_scores_to_df(scores, subtype, host):
    df = pd.DataFrame(columns=['cai_score'], data=scores)
    df['subtype'] = subtype
    df['host'] = host
    return df

# Takes list of dataframes and concatenates into single datafarame
def concat_df(list_df):
    return pd.concat(list_df, ignore_index=True)

# Standardize host name.
# Supported host types: duck, chicken, swine, bovine, and human.
def standarize_host(record_description, 
                    host_name,
                    duck_keywords=duck_keywords,
                    duck_host_names=duck_host_names,
                    chicken_keywords=chicken_keywords,
                    chicken_host_names=chicken_host_names,
                    bovine_keywords=bovine_keywords,
                    bovine_host_names=bovine_host_names,
                    swine_host_names=swine_host_names,
                    human_host_names=human_host_names):
    
    # Check if host name matches ANY substring in duck_keywords
    # Change section of record descrition from 'Influenza A virus' to 'Influenza A virus duck'
    for keyword in duck_keywords:
        if keyword.lower() in host_name.lower():
            record_description = re.sub(r'Influenza A virus', 'Influenza A virus duck', record_description)
            return record_description
        
    # Check if host name matches any name in duck_host_names
    # Change section of record descrition from 'Influenza A virus' to 'Influenza A virus duck'
    for name in duck_host_names:
        if name.lower() == host_name.lower():
            record_description = re.sub(r'Influenza A virus', 'Influenza A virus duck', record_description)
            return record_description
        
    # Check if host name matches ANY substring in chicken_keywords
    # Change section of record descrition from 'Influenza A virus' to 'Influenza A virus chicken'
    for keyword in chicken_keywords:
        if keyword.lower() in host_name.lower():
            record_description = re.sub(r'Influenza A virus', 'Influenza A virus chicken', record_description)
            return record_description
        
    # Check if host name matches any name in chicken_host_names
    # Change section of record descrition from 'Influenza A virus' to 'Influenza A virus chicken'
    for name in chicken_host_names:
        if name.lower() == host_name.lower():
            record_description = re.sub(r'Influenza A virus', 'Influenza A virus chicken', record_description)
            return record_description

    # Check if host name matches ANY substring in bovine_keywords
    # Change section of record descrition from 'Influenza A virus' to 'Influenza A virus bovine'
    for keyword in bovine_keywords:
        if keyword.lower() in host_name.lower():
            record_description = re.sub(r'Influenza A virus', 'Influenza A virus bovine', record_description)
            return record_description
        
    # Check if host name matches any name in bovine_host_names
    # Change section of record descrition from 'Influenza A virus' to 'Influenza A virus bovine'
    for name in bovine_host_names:
        if name.lower() == host_name.lower():
            record_description = re.sub(r'Influenza A virus', 'Influenza A virus bovine', record_description)
            return record_description

    # Check if host name matches any name in swine_host_names
    # Change section of record descrition from 'Influenza A virus' to 'Influenza A virus swine'
    for name in swine_host_names:
        if name.lower() == host_name.lower():
            record_description = re.sub(r'Influenza A virus', 'Influenza A virus swine', record_description)
            return record_description

    # Check if host name matches any name in human_host_names
    # Change section of record descrition from 'Influenza A virus' to 'Influenza A virus human'
    for name in human_host_names:
        if name.lower() == host_name.lower():
            record_description = re.sub(r'Influenza A virus', 'Influenza A virus human', record_description)
            return record_description

    