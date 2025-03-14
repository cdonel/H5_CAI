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

# Write to text file.
def write_txt(path, data):
    with open(path, 'w') as text_file:
        text_file.write(data)

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
    return pd.read_csv(path, keep_default_na=False, na_values=[None])

# Scores to dataframe
def cai_scores_to_df(scores, subtype, host, genes):
    df = pd.DataFrame(columns=['cai_score'], data=scores)
    df['subtype'] = subtype
    df['host'] = host
    df['gene'] = genes
    return df

# Takes list of dataframes and concatenates into single datafarame
def concat_df(list_df):
    return pd.concat(list_df, ignore_index=True)

# Extracts gene name from sequence description and returns list 
def get_H5_gene_names(seq_records):
    genes = [] # empty list for holding gene name ex. HA
    for record in seq_records:
        x =  re.split(r".*?\(.*?\(.*?\)\)\s", record.description) # split everything before gene name
        x = list(filter(None, x)) # remove empty strings from list
        y = re.split(r',\s(?:partial|complete)\scds', x[0]) # split everything after gene name
        y = list(filter(None, y)) # remove empty strings from list
        gene = y[0]
        
        if re.search('hemagglutinin', gene, re.IGNORECASE): # if hemagglutinin is found in gene append HA to genes list
            genes.append('HA')

        elif re.search('(HA)', gene, re.IGNORECASE): # if (HA) is found in gene append HA to genes list
            genes.append('HA')

        elif re.search('neuraminidase', gene, re.IGNORECASE): # if neuraminidase is found in gene append NA to genes list
            genes.append('NA')

        elif re.search('(NA)', gene, re.IGNORECASE): # if (NA) is found in gene append NA to genes list
            genes.append('NA')

        elif re.search('matrix protein', gene, re.IGNORECASE): # if matrix protein is found in gene append MP to genes list
            genes.append('MP')

        elif re.search('structural', gene, re.IGNORECASE): # if structural is found in gene append NS to genes list
            genes.append('NS')

        elif re.search('PB1', gene, re.IGNORECASE): # if PB1 is found in gene append PB1 to genes list
            genes.append('PB1')

        elif re.search('PB2', gene, re.IGNORECASE): # if PB2 is found in gene append PB2 to genes list
            genes.append('PB2')

        elif re.search('(PA)', gene, re.IGNORECASE): # if PA is found in gene append PA to genes list
            genes.append('PA')

        elif re.search('(NP)', gene, re.IGNORECASE): # if NP is found in gene append NP to genes list
            genes.append('NP')
            
        else:
            genes.append("No gene") # if no gene was found above add No gene to genes list

    return genes

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

# Goes through each cleaned fasta file for each H5 subtypes and extracts gene name, subtype, and accession
# Writes dataframe to csv
def write_subtype_genes():
    subtype_genes = pd.DataFrame(columns=['subtype', 'gene', 'accession']) # empty dataframe to hold subtype, gene, and accession.
    write_path = "data/subtype_sequences/meta/subtype_genes.csv" # write file destination

    for H5 in subtypes: # for each H5 subtype
        read_path = "data/subtype_sequences/{0}/01_{0}.fasta".format(H5) # input file location
        seq_records = read_fasta(read_path) # read fasta file

        for record in seq_records:
            x =  re.split(r".*?\(.*?\(.*?\)\)\s", record.description) # split everything before gene name
            x = list(filter(None, x)) # remove empty strings from list
            y = re.split(r',\s(?:partial|complete)\scds', x[0]) # split everything after gene name
            y = list(filter(None, y)) # remove empty strings from list
            df = pd.DataFrame([{'subtype':H5, 'gene':y[0], 'accession':record.name}])
            subtype_genes = pd.concat([df, subtype_genes], ignore_index=True) # add subtype, gene, and accession

    df_to_csv(subtype_genes, write_path) # write to dataframe to csv