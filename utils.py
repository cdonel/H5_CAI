from Bio import SeqIO
from Bio import Entrez
import pandas as pd
import re
import config
from scipy.stats import mannwhitneyu

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

# Write dataframe to excel.
def df_to_xlsx(df, path):
    df.to_excel(path, index=True)

# Read CSV to datafrmae
def csv_to_df(path):
    return pd.read_csv(path, keep_default_na=False, na_values=[None])

# Scores to dataframe
def cai_scores_to_df(scores, seq_records):
    subtype = []
    gene = []
    host = []
    if len(scores) == len(seq_records):
        for record in seq_records:
            subtype.append(re.search(r'subtype=(.*?)\|', record.description).group(1))
            host.append(re.search(r'host=(.*?)\|', record.description).group(1))
            gene.append(re.search(r'gene=(.*?)\|', record.description).group(1))

    df = pd.DataFrame({'cai_score':scores, 'gene':gene, 'host':host, 'subtype':subtype})
    return df

# Vectors to dataframe
def cai_vectors_to_df(vectors, seq_records):
    subtype = []
    gene = []
    host = []
    for record in seq_records:
        subtype.append(re.search(r'subtype=(.*?)\|', record.description).group(1))
        host.append(re.search(r'host=(.*?)\|', record.description).group(1))
        gene.append(re.search(r'gene=(.*?)\|', record.description).group(1))

    df = pd.DataFrame({'cai_vector':vectors, 'gene':gene, 'host':host, 'subtype':subtype})
    return df

# Takes list of dataframes and concatenates into single datafarame
def concat_df(list_df):
    return pd.concat(list_df, ignore_index=True)

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

# Subtypes: H5N1, H5N2, H5N6, H5N8
# Host: chicken
# Test: Mann-Whitney U rank
def test_chicken(dataframe):
    
    H5N1_sample = dataframe[dataframe["subtype"].str.contains("H5N1") == True] # filter for H5N1 
    H5N2_sample = dataframe[dataframe["subtype"].str.contains("H5N2") == True] # filter for H5N2
    H5N6_sample = dataframe[dataframe["subtype"].str.contains("H5N6") == True] # filter for H5N6 
    H5N8_sample = dataframe[dataframe["subtype"].str.contains("H5N8") == True] # filter for H5N8 

    H5N1_H5N2_U1, H5N1_H5N2_p = mannwhitneyu(H5N1_sample['cai_score'].to_list(), # perform rank sum test between H5N1 and H5N2
                                             H5N2_sample['cai_score'].to_list(), 
                                             method="exact")
    
    H5N1_H5N6_U1, H5N1_H5N6_p = mannwhitneyu(H5N1_sample['cai_score'].to_list(), # perform rank sum test between H5N1 and H5N6
                                             H5N6_sample['cai_score'].to_list(), 
                                             method="exact")
    
    H5N1_H5N8_U1, H5N1_H5N8_p = mannwhitneyu(H5N1_sample['cai_score'].to_list(), # perform rank sum test between H5N1 and H5N8
                                             H5N8_sample['cai_score'].to_list(), 
                                             method="exact")

    H5N2_H5N6_U1, H5N2_H5N6_p = mannwhitneyu(H5N2_sample['cai_score'].to_list(), # perform rank sum test between H5N2 and H5N6
                                             H5N6_sample['cai_score'].to_list(), 
                                             method="exact")
    
    H5N2_H5N8_U1, H5N2_H5N8_p = mannwhitneyu(H5N2_sample['cai_score'].to_list(), # perform rank sum test between H5N2 and H5N8
                                             H5N8_sample['cai_score'].to_list(), 
                                             method="exact")

    H5N6_H5N8_U1, H5N6_H5N8_p = mannwhitneyu(H5N6_sample['cai_score'].to_list(), # perform rank sum test between H5N6 and H5N8
                                             H5N8_sample['cai_score'].to_list(), 
                                             method="exact")
    
    gene = pd.unique(dataframe['gene'])[0]

    stats = [                           # list containing tupes with each subtype sample and p value
        ('H5N1', 'H5N2', H5N1_H5N2_p),
        ('H5N1', 'H5N6', H5N1_H5N6_p),
        ('H5N1', 'H5N8', H5N1_H5N8_p),
        ('H5N2', 'H5N6', H5N2_H5N6_p),
        ('H5N2', 'H5N8', H5N2_H5N8_p),
        ('H5N6', 'H5N8', H5N6_H5N8_p)
        ]

    return stats

# Subtypes: H5N1 and H5N6
# Host: human
# Test: Mann-Whitney U rank
def test_human(dataframe):
    
    H5N1_sample = dataframe[dataframe["subtype"].str.contains("H5N1") == True] # filter for H5N1 
    H5N6_sample = dataframe[dataframe["subtype"].str.contains("H5N6") == True] # filter for H5N6 

    H5N1_H5N6_U1, H5N1_H5N6_p = mannwhitneyu(H5N1_sample['cai_score'].to_list(), # perform rank sum test between H5N1 and H5N6
                                             H5N6_sample['cai_score'].to_list(), 
                                             method="exact")
    
    
    gene = pd.unique(dataframe['gene'])[0]

    stats = [('H5N1', 'H5N6', H5N1_H5N6_p)] # list containing tupes with each subtype sample and p value

    return stats

# Subtypes: H5N1, H5N2, H5N3, H5N6, H5N8
# Host: duck
# Test: Mann-Whitney U rank
def test_duck(dataframe):
    
    H5N1_sample = dataframe[dataframe["subtype"].str.contains("H5N1") == True] # filter for H5N1 
    H5N2_sample = dataframe[dataframe["subtype"].str.contains("H5N2") == True] # filter for H5N2
    H5N3_sample = dataframe[dataframe["subtype"].str.contains("H5N3") == True] # filter for H5N3
    H5N6_sample = dataframe[dataframe["subtype"].str.contains("H5N6") == True] # filter for H5N6 
    H5N8_sample = dataframe[dataframe["subtype"].str.contains("H5N8") == True] # filter for H5N8

    H5N1_H5N2_U1, H5N1_H5N2_p = mannwhitneyu(H5N1_sample['cai_score'].to_list(), # perform rank sum test between H5N1 and H5N2
                                             H5N2_sample['cai_score'].to_list(), 
                                             method="exact")
    
    H5N1_H5N3_U1, H5N1_H5N3_p = mannwhitneyu(H5N1_sample['cai_score'].to_list(), # perform rank sum test between H5N1 and H5N3
                                             H5N3_sample['cai_score'].to_list(), 
                                             method="exact")
    
    H5N1_H5N6_U1, H5N1_H5N6_p = mannwhitneyu(H5N1_sample['cai_score'].to_list(), # perform rank sum test between H5N1 and H5N6
                                             H5N6_sample['cai_score'].to_list(), 
                                             method="exact")
    
    H5N1_H5N8_U1, H5N1_H5N8_p = mannwhitneyu(H5N1_sample['cai_score'].to_list(), # perform rank sum test between H5N1 and H5N2
                                             H5N8_sample['cai_score'].to_list(), 
                                             method="exact")
    
    H5N2_H5N3_U1, H5N2_H5N3_p = mannwhitneyu(H5N2_sample['cai_score'].to_list(), # perform rank sum test between H5N1 and H5N2
                                             H5N3_sample['cai_score'].to_list(), 
                                             method="exact")
    
    H5N2_H5N6_U1, H5N2_H5N6_p = mannwhitneyu(H5N2_sample['cai_score'].to_list(), # perform rank sum test between H5N1 and H5N2
                                             H5N6_sample['cai_score'].to_list(), 
                                             method="exact")
    
    H5N2_H5N8_U1, H5N2_H5N8_p = mannwhitneyu(H5N2_sample['cai_score'].to_list(), # perform rank sum test between H5N1 and H5N2
                                             H5N8_sample['cai_score'].to_list(), 
                                             method="exact")
    
    H5N3_H5N6_U1, H5N3_H5N6_p = mannwhitneyu(H5N3_sample['cai_score'].to_list(), # perform rank sum test between H5N1 and H5N2
                                             H5N6_sample['cai_score'].to_list(), 
                                             method="exact")
    
    H5N3_H5N8_U1, H5N3_H5N8_p = mannwhitneyu(H5N3_sample['cai_score'].to_list(), # perform rank sum test between H5N1 and H5N2
                                             H5N8_sample['cai_score'].to_list(), 
                                             method="exact")
    
    H5N6_H5N8_U1, H5N6_H5N8_p = mannwhitneyu(H5N6_sample['cai_score'].to_list(), # perform rank sum test between H5N1 and H5N2
                                             H5N8_sample['cai_score'].to_list(), 
                                             method="exact")
    
    gene = pd.unique(dataframe['gene'])[0]

    # list containing tupes with each subtype sample and p value
    stats = [
        ('H5N1', 'H5N2', H5N1_H5N2_p),
        ('H5N1', 'H5N3', H5N1_H5N3_p),
        ('H5N1', 'H5N6', H5N1_H5N6_p),
        ('H5N1', 'H5N8', H5N1_H5N8_p),
        ('H5N2', 'H5N3', H5N2_H5N3_p),
        ('H5N2', 'H5N6', H5N2_H5N6_p),
        ('H5N2', 'H5N8', H5N2_H5N8_p),
        ('H5N3', 'H5N6', H5N3_H5N6_p),
        ('H5N3', 'H5N8', H5N3_H5N8_p),
        ('H5N6', 'H5N8', H5N6_H5N8_p)
        ]

    return stats

# Subtypes: H5N1, H5N2, and H5N6
# Host: swine
# Test: Mann-Whitney U rank
def test_swine(dataframe):
    
    H5N1_sample = dataframe[dataframe["subtype"].str.contains("H5N1") == True] # filter for H5N1 
    H5N6_sample = dataframe[dataframe["subtype"].str.contains("H5N6") == True] # filter for H5N6 
    
    H5N1_H5N6_U1, H5N1_H5N6_p = mannwhitneyu(H5N1_sample['cai_score'].to_list(), # perform rank sum test between H5N1 and H5N6
                                             H5N6_sample['cai_score'].to_list(), 
                                             method="exact")
    
    gene = pd.unique(dataframe['gene'])[0]

    # list containing tupes with each subtype sample and p value
    stats = [
        ('H5N1', 'H5N6', H5N1_H5N6_p)
        ]

    return stats

# Writes dataframe to xlsx detailing sequence counts for subtype hosts and genes.
def seq_count_to_xlsx(cai_results): 
    write_path = "data/plots/sequence_table.xlsx"
    cai_results['subtype_host'] = cai_results['subtype'] + " " + cai_results['host'] # Combine subtype and hosts into one column
    cai_results = cai_results.sort_values(['subtype', 'host', 'gene']) # sort subtype, host then gene in ascending order
    seq_count = pd.DataFrame(columns=pd.unique(cai_results['gene']), index=pd.unique(cai_results['subtype_host']))

    for gene in pd.unique(cai_results['gene']): # loop through list of unique genes

        for subtype_host in pd.unique(cai_results['subtype_host']): # loop through list of unique subtype_hosts
            df = cai_results[cai_results["gene"].str.contains(gene) == True] # filter for gene
            df = df[df["subtype_host"].str.contains(subtype_host) == True] # filter for subtype_host
            seq_count.loc[subtype_host, gene] = df.index.size # add sequence count to dataframe
            
    df_to_xlsx(seq_count, write_path)

# searches and extracts location in sequence record
def extract_location(record):
    try:
        pattern = r'.*?\(A/.*?/(.*?)/'
        location = re.search(pattern, record.description).group(1)
    
    except:
        location = 'None'
    
    return location

# searches and extracts host in sequence record
def extract_host(record):
    try:
        pattern = r'A/(.*?)/'
        host = re.search(pattern, record.description).group(1)

    except:
        host = 'None'

    return host

# searches and extracts subtype in sequence record
def extract_subtype(record):
    try:
        pattern = r'(H\dN\d)'
        subtype = re.search(pattern, record.description).group(1)

    except:
        subtype = 'None'
        
    return subtype

# searches and extracts gene in sequence record
def extract_gene(record):
    try:
        x =  re.split(r".*?\(.*?\(.*?\)\)\s", record.description) # split everything before gene name
        x = list(filter(None, x)) # remove empty strings from list
        y = re.split(r',\s(?:partial|complete)\scds', x[0]) # split everything after gene name
        y = list(filter(None, y)) # remove empty strings from list
        gene = y[0]

        if re.search('hemagglutinin', gene, re.IGNORECASE): # if hemagglutinin is found in gene append HA to genes list
            return 'HA'

        elif re.search('(HA)', gene, re.IGNORECASE): # if (HA) is found in gene append HA to genes list
            return 'HA'
        
        elif re.search('segment 4', gene, re.IGNORECASE): # if (HA) is found in gene append HA to genes list
            return 'HA'

        elif re.search('neuraminidase', gene, re.IGNORECASE): # if neuraminidase is found in gene append NA to genes list
            return 'NA'

        elif re.search('(NA)', gene, re.IGNORECASE): # if (NA) is found in gene append NA to genes list
            return 'NA'
        
        elif re.search('sgement 6', gene, re.IGNORECASE): # if (NA) is found in gene append NA to genes list
            return 'NA'

        elif re.search('matrix protein', gene, re.IGNORECASE): # if matrix protein is found in gene append MP to genes list
            return 'MP'
        
        elif re.search('segment 7', gene, re.IGNORECASE): # if matrix protein is found in gene append MP to genes list
            return 'MP'

        elif re.search('structural', gene, re.IGNORECASE): # if structural is found in gene append NS to genes list
            return 'NS'
        
        elif re.search('segment 8', gene, re.IGNORECASE): # if structural is found in gene append NS to genes list
            return 'NS'

        elif re.search('PB1', gene, re.IGNORECASE): # if PB1 is found in gene append PB1 to genes list
            return 'PB1'
        
        elif re.search('segment 2', gene, re.IGNORECASE): # if PB2 is found in gene append PB2 to genes list
            return 'PB1'

        elif re.search('PB2', gene, re.IGNORECASE): # if PB2 is found in gene append PB2 to genes list
            return 'PB2'
        
        elif re.search('segment 1', gene, re.IGNORECASE): # if PB2 is found in gene append PB2 to genes list
            return 'PB2'

        elif re.search('(PA)', gene, re.IGNORECASE): # if PA is found in gene append PA to genes list
            return 'PA'
        
        elif re.search('segment 3', gene, re.IGNORECASE): # if PA is found in gene append PA to genes list
            return 'PA'

        elif re.search('(NP)', gene, re.IGNORECASE): # if NP is found in gene append NP to genes list
            return 'NP'
        
        elif re.search('segment 5', gene, re.IGNORECASE): # if NP is found in gene append NP to genes list
            return 'NP'
        
        else:
            return 'None'

    except:
        gene = 'None'

    return gene

# creates new sequence description containing only accession, subtype, gene, host, and location.
def build_sequence_description(seq_records):
    for record in seq_records:
        subtype = extract_subtype(record)
        host = extract_host(record)
        length = len(record.seq)
        standardized_host = standarize_host(host)
        if standardized_host == 'human':
            location = host
        else:
            location = extract_location(record)
        gene = extract_gene(record)

        record.description = '[subtype={0}|gene={1}|host={2}|location={3}|length={4}]'.format(subtype, gene, standardized_host, location, length)

    return seq_records

# Standardize host name.
# Supported host types: duck, chicken, swine, bovine, and human.
def standarize_host(host,
                    duck_keywords=duck_keywords,
                    duck_host_names=duck_host_names,
                    chicken_keywords=chicken_keywords,
                    chicken_host_names=chicken_host_names,
                    bovine_keywords=bovine_keywords,
                    bovine_host_names=bovine_host_names,
                    swine_host_names=swine_host_names,
                    human_host_names=human_host_names):
    
    # Check if host name matches ANY substring in duck_keywords
    # Return duck
    for keyword in duck_keywords:
        if keyword.lower() in host.lower():
            return 'duck'
        
    # Check if host name matches any name in duck_host_names
    # Return duck
    for name in duck_host_names:
        if name.lower() == host.lower():
            return 'duck'
        
    # Check if host name matches ANY substring in chicken_keywords
    # Return chicken
    for keyword in chicken_keywords:
        if keyword.lower() in host.lower():
            return 'chicken'
        
    # Check if host name matches any name in chicken_host_names
    # Return chicken
    for name in chicken_host_names:
        if name.lower() == host.lower():
            return 'chicken'

    # Check if host name matches ANY substring in bovine_keywords
    # Return bovine
    for keyword in bovine_keywords:
        if keyword.lower() in host.lower():
            return 'bovine'
        
    # Check if host name matches any name in bovine_host_names
    # Return bovine
    for name in bovine_host_names:
        if name.lower() == host.lower():
            return 'bovine'

    # Check if host name matches any name in swine_host_names
    # Return swine
    for name in swine_host_names:
        if name.lower() == host.lower():
            return 'swine'

    # Check if host name matches any name in human_host_names
    # Return human
    for name in human_host_names:
        if name.lower() == host.lower():
            return 'human'
        
# If subtype = None and host = None then do not add to new seq records.
def remove_sequences(seq_records):
    new_seq_records = []
    for record in seq_records:
        if re.search(r'host=None', record.description) or re.search(r'subtype=None', record.description):
            pass
        else:
            new_seq_records.append(record)

    return new_seq_records
