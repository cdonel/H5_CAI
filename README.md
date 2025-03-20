# H5_CAI
## BIOT670I group alpha H5 influenza codon adaptation index (CAI)

The script has three distinct sections:

1. Retrieving and processing host reference genes into FASTA format and is denoted by scripts beginning with the number 1.
2. Retrieving and processing H5Nx gene sequences into FASTA format and is denoted by scripts beginning with the number 2.
3. Scoring CAI and plotting results  and is denoted by scripts beginning with the number 3.

### Section 1: Host reference genes
1. 10_retrieve_codon_usage_index.py parses HTML from http://www.kazusa.or.jp/codon/ which contains codon usage indexes for many species. 
URL for human codon usage index: [http://www.kazusa.or.jp/codon/cgi-bin/showcodon.cgi?species=9606](http://www.kazusa.or.jp/codon/current/species/9606)
2. 11_parse_codon_usage_index.py extracts the accession and protein ids from the codon usage index.
3. 12_retrieve_genbank.py uses the accession id to retrieve genbank records from NCBI.
4. 13_genbank_to_fasta.py uses the protein id to extract the correct cds from each genbank record.

### Section 2: H5Nx gene sequences
1. 20_retrieve_sequences.py retrieves sequences from nucleotide database in NCBI based on passed search terms.
2. 21_clean_sequences.py filters sequences on selected criteria.
3. 22_split_sequences.py splits the FASTA files into multiple FASTA files that are host species specific.

### Section 3. Analysis
1. 30_cai.py scores H5Nx sequences based on host refernce genes.
2. 31_figures.py generates boxplots and caculates statistics using Mann-Whiteny U test.

**Plot from analysis:**

CAI scores for hemagglutinin gene in chicken hosts for H5N1, H5N2, H5N6, and H5N8.
![alt text](https://github.com/cdonel/H5_CAI/blob/main/readme_images/chicken_HA.jpeg)
