import utils
import re
import pandas as pd

def main():
    for host in utils.host_names:
        try:
            read_path = "data/codon_usage_database/{0}.txt".format(host) # Read text file
            write_path = "data/genbank/host_info/{0}.csv".format(host) # Write file destination
            seq_info = get_codon_db_info(read_path) # Dataframe with accession ids and protein ids
            utils.df_to_csv(seq_info, write_path) # Write dataframe to CSV
        except:
            print("No {0}.txt file found.".format(host))

# Parse accession id and protein id from codon usage database text file.
# Opens text file and reads by line. Uses regex to find IDs.
def get_codon_db_info(read_path):
    accessions = [] # Empty list for accession ids
    protein_ids = [] # Empty list for protein ids
    with open(read_path) as f:
        for line in f:
            if line[0] == '>': # If line begins with >
                accession_pattern = r'>(\w*)' # Accession id pattern
                protein_id_pattern = r'protein_id="(.*?)"' # Protein id pattern
                try:
                    accession = re.search(accession_pattern, line).group(1) # accession id
                    protein_id = re.search(protein_id_pattern, line).group(1) # protein id
                    accessions.append(accession) # add accession id to list
                    protein_ids.append(protein_id) # add protein id to list
                except:
                    print('Did not extract for {0}'.format(accession))

    seq_info = pd.DataFrame({'accession':accessions, 'protein_id':protein_ids}) # Create dataframe with accession and protein ids
    return seq_info # Return dataframe
                
if __name__ == '__main__':
    main()
