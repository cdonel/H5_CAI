import utils
import pandas as pd
import re

def main():
    for H5 in utils.subtypes:
        read_path = "data/subtype_sequences/{0}/01_{0}.fasta".format(H5) # read file location
        write_path = "data/subtype_sequences/H5_hosts/{0}_hosts.csv".format(H5) # write file destination
        seq_records = utils.read_fasta(read_path) # read fasta sequences
        host_df = get_hosts(seq_records) # Extract host name and number of times it occurs into dataframe.
        # seq_records = standardize_seq_records_hosts(seq_records) # Standardize host names in seq record description
        utils.df_to_csv(host_df, write_path) # Write dataframe to csv

# Get host names and occurences from sequence description and return in dataframe.
def get_hosts(seq_records):
    temp_list = []
    for i in range(0, len(seq_records)):
        record = seq_records[i]

        # Extract host name from sequence description using regular expression.
        # Standardize host name for duck, chicken, swine, bovine, and human in record_description.
        try:
            host_name = re.search(r'A/(.*?)/', record.description).group(1)
            temp_list.append(host_name)
            
        except:
            print("No host extracted for {0}.".format(record.name))

    host_df = pd.DataFrame({'host_organism':list(set(temp_list))})


    # Count occurence of host names extract to get number of sequences for host.
    for i in range(0, host_df['host_organism'].size):
        host_df.loc[i, 'seq_count'] = int(temp_list.count(host_df.loc[i, 'host_organism']))
    
    return host_df

if __name__ == "__main__":
    main()