import utils
import re

def main():
    for H5 in utils.subtypes:
        read_path = 'data/subtype_sequences/{0}/01_{0}.fasta'.format(H5) # Input file location
        seq_records = utils.read_fasta(read_path) # Read fasta file
        seq_records = standardize_seq_records_hosts(seq_records) # Add standarized host name to record description

        for host in utils.host_names:
            write_path = 'data/subtype_sequences/{0}/02_{0}_{1}.fasta'.format(H5, host) # Output file destination
            new_seq_records = split_sequences(seq_records, host) # Make new sequence records for each host

            if len(new_seq_records) != 0: # If length of seq_records is not zero
                utils.write_fasta(new_seq_records, write_path) # write fasta file
            
# Stadardized host names in record description. 
def standardize_seq_records_hosts(seq_records):

    for i in range(0, len(seq_records)):
        record = seq_records[i]

        # Extract host name from sequence description using regular expression.
        # Standardize host name for duck, chicken, swine, bovine, and human in record_description.
        try:
            host_name = re.search(r'A/(.*?)/', record.description).group(1)
            new_description = utils.standarize_host(record.description, host_name)
            if new_description != None:
                seq_records[i].description = new_description
            
        except:
            pass
        
    return seq_records

# Checks for host name name and if found puts into new sequence records list.
def split_sequences(seq_records, host_name):
    new_seq_records = []
    for i in range(0, len(seq_records)):
        record = seq_records[i]
        try:

            group = re.search(rf'Influenza A virus ({host_name})', record.description).group(1)
            if group == host_name:
                new_seq_records.append(seq_records[i])

        except:
            pass

    return new_seq_records


if __name__ == "__main__":
    main()