import test_utils as utils
import re

def main():
    for H5 in utils.subtypes:
        read_path = 'data/subtype_sequences/{0}/01_{0}.fasta'.format(H5) # Input file location
        seq_records = utils.read_fasta(read_path) # Read fasta file

        for host in utils.host_names:
            write_path = 'data/subtype_sequences/{0}/02_{0}_{1}.fasta'.format(H5, host) # Output file destination
            new_seq_records = split_sequences(seq_records, host) # Make new sequence records for each host

            if len(new_seq_records) != 0: # If length of seq_records is not zero
                utils.write_fasta(new_seq_records, write_path) # write fasta file
            

# Checks for host name name and if found puts into new sequence records list.
def split_sequences(seq_records, host):
    new_seq_records = []
    for record in seq_records:
        try:
            if re.search(rf'host={host}', record.description):
                new_seq_records.append(record)

        except:
            pass

    return new_seq_records


if __name__ == "__main__":
    main()