import re
import utils

def run():
    print('Running: Cleaning H5Nx sequences.')

    for H5 in utils.subtypes:
        read_path = "data/subtype_sequences/{0}/00_{0}_raw.fasta".format(H5) # read file location
        write_path = "data/subtype_sequences/{0}/01_{0}.fasta".format(H5) # output file destination
        seq_records = utils.read_fasta(read_path) # Read fasta sequences

        # List size equal to number of sequences. Filled with False values by default. If sequence
        # meets criteria then false will be changed to True.
        keep_records = [False] * len(seq_records) 
        keep_records = check_sequence_start(seq_records, keep_records) # Check for start codon.
        keep_records = check_sequence_end(seq_records, keep_records) # Check for stop codon.
        keep_records = check_sequence_triplet(seq_records, keep_records) # Check sequence length divisible by 3.

        # Where keep_records element is True add sequence to new_seq_records.
        new_seq_records = filter_seq_records(seq_records, keep_records) 

        # Build new sequence descriptions
        # Extracts subtype, host, gene, location, and sequence length
        new_seq_records = utils.build_sequence_description(new_seq_records)

        # Remove sequences that do not contain subtype or host.
        new_seq_records = utils.remove_sequences(new_seq_records)

        utils.write_fasta(new_seq_records, write_path) # Write fasta file.
        print("Completed: {0} {1} sequences meet criteria.".format(len(new_seq_records), H5)) # Print number of sequences kept.

# Generate new list of data that is based on keep_seqs. If true then add to new list.
def filter_seq_records(seq_records, keep_records):
    kept_seq_records = []
    for i in range(0, len(seq_records)):
        if keep_records[i] == True:
            kept_seq_records.append(seq_records[i])
        else:
            pass
    return kept_seq_records

# Filters for term in sequence record description.
def keep_term(seq_records, keep_records, term):
    for i in range(0, len(seq_records)):
        record = seq_records[i]
        if keep_records[i] == True:
            if re.search(term, record.description):
                pass
            else:
                keep_records[i] = False
    return keep_records

# Filters for term in sequence record description.
def remove_term(seq_records, keep_records, term):
    for i in range(0, len(seq_records)):
        record = seq_records[i]
        if keep_records[i] == True:
            if re.search(term, record.description):
                keep_records[i] = False
            else:
                pass
    return keep_records

# If first three nucleotides are ATG change keep_record element from false to true. Else pass.
def check_sequence_start(seq_records, keep_records):
    for i in range(0, len(seq_records)):
        record = seq_records[i]
        if record.seq[0:3] == "ATG":
            keep_records[i] = True
        else:
            pass
    return keep_records

# If last three nucleotides are TAA, TGA, and TAG keep sequences. Else do not.
def check_sequence_end(seq_rcords, keep_records):
    for i in range(0, len(seq_rcords)):
        record = seq_rcords[i]
        if keep_records[i] == True:
            if record.seq[-3:] == "TAA":
                pass
            elif record.seq[-3:] == "TGA":
                pass
            elif record.seq[-3:] == "TAG":
                pass
            else:
                keep_records[i] = False
        else:
            pass
    return keep_records

# Keep sequence if length is divisible by 3.
def check_sequence_triplet(seq_records, keep_records):
    for i in range(0, len(seq_records)):
        record = seq_records[i]
        if keep_records[i] == True:
            if len(record.seq) % 3 == 0:
                pass
            else:
                keep_records[i] = False
        else:
            pass
    return keep_records
