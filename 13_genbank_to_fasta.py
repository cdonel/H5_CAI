from Bio.SeqRecord import SeqRecord
import utils

def run(host):
    print('Running: Converting genbank to fasta.')
    try:
        read_path_gb = "data/genbank/{0}.gb".format(host) # read file location for genbank
        read_path_csv = "data/genbank/host_info/{0}.csv".format(host) # read file location for csv
        write_path = "data/host_sequences/{0}.fasta".format(host) # Write file destination
        gb_records = utils.read_genbank(read_path_gb) # read genbank
        seq_info = utils.csv_to_df(read_path_csv) # read csv
        seq_records = get_genbank_cds(gb_records, seq_info) # List containing sequences retrieve from genbank records
        utils.write_fasta(seq_records, write_path) # Write fasta file
    
    except:
        print("No genbank file found for {0}.".format(host))

# Takes genbank records, extracts cds, and builds sequence record.
# Looks for genbank records using CSV file built with 6_get_genbank.py
# Uses protein_id in seq_info to make sure correct cds is being extracted.
def get_genbank_cds(gb_records, seq_info):
    seq_records = [] # Empty list to hold sequence records

    for i in range(0, len(gb_records)): # loop through length of gb_records
        record = gb_records[i] # single genbank record

        if record.features: # If record attribute is features

            for feature in record.features: # for each feature in genbank record

                if feature.type == "CDS": # If feature type is equal to coding sequence (CDS)

                    try:
                        protein_id = feature.qualifiers['protein_id'][0] # parse protein id from genbank record

                    except:
                        print("{0} not added to sequence records. Does not have protein id.".format(record.name))

                    if protein_id in seq_info['protein_id'].values: # If protein id from genbank record matches protein id in <host>.csv

                        try:
                            seq = feature.location.extract(record).seq # Extract CDS sequence from genbank record
                            id = record.id # Extract id from genbank record
                            description = record.description # Extract description from genbank record
                            name = record.name # Extract name from genbank record
                            seq_record = SeqRecord(seq=seq, id=id, description=description, name=name) # Build sequence record from extracted data
                            seq_records.append(seq_record) # Add sequence record to seq_records list

                        except:
                            print("{0} was not added to sequence records".format(protein_id))
    
    return seq_records # Return list of sequence records
