import utils

def run():
    for host in utils.host_names:
        try:
            read_path = "data/genbank/host_info/{0}.csv".format(host) # Input file location
            write_path = "data/genbank/{0}.gb".format(host) # Write file destination
            seq_info = utils.csv_to_df(read_path) # Dataframe with accession and protein ids
            accession_ids = seq_info['accession'].to_list() # List of accession ids
            handle = utils.get_genbank_handle(accession_ids=accession_ids) # Retrieves genbank handle using list of accession ids
            utils.write_genbank(handle, write_path) # Write genbank file

        except:
            print("No host csv file for {0}.".format(host))

if __name__ == '__main__':
    run()