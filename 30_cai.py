import utils
import codonbias as cb

def main():
    write_path = "data/cai_results/cai_results.csv"
    results = [] # Empty list for adding dataframe containing scores from each host and each subtype.
    for host in utils.host_names:
        read_path_host = "data/host_sequences/{0}.fasta".format(host) # input file location
        host_seq_records = utils.read_fasta(read_path_host) # read host fasta file
        cai_model = get_cai_model(host_seq_records) # Build CAI model using host sequences

        for H5 in utils.subtypes:
            try:
                read_path_H5 = "data/subtype_sequences/{0}/02_{0}_{1}.fasta".format(H5, host) # read file location
                seq_records = utils.read_fasta(read_path_H5) # List of sequence records
                scores = score_seqs_cai(cai_model, seq_records) # List containing scores for each sequence
                scores_df = utils.cai_scores_to_df(scores, H5, host) # Dataframe containing scores, H5 subtype, and host name
                results.append(scores_df) # Add dataframe to results list
            
            except:
                print("02_{0}_{1}.fasta file found.".format(H5, host))

    cai_results = utils.concat_df(results) # Combine dataframes with cai scores into one.
    utils.df_to_csv(cai_results, write_path) # Write dataframe to csv

# Return CAI model using host sequences
def get_cai_model(host_seq_records):
    host_seqs = [str(record.seq) for record in host_seq_records]
    cai_model = cb.scores.CodonAdaptationIndex(ref_seq=host_seqs)
    return cai_model

# Scores each sequence in list using CAI
def score_seqs_cai(cai_model, seq_records):
    seqs = [str(record.seq) for record in seq_records]
    scores = cai_model.get_score(seqs)
    return scores

if __name__ == "__main__":
    main()









