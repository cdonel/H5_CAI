import utils
import codonbias as cb
import numpy as np

def run():
    standard_cai_analysis()

# Measures CAI of gene found in host to the host's reference genes.
# Example: H5N1 found in duck is score against duck reference genes.
def standard_cai_analysis():
    score_write_path = "data/cai_results/scores_cai_results.csv"
    vector_write_path = "data/cai_results/vector_cai_results.csv"
    mean_results = [] # Empty list for adding dataframe containing scores from each host and each subtype.
    vector_results = []
    for host in utils.host_names:
        read_path_host = "data/host_sequences/{0}.fasta".format(host) # input file location
        host_seq_records = utils.read_fasta(read_path_host) # read host fasta file
        cai_model = get_cai_model(host_seq_records) # Build CAI model using host sequences

        for H5 in utils.subtypes:
            try:
                read_path_H5 = "data/subtype_sequences/{0}/02_{0}_{1}.fasta".format(H5, host) # read file location
                seq_records = utils.read_fasta(read_path_H5) # List of sequence records
                scores, vectors = score_seqs_cai(cai_model, seq_records) # List containing scores for each sequence
                scores_df = utils.cai_scores_to_df(scores, seq_records) # Dataframe containing scores, H5 subtype, and host name
                vectors_df = utils.cai_vectors_to_df(vectors, seq_records)
                mean_results.append(scores_df) # Add dataframe to results list
                vector_results.append(vectors_df)
            
            except Exception as e:
                print(e)

    mean_cai_results = utils.concat_df(mean_results) # Combine dataframes with cai scores into one.
    vector_cai_results = utils.concat_df(vector_results)
    utils.df_to_csv(mean_cai_results, score_write_path) # Write dataframe to csv
    utils.df_to_csv(vector_cai_results, vector_write_path) # Write dataframe to csv

# Return CAI model using host sequences
def get_cai_model(host_seq_records):
    host_seqs = [str(record.seq) for record in host_seq_records]
    cai_model = cb.scores.CodonAdaptationIndex(ref_seq=host_seqs)
    return cai_model

# Scores each sequence in list using CAI
def score_seqs_cai(cai_model, seq_records):
    seqs = [str(record.seq) for record in seq_records]
    scores = cai_model.get_score(seqs)
    vectors = cai_model.get_vector(seqs)
    vectors = [vector[~np.isnan(vector)] for vector in vectors]
    return scores, vectors

if __name__ == "__main__":
    run()









