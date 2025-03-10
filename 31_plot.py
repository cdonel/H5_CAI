import seaborn as sns
import matplotlib.pyplot as plt
import utils

def main():
    read_path = 'data/cai_results/cai_results.csv' # read file location
    cai_results = utils.csv_to_df(read_path) # csv to dataframe
    for H5 in utils.subtypes:
        write_path = 'data/plots/{0}.png'.format(H5) # destination of write file
        results_host = cai_results.loc[cai_results['subtype'] == H5] # select data from dataframe of specific host
        write_boxplot(results_host, write_path) # Make and save boxplot

# Save boxplot produced to read path
def write_boxplot(dataframe, path):
    sns.boxplot(dataframe, x='subtype', y='cai_score', hue='host')
    plt.xlabel('Influenza A virus (H5)')
    plt.ylabel('Codon adaptation index (CAI)')
    plt.title('Codon adaptation index for {0} across hosts'.format(dataframe['subtype'].unique()[0]))
    plt.legend(loc='lower right')
    plt.savefig(path)
    plt.clf()

if __name__ == "__main__":
    main()