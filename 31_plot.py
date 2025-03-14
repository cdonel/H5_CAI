import seaborn as sns
import matplotlib.pyplot as plt
import utils

def main():
    read_path = 'data/cai_results/cai_results.csv' # read file location
    cai_results = utils.csv_to_df(read_path) # csv to dataframe
    plot_1(cai_results)

# H5NX subtypes: H5N1 and H5N6
# Hosts: chicken, human, swine
# Genes: HA, NA, NP, PA, PB1, and PB2
def plot_1(cai_results):
    write_path = 'data/plots/H5N1_H5N6_chicken_human_swine.png' # file destination
    cai_results = cai_results[cai_results["gene"].str.contains("No gene|MP|NS") == False] # Exclude genes: No gene, MP, and NS. Not enough data to include.
    cai_results = cai_results[cai_results["host"].str.contains("human|chicken|swine") == True] # Select chicken, human, and swine hosts
    cai_results = cai_results[cai_results["subtype"].str.contains("H5N1|H5N6") == True] # select H5N1 and H5N6 subtypes
    cai_results['legend'] = cai_results['subtype'] + " " + cai_results['host'] # Combine subtype and hosts into one column
    cai_results = cai_results.sort_values(['host', 'gene']) # sort host then gene in ascending order
    width = 0.8
    fig_x_size = 17.5
    fig_y_size = 5
    write_boxplot(cai_results, write_path, fig_x_size, fig_y_size, width) # plot

# Save boxplot produced to read path
def write_boxplot(dataframe, path, fig_x_size, fig_y_size, width):
    sns.set_theme(rc={'figure.figsize':(fig_x_size,fig_y_size)}) # set plot figure size

    # plot gene vs cai score, using subtype and host as hue, remove outliers.
    ax = sns.boxplot(dataframe, x='gene', y='cai_score', hue='legend', showfliers=False, palette="Set1", width=width)

    sns.move_legend(ax, "upper left", bbox_to_anchor=(1, 1)) # move legend outside of plot to the center right.
    plt.xlabel('H5N1 and H5N6 genes') # x label
    plt.ylabel('CAI score') # y label
    plt.title('CAI for H5N1 and H5N6 genes in chickens, humans, and swine') # plot title
    plt.savefig(path) # save plot
    plt.clf() # clear plot

if __name__ == "__main__":
    main()