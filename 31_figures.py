import seaborn as sns
import matplotlib.pyplot as plt
import utils
import pandas as pd
import starbars


def main():
    plot_standard_cai_analysis()


# Plots standard CAI analysis. Host reference genes used are the host that the H5 strain was found in.
def plot_standard_cai_analysis():
    read_path = 'data/cai_results/scores_cai_results.csv' # read file location
    cai_results = utils.csv_to_df(read_path) # csv to dataframe
    plot_genes(cai_results, host='chicken', subtypes=['H5N1', 'H5N2', 'H5N6', 'H5N8'], min_y=0.7, max_y=0.9)
    plot_genes(cai_results, host='duck', subtypes=['H5N1', 'H5N2', 'H5N3', 'H5N6', 'H5N8'], min_y=0.5, max_y=0.8)
    plot_genes(cai_results, host='human', subtypes=['H5N1', 'H5N6'], min_y=0.6, max_y=0.75)
    plot_genes(cai_results, host='swine', subtypes=['H5N1', 'H5N6'], min_y=0.5, max_y=0.65)


# Plots boxplots for each gene in a specific host.
# Genes: HA, NA, NP, PA, PB1, and PB2.
def plot_genes(cai_results, host, subtypes, min_y, max_y, strain_host=None):
    cai_results = cai_results[cai_results["gene"].str.contains("None|MP|NS") == False] # exlcude no gene, MP, and NS genes.
    cai_results = cai_results[cai_results["host"].str.contains(host) == True] # Select host 
    cai_results = cai_results.sort_values(['gene']) # sort by gene
    subtype_condition = ""

    for i in range(0, len(subtypes)): # Builds dataframe condition for selecting H5 subtypes

        if i != len(subtypes) - 1:
            subtype_condition += subtypes[i] + "|"

        else:
            subtype_condition += subtypes[i]

    for gene in pd.unique(cai_results['gene']): # for each gene

        if strain_host == None:
            write_path = 'data/plots/{0}_strain_to_{0}_host/{0}_{1}.jpeg'.format(host, gene) # file destination
        
        else: 
            write_path = 'data/plots/{0}_strain_to_{1}_host/{0}_{1}_{2}.jpeg'.format(strain_host, host, gene) # file destination

        data = cai_results[cai_results["gene"].str.contains(gene) == True] # Filter for gene
        data = data[data["subtype"].str.contains(subtype_condition) == True] # Filter for subtypes
        data = data.sort_values(['subtype']) # sort subtype

        match host: # Chooses host specific method for calculating mann-whitney u test
            case 'chicken':
                host_test = utils.test_chicken(data)
            case 'duck':
                host_test = utils.test_duck(data)
            case 'human':
                host_test = utils.test_human(data)
            case 'swine':
                host_test = utils.test_swine(data)

        boxplot_genes(dataframe=data, x='subtype', y='cai_score', path=write_path, 
                      stats=host_test, fig_x_size=7.5, fig_y_size=7.5, 
                      title=gene, min_y=min_y, max_y=max_y,
                      hue='subtype')

# Plot each gene for a given host.
def boxplot_genes(dataframe, path, stats, x, y,
                  fig_x_size, fig_y_size, title, min_y=0,max_y=1,
                  hue='subtype'):

    sns.set_theme(rc={'figure.figsize':(fig_x_size,fig_y_size)}) # set plot figure size

    sns.boxplot(dataframe, x=x, y=y, hue=hue, showfliers=False, 
                     palette="Set1", width=0.5)
    
    starbars.draw_annotation(stats, ns_show=False, fontsize=20, bar_gap=0.15) # statistical annotations (p values)
    plt.ylabel('CAI score') # y label
    plt.ylim(min_y,max_y) # y-axix limits
    plt.title(title) # plot title
    plt.savefig(path)
    plt.clf() # clear plot data

if __name__ == "__main__":
    main()