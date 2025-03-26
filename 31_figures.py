import seaborn as sns
import matplotlib.pyplot as plt
import utils
import pandas as pd
import starbars


def run():
    plot_standard_cai_analysis()


# Plots standard CAI analysis. Host reference genes used are the host that the H5 strain was found in.
def plot_standard_cai_analysis(style=None):
    read_path = 'data/cai_results/scores_cai_results.csv' # read file location
    cai_results = utils.csv_to_df(read_path) # csv to dataframe

    # Plot specific subtype selecting which hosts.
    plot_subtype(cai_results, 'H5N1', ['bovine', 'chicken', 'human', 'swine'], min_y=0.55, max_y=1, bar_gap=0.01)
    plot_subtype(cai_results, 'H5N6', ['chicken', 'human', 'swine'], min_y=0.5, max_y=0.9, bar_gap=0.05)

    # Plot specifc hosts for different subtypes
    plot_genes(cai_results, host='chicken', subtypes=['H5N1', 'H5N2', 'H5N6', 'H5N8'], min_y=0.7, max_y=0.9, bar_gap=0.15)
    plot_genes(cai_results, host='duck', subtypes=['H5N1', 'H5N2', 'H5N3', 'H5N6', 'H5N8'], min_y=0.5, max_y=0.8, bar_gap=0.15)
    plot_genes(cai_results, host='human', subtypes=['H5N1', 'H5N6'], min_y=0.6, max_y=0.75, bar_gap=0.15)
    plot_genes(cai_results, host='swine', subtypes=['H5N1', 'H5N6'], min_y=0.5, max_y=0.65, bar_gap=0.15)

# Plot each gene for a given host.
def boxplot_genes(dataframe, path, stats, x, y,
                  fig_x_size, fig_y_size, title, min_y=0,max_y=1,
                  hue='subtype', bar_gap=0.03):

    sns.set_theme(rc={'figure.figsize':(fig_x_size,fig_y_size)}) # set plot figure size

    ax = sns.boxplot(dataframe, x=x, y=y, hue=hue, showfliers=False, 
                     palette="Set1", width=0.5)
    
    starbars.draw_annotation(stats, ns_show=False, fontsize=20, bar_gap=bar_gap) # statistical annotations (p values)
    plt.ylim(min_y,max_y) # y-axix limits
    plt.xlabel('Hosts', fontsize=14)
    plt.ylabel('CAI score', fontsize=14)
    plt.title(title, fontsize=16) # plot title
    plt.tick_params(axis='both', which='major', labelsize=12)
    plt.savefig(path)
    plt.clf() # clear plot data

def plot_subtype(cai_results, subtype, hosts, min_y, max_y, bar_gap):
    cai_results = cai_results[cai_results["gene"].str.contains("None|MP|NS") == False] # exlcude no gene, MP, and NS genes.
    cai_results = cai_results[cai_results["subtype"].str.contains(subtype) == True] # filter for subtype.
    host_condition = ""

    for i in range(0, len(hosts)): # Builds dataframe condition for selecting H5 subtypes

        if i != len(hosts) - 1:
            host_condition += hosts[i] + "|"

        else:
            host_condition += hosts[i]

    for gene in pd.unique(cai_results['gene']): # for each gene
        write_path = 'data/plots/{0}/{0}_{1}.jpeg'.format(subtype, gene)
        data = cai_results[cai_results["gene"].str.contains(gene) == True] # Filter for gene
        data = data[data["host"].str.contains(host_condition) == True] # Filter for hosts
        data = data.sort_values(['gene','host']) # sort gene, hosts

        match subtype: # choses which test to employ (all tests are Mann-Whitney U but generate different samples)
            case "H5N1":
                subtype_test = utils.test_H5N1(data)
            case "H5N6":
                subtype_test = utils.test_H5N6(data)

        match gene: # Set title for each plot
            case 'HA':
                title = 'Subfigure A: Hemagglutinin (HA)'
            case 'NA':
                title = 'Subfigure B: Neuraminidase (NA)'
            case 'NP':
                title = 'Subfigure C: Nucleoprotein (NP)'
            case 'PA':
                title = 'Subfigure D: RNA-dependent RNA polymerase (PA)'
            case 'PB1':
                title = 'Subfigure E: RNA-dependent RNA polymerase subunit PB1'
            case 'PB2':
                title = 'Subfigure F: RNA-dependent RNA polymerase subunit PB2'

        boxplot_genes(dataframe=data, x='host', y='cai_score', path=write_path, 
                      stats=subtype_test, fig_x_size=7.5, fig_y_size=7.5, 
                      title=title, min_y=min_y, max_y=max_y,
                      hue='host', bar_gap=bar_gap)

# Plots boxplots for each gene in a specific host.
# Genes: HA, NA, NP, PA, PB1, and PB2.
def plot_genes(cai_results, host, subtypes, min_y, max_y, bar_gap, strain_host=None):
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
                      hue='subtype', bar_gap=bar_gap)

if __name__ == "__main__":
    run()
