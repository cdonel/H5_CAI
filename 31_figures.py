import seaborn as sns
import matplotlib.pyplot as plt
import utils
import pandas as pd


def main():
    read_path = 'data/cai_results/cai_results.csv' # read file location
    cai_results = utils.csv_to_df(read_path) # csv to dataframe
    fig_1(cai_results)
    fig_2(cai_results)
    fig_3(cai_results)


# Writes dataframe to xlsx detailing sequence counts for subtype hosts and genes.
def seq_count_to_xlsx(cai_results): 
    write_path = "data/plots/sequence_table.xlsx"
    cai_results['subtype_host'] = cai_results['subtype'] + " " + cai_results['host'] # Combine subtype and hosts into one column
    cai_results = cai_results[cai_results["gene"].str.contains("No gene") == False] # Exclude genes: No gene. Not enough data to include.
    cai_results = cai_results.sort_values(['subtype', 'host', 'gene']) # sort subtype, host then gene in ascending order
    seq_count = pd.DataFrame(columns=pd.unique(cai_results['gene']), index=pd.unique(cai_results['subtype_host']))
    for gene in pd.unique(cai_results['gene']): # loop through list of unique genes
        for subtype_host in pd.unique(cai_results['subtype_host']): # loop through list of unique subtype_hosts
            df = cai_results[cai_results["gene"].str.contains(gene) == True] # filter for gene
            df = df[df["subtype_host"].str.contains(subtype_host) == True] # filter for subtype_host
            seq_count.loc[subtype_host, gene] = df.index.size # add sequence count to dataframe
    utils.df_to_xlsx(seq_count, write_path)

# H5NX subtypes: H5N1 and H5N6
# Hosts: chicken, human, swine
# Genes: HA, NA, NP, PA, PB1, and PB2
def fig_1(cai_results):
    write_path = 'data/plots/H5N1_H5N6_chicken_human_swine.png' # file destination
    cai_results = cai_results[cai_results["gene"].str.contains("No gene|MP|NS") == False] # Exclude genes: No gene, MP, and NS. Not enough data to include.
    cai_results = cai_results[cai_results["host"].str.contains("human|chicken|swine") == True] # Select chicken, human, and swine hosts
    cai_results = cai_results[cai_results["subtype"].str.contains("H5N1|H5N6") == True] # select H5N1 and H5N6 subtypes
    cai_results['legend'] = cai_results['subtype'] + " " + cai_results['host'] # Combine subtype and hosts into one column
    cai_results = cai_results.sort_values(['host', 'gene']) # sort host then gene in ascending order
    title = 'CAI for H5N1 and H5N6 genes in chickens, humans, and swine'
    x_label = 'H5N1 and H5N6 genes'
    width = 0.8
    fig_x_size = 17.5
    fig_y_size = 5
    write_boxplot(cai_results, 
                  write_path, 
                  fig_x_size, 
                  fig_y_size, 
                  width, 
                  title, 
                  x_label)

# H5NX subtypes: H5N1
# Hosts: bovine, chicken, duck, human, and swine
# Genes: HA, NA, NP, PA, PB1, and PB2
def fig_2(cai_results):
    write_path = 'data/plots/H5N1_bovine_chicken_duck_human_swine.png' # file destination
    cai_results = cai_results[cai_results["gene"].str.contains("No gene|MP|NS") == False] # Exclude genes: No gene, MP, and NS. Not enough data to include.
    cai_results = cai_results[cai_results["subtype"].str.contains("H5N1") == True] # select H5N1 and H5N6 subtypes
    cai_results['legend'] = cai_results['subtype'] + " " + cai_results['host'] # Combine subtype and hosts into one column
    cai_results = cai_results.sort_values(['host', 'gene']) # sort host then gene in ascending order
    title = 'CAI for H5N1 in bovine, chicken, duck, swine, and human'
    x_label = 'H5N1 genes'
    width = 0.8
    fig_x_size = 17.5
    fig_y_size = 5
    order = ['H5N1 chicken', 'H5N1 human', 'H5N1 bovine', 'H5N1 duck', 'H5N1 swine']
    write_boxplot(cai_results, 
                  write_path, 
                  fig_x_size, 
                  fig_y_size, 
                  width,
                  title,
                  x_label,
                  order) 
    
# H5NX subtypes: H5N1, H5N2, H5N6, H5N8
# Hosts: chicken
# Genes: HA, NA, NP, PA, PB1, and PB2
def fig_3(cai_results):
    write_path = 'data/plots/H5NX_chicken.png' # file destination
    cai_results = cai_results[cai_results["gene"].str.contains("No gene|MP|NS") == False] # Exclude genes: No gene, MP, and NS. Not enough data to include.
    cai_results = cai_results[cai_results["host"].str.contains("chicken") == True] # select chicken host
    cai_results = cai_results[cai_results["subtype"].str.contains("H5N3|H5N5") == False] # Exclude H3N3 and H5N5. Not enough data
    cai_results['legend'] = cai_results['subtype'] + " " + cai_results['host'] # Combine subtype and hosts into one column
    cai_results = cai_results.sort_values(['subtype', 'gene']) # sort subtype then gene in ascending order
    title = 'CAI for H5N1, H5N2, H5N6, and H5N8 genes in chickens'
    x_label = 'H5N1 genes'
    width = 0.8 # width of box in boxplot
    fig_x_size = 17.5
    fig_y_size = 5
    write_boxplot(cai_results, 
                  write_path, 
                  fig_x_size, 
                  fig_y_size, 
                  width,
                  title,
                  x_label) 

# Save boxplot produced to read path
def write_boxplot(dataframe, path, fig_x_size, fig_y_size, width, title, x_label, order=None):
    sns.set_theme(rc={'figure.figsize':(fig_x_size,fig_y_size)}) # set plot figure size

    # plot gene vs cai score, using subtype and host as hue, remove outliers.
    if order != None:
        ax = sns.boxplot(dataframe, x='gene', y='cai_score', hue='legend', showfliers=False, palette="Set1", width=width, hue_order=order)

    else:
        ax = sns.boxplot(dataframe, x='gene', y='cai_score', hue='legend', showfliers=False, palette="Set1", width=width)

    sns.move_legend(ax, "upper left", bbox_to_anchor=(1, 1)) # move legend outside of plot to the center right.
    plt.xlabel(x_label) # x label
    plt.ylabel('CAI score') # y label
    plt.ylim(0.5,0.9)
    plt.title(title) # plot title
    plt.savefig(path) # save plot
    plt.clf() # clear plot

if __name__ == "__main__":
    main()