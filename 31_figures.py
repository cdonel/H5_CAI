import seaborn as sns
import matplotlib.pyplot as plt
import utils
import pandas as pd
import starbars


def main():
    read_path = 'data/cai_results/cai_results.csv' # read file location
    cai_results = utils.csv_to_df(read_path) # csv to dataframe
    seq_count_to_xlsx(cai_results)
    cai_results = cai_results[cai_results["gene"].str.contains("No gene|MP|NS") == False] # Exclude genes: No gene, MP, and NS. Not enough data to include.
    #chicken(cai_results) # Mann-Whitney U rank for H5NX chicken genes
    #human(cai_results) # Mann-Whitney U rank for H5NX human genes
    #duck(cai_results) # Mann-Whitney U rank for H5NX duck genes
    #swine(cai_results) # Mann-Whitney U rank for H5NX swine genes

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

# H5NX subtypes: H5N1, H5N2, H5N6, H5N8
# Hosts: chicken
# Genes: HA, NA, NP, PA, PB1, and PB2
def chicken(cai_results):
    cai_results = cai_results[cai_results["host"].str.contains("chicken") == True] # Select chicken
    cai_results = cai_results.sort_values(['gene'])

    for gene in pd.unique(cai_results['gene']):
        write_path = 'data/plots/chicken_{0}.jpeg'.format(gene) # file destination
        data = cai_results[cai_results["gene"].str.contains(gene) == True] # Filter for gene
        data = data[data["subtype"].str.contains("H5N1|H5N2|H5N6|H5N8") == True] # Filter for subtypes
        data = data.sort_values(['subtype']) # sort subtype

        title = gene
        width = 0.5
        fig_x_size = 7.5
        fig_y_size = 7.5

        stats = utils.test_chicken(data)

        write_boxplot(
                        dataframe=data, 
                        path=write_path, 
                        stats=stats,
                        fig_x_size=fig_x_size, 
                        fig_y_size=fig_y_size, 
                        width=width, 
                        title=title, 
                        min_y=0.70,
                        max_y=0.9,
                        hue='subtype'
                        )
        
# H5NX subtypes: H5N1 and H5N6
# Hosts: human
# Genes: HA, NA, NP, PA, PB1, and PB2
def human(cai_results):
    cai_results = cai_results[cai_results["host"].str.contains("human") == True] # Select human

    for gene in pd.unique(cai_results['gene']):
        write_path = 'data/plots/human_{0}.jpeg'.format(gene) # file destination
        data = cai_results[cai_results["gene"].str.contains(gene) == True] # Filter for gene
        data = data[data["subtype"].str.contains("H5N1|H5N6") == True] # Filter for subtypes
        data = data.sort_values(['subtype']) # sort subtype

        title = gene
        width = 0.5
        fig_x_size = 7.5
        fig_y_size = 7.5

        stats = utils.test_human(data)

        write_boxplot(  
                        dataframe=data, 
                        path=write_path, 
                        stats=stats,
                        fig_x_size=fig_x_size, 
                        fig_y_size=fig_y_size, 
                        width=width, 
                        title=title, 
                        min_y=0.65,
                        max_y=0.80,
                        hue='subtype'
                        )
        
# H5NX subtypes: H5N1, H5N2, H5N3, H5N6, and H5N8
# Hosts: duck
# Genes: HA, NA, NP, PA, PB1, and PB2
def duck(cai_results):
    cai_results = cai_results[cai_results["host"].str.contains("duck") == True] # Select chicken

    for gene in pd.unique(cai_results['gene']):
        write_path = 'data/plots/duck_{0}.jpeg'.format(gene) # file destination
        data = cai_results[cai_results["gene"].str.contains(gene) == True] # Filter for gene
        data = data[data["subtype"].str.contains("H5N1|H5N2|H5N3|H5N6|H5N8") == True] # Filter for subtypes
        data = data.sort_values(['subtype']) # sort subtype

        title = gene
        width = 0.5
        fig_x_size = 7.5
        fig_y_size = 7.5

        stats = utils.test_duck(data)

        write_boxplot(  
                        dataframe=data, 
                        path=write_path, 
                        stats=stats,
                        fig_x_size=fig_x_size, 
                        fig_y_size=fig_y_size, 
                        width=width, 
                        title=title, 
                        min_y=0.5,
                        max_y=0.8,
                        hue='subtype'
                        )
        
# H5NX subtypes: H5N1 and H5N6
# Hosts: swine
# Genes: HA, NA, NP, PA, PB1, and PB2
def swine(cai_results):
    cai_results = cai_results[cai_results["host"].str.contains("swine") == True] # Select chicken

    for gene in pd.unique(cai_results['gene']):
        write_path = 'data/plots/swine_{0}.jpeg'.format(gene) # file destination
        data = cai_results[cai_results["gene"].str.contains(gene) == True] # Filter for gene
        data = data[data["subtype"].str.contains("H5N1|H5N6") == True] # Filter for subtypes
        data = data.sort_values(['subtype']) # sort subtype

        title = gene
        width = 0.5
        fig_x_size = 7.5
        fig_y_size = 7.5

        stats = utils.test_swine(data)

        write_boxplot(  
                        dataframe=data, 
                        path=write_path, 
                        stats=stats,
                        fig_x_size=fig_x_size, 
                        fig_y_size=fig_y_size, 
                        width=width, 
                        title=title, 
                        min_y=0.5,
                        max_y=0.65,
                        hue='subtype'
                        )

# Save boxplot produced to read path
def write_boxplot(
                    dataframe, 
                    path, 
                    stats,
                    fig_x_size, 
                    fig_y_size, 
                    width, 
                    title, 
                    min_y=0,
                    max_y=1,
                    hue='subtype'
                    ):

    sns.set_theme(rc={'figure.figsize':(fig_x_size,fig_y_size)}) # set plot figure size

    ax = sns.boxplot(
                        dataframe, 
                        x='subtype', 
                        y='cai_score', 
                        hue=hue, 
                        showfliers=False, 
                        palette="Set1", 
                        width=width
                        )
    
    starbars.draw_annotation(stats, ns_show=False, fontsize=20, bar_gap=0.15)
    plt.ylabel('CAI score') # y label
    plt.ylim(min_y,max_y)
    plt.title(title) # plot title
    plt.savefig(path)
    plt.clf()

if __name__ == "__main__":
    main()