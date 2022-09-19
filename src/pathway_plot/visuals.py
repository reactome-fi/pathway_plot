import os
import pandas as pd
import numpy as np
from reactome2py import analysis
import plotly.express as px

from . import hierarchy

module_dir = os.path.dirname(__file__)
data_dir = "data/"

PATHWAYS = hierarchy.get_reactome_pathways()

# example test data taken from reactome analysis example input "Gene name list"
# list saved as .txt with header in example
TEST_INPUT_F = os.path.join(module_dir, data_dir, "reactome_example_gene_list.txt")

def gene_list_from_txt(fp):
    # get comma-separated gene list as string from example text file
    df = pd.read_csv(fp)
    glist = list(df[df.columns[0]])
    ids = ",".join(glist)
    return ids

def merge_enrichment_results(df):
    # merge analysis-service results with reactome pathways for plotting by top-level pathways
    # add column for additional plotting 
    df['-log10(pValue)'] = -np.log10(df['Entities pValue']) 
    return PATHWAYS.merge(df, how='left', left_on='Pathway', right_on='Pathway name')

def make_pathway_plot(df):

    fig = px.scatter(
        x=df['Pathway'],
        y=df['-log10(pValue)'],
        color=df['Top Level Pathway'],
        labels={
            "y":"-log10(pValue)",
            "x":"Pathway"},
        custom_data=[df['Pathway'],df['Stable ID'],df['Top Level Pathway'], df['Entities pValue']]
    )

    # update hover labels
    fig.update_traces(selector={'mode':'markers'},
        hovertemplate="<br>".join([
            "-log10(pValue): %{y}",
            "<b>Pathway: %{customdata[0]}</b>",
            "Stable ID: %{customdata[1]}",
            "Top Pathway: %{customdata[2]}",
            "P-value: %{customdata[3]}"]))
            #add <extra></extra> to modify box-adjacent text)
    
    # rename legend
    fig.update_layout(
        legend_title="Top Level Pathways",
        plot_bgcolor="#FFFFFF")

    fig.update_xaxes(
        gridcolor="#D8D8D8",
        showdividers=False,
        showticklabels=False,
        title='Pathways')

    fig.update_yaxes(
        rangemode='tozero',
        gridcolor="#D8D8D8",
        zerolinecolor="#000000",
        title='-log10(pValue)')

    config = {'displaylogo': False}
    fig.show(config=config)

