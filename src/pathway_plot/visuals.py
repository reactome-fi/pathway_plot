"""
Using reactome analysis example gene list, submit for pathway analysis using reactome2py and visualize with plotly
"""

import os
import pandas as pd
import numpy as np
from reactome2py import analysis
import plotly.express as px

script_dir = os.path.dirname(__file__) 
data_dir = "../../tests/"

# reactome pathway metadata
pathways_f = os.path.join(script_dir, data_dir, "Pathway_List_In_Hierarchy_Release_79.txt")
pathways = pd.read_csv(pathways_f, sep='\t')

# example test data taken from reactome analysis example input "Gene name list"
# list saved as .txt with header like example_data
test_input_f = os.path.join(script_dir, data_dir, "reactome_example_gene_list.txt")

def gene_list_from_txt(fp):
    # get comma-separated gene list as string from example text file
    df = pd.read_csv(fp)
    glist = list(df[df.columns[0]])
    ids = ",".join(glist)
    return ids

def merge_results(df):
    # merge analysis-service results with reactome pathways for plotting by top-level pathways
    # add column for additional plotting 
    df['-log10(pValue)'] = -np.log10(df['Entities pValue']) 
    return pathways.merge(df, how='left', left_on='Pathway', right_on='Pathway name')

def make_scatter_plot(df):

    fig = px.scatter(
        x=df['Pathway'],
        y=df['-log10(pValue)'],
        color=df.Top_Level_Pathway,
        labels={
            "y":"-log10(pValue)",
            "x":"Pathway"},
        custom_data=[df['Pathway'],df['Stable_ID'],df['Top_Level_Pathway'], df['Entities pValue']]
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

from requests.exceptions import ConnectionError
import requests
NumberTypes = (int, float, complex)

def identifiers_form(path, interactors=False, page_size='1', page='1', species='Homo Sapiens', sort_by='ENTITIES_FDR',
                     order='ASC', resource='TOTAL', p_value='1', include_disease=True, min_entities=None, max_entities=None,
                     projection=False):
    """
    Given a file path with a list of identifiers conducts reactome pathway enrichment analysis

    :param path: absolute path to the the txt file with identifier symbols to be analysed - refer to https://reactome.org/dev/analysis for format.
    :param interactors: boolean value if set to false, your query will consider only manually curated Reactome
        pathways with known biological significance. if true, your query will consider Reactome pathways that
        have been expanded by including all available protein-protein interactors from the IntAct database.
    :param page_size: page size
    :param page: number of pages
    :param species: list of species to filter the result (accepts taxonomy ids, species names and dbId)
    :param sort_by: how to sort the result. Available filters: TOTAL_ENTITIES, TOTAL_REACTIONS, TOTAL_INTERACTIONS, FOUND_ENTITIES,
        FOUND_INTERACTIONS, FOUND_REACTIONS, ENTITIES_RATIO, ENTITIES_PVALUE, ENTITIES_FDR, REACTIONS_RATIO
    :param order: order ASC or DESC
    :param resource: the resource to sort TOTAL, UNIPORT, ENSEMBLE, CHEMBI, IUPHAR, MIRBASE, NCBI_PROTEIN, EMBL, COMPOUND, PUBCEM_COMPOUND
    :param p_value: defines the pValue threshold. Only hit pathway with pValue equals or below the threshold will be returned
    :param include_disease: set to ‘false’ to exclude the disease pathways from the result (it does not alter the statistics)
    :param projection: if true, projects the identifiers to human and only shows the result in this species
    :param max_entities: maximum number of contained entities per pathway (takes into account the resource)
    :param min_entities: minimum number of contained entities per pathway (takes into account the resource)
    :return:
    """

    if isinstance(page_size, NumberTypes):
        page_size = str(page_size)

    if isinstance(page, NumberTypes):
        page = str(page)

    if isinstance(p_value, NumberTypes):
        p_value = str(p_value)

    if isinstance(min_entities, NumberTypes):
        min_entities = str(min_entities)

    if isinstance(max_entities, NumberTypes):
        max_entities = str(max_entities)

    if interactors:
        interactors = 'true'
    else:
        interactors = 'false'

    if include_disease:
        include_disease = 'true'
    else:
        include_disease = 'false'

    headers = {
        'content-type': 'text/plain',
    }

    params = (
        ('interactors', interactors),
        ('pageSize', page_size),
        ('page', page),
        ('sortBy', sort_by),
        ('order', order),
        ('species',  species),
        ('resource', resource),
        ('pValue', p_value),
        ('includeDisease', include_disease),
        ('min', min_entities),
        ('max', max_entities),
    )

    if projection:
        url = 'https://reactome.org/AnalysisService/identifiers/form/projection'
    else:
        url = 'https://reactome.org/AnalysisService/identifiers/'

    #data = open(path, 'rb').read()
    files = {'file': open(path, 'rb')}
    with open(path, 'rb') as file:
        try:
            response = requests.post(url=url, headers=headers, params=params, data=file)
        except ConnectionError as e:
            print(e)

        if response.status_code == 200:
            return response.json()
        else:
            print('Status code returned a value of %s' % response.status_code)

# parse input, submit to analysis service, get results as df
#ids = gene_list_from_txt(test_input_f)

# ids = gene_list_from_txt(open(test_input_f))
response = analysis.identifiers(open(test_input_f,'rb'))
#response = identifiers_form(test_input_f)
token = response['summary']['token']
results = analysis.pathway2df(token)

# convert dtypes - results df returned with all string dtypes
convert = {
    "#Entities found": "int64",
    "#Entities total": "int64",
    "Entities ratio": "float64",
    "Entities pValue": "float64",
    "Entities FDR": "float64",
    "#Reactions found": "int64",
    "#Reactions total": "int64",
    "Reactions ratio": "float64"
}
results = results.astype(convert)

# merge results with complete reactome pathways 
#df = merge_results(results)
#make_scatter_plot(df)
#results.to_csv('example_ouput.csv')