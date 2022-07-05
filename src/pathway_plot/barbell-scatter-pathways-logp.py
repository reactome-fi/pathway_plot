from importlib.resources import path
import plotly.express as px
import plotly.graph_objects as go
import pandas as pd
import numpy as np

import os
script_dir = os.path.dirname(__file__) 
data_dir = "../Data/"

# test input from reactome analysis
df1_f = os.path.join(script_dir, data_dir, "result_reactome_web_example.csv")
df2_f = os.path.join(script_dir, data_dir, "result_reactome_web_example_randomized.csv")

df1 = pd.read_csv(df1_f)
df2 = pd.read_csv(df2_f)

df1['-logp_sample1'] = -np.log10(df1['Entities pValue'])
df2['-logp_sample2'] = np.abs(-np.log10(df2['Entities pValue'])) # np.abs to account for negative values possible from p-val randomization

# reactome pathway metadata
pathways_f = os.path.join(script_dir, data_dir, "Pathway_List_In_Hierarchy_Release_79.txt")
pathways = pd.read_csv(pathways_f, sep='\t')

# combine results with reactome pathway metadata
dfa = pathways.merge(df1, how='left', left_on='Pathway', right_on='Pathway name')
dfb = pathways.merge(df2, how='left', left_on='Pathway', right_on='Pathway name')

def make_barbell_plot(df, df2):
    
    # # This is the default plotly categorical color seq
    # color_seq = ['#636EFA', '#EF553B', '#00CC96', '#AB63FA', '#FFA15A', '#19D3F3', '#FF6692', '#B6E880', '#FF97FF', '#FECB52']

    fig1 = px.scatter(
        x=df['Pathway'],
        y=df['-logp_sample1'],
        color=df.Top_Level_Pathway + ' (a)',
        labels={
            "y":"-log10(pValue)",
            "x":"Pathway"},
        custom_data=[df['Pathway'],df['Stable_ID'],df['Top_Level_Pathway'], df['Entities pValue']]
    )

    fig2 = px.scatter(
        x=df2['Pathway'],
        y=df2['-logp_sample2'],
        color=df2.Top_Level_Pathway + ' (b)',
        labels={
            "y":"-log10(pValue)",
            "x":"Pathway"},
        custom_data=[df2['Pathway'],df2['Stable_ID'],df2['Top_Level_Pathway'], df2['Entities pValue']],)

    # # Adding traces in categorical order for both samples fixes x-axis range issues with selection
    from plotly.subplots import make_subplots
    fig = make_subplots()
    for i in range(29):
        fig.add_trace(fig1.data[i])
        fig.add_trace(fig2.data[i])


    # Merging results to identify where barbells can be drawn
    results = df[['Pathway','-logp_sample1','Top_Level_Pathway']].merge(df2[['Pathway','-logp_sample2','Top_Level_Pathway']])
    results['score_diff'] = np.abs(results['-logp_sample1'] - results['-logp_sample2'])

    line_idx = results[~results['score_diff'].isna()].index

    # # adding bars as shapes (issues with categorical axes https://github.com/plotly/plotly.js/issues/1720 and persistence with numeric axes)
    # for i in line_idx:
    #     fig.add_shape(
    #         type="line", 
    #         #x0 = results.iloc[i].Pathway, x1 = results.iloc[i].Pathway,
    #         x0 = i, x1 = i,
    #         y0 = min(results.iloc[i]['-logp_sample1'], results.iloc[i]['-logp_sample2']),
    #         y1 = max(results.iloc[i]['-logp_sample1'], results.iloc[i]['-logp_sample2']),
    #         opacity=0.5)

    # adding bars as traces (somewhat better behavior but display is tied to a legendgroup) - likely only possible in Dash app w/ callbacks
    # can also duplicate lines so bar is display for each sample sample
    for i in line_idx:
        fig.add_traces([
            # first set of bars
            go.Scatter(
            x = [results.iloc[i].Pathway, results.iloc[i].Pathway],
            y = [
                min(results.iloc[i]['-logp_sample1'], results.iloc[i]['-logp_sample2']), 
                max(results.iloc[i]['-logp_sample1'], results.iloc[i]['-logp_sample2'])
            ],
            mode='lines',
            opacity=0.5,
            line = dict(color='black'),
            hoverinfo='skip',
            showlegend=False,
            legendgroup=results.iloc[i]['Top_Level_Pathway'] + ' (a)'),
            # duplicate bars for other legend group
            go.Scatter(
            x = [results.iloc[i].Pathway, results.iloc[i].Pathway],
            y = [
                min(results.iloc[i]['-logp_sample1'], results.iloc[i]['-logp_sample2']), 
                max(results.iloc[i]['-logp_sample1'], results.iloc[i]['-logp_sample2'])
            ],
            mode='lines',
            opacity=0.5,
            line = dict(color='black'),
            hoverinfo='skip',
            showlegend=False,
            legendgroup=results.iloc[i]['Top_Level_Pathway'] + ' (b)')])

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
        legend_title="GO Biological Process",
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


def make_scatter_plot(df):

    fig = px.scatter(
        x=df['Pathway'],
        y=df['-logp_sample1'],
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
        legend_title="GO Biological Process",
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

#make_barbell_plot(dfa, dfb)
make_scatter_plot(dfa)
