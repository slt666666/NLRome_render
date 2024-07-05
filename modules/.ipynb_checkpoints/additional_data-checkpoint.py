import pandas as pd
import numpy as np
import plotly.graph_objs as go
from plotly.subplots import make_subplots
import pickle
import pathlib

                                            
def make_expression_figure(exp_data, check_chosen, target_NLRs):
    fig = go.Figure()
    for exp_name in check_chosen:
        fig.add_trace(
            go.Bar(
                y=exp_data.loc[target_NLRs, :].loc[:, exp_name].values, 
                x=target_NLRs,
                name=exp_name
            ),
        )
    fig.update_xaxes(
        autotickangles=[0, 75, 90],
        tickfont=dict(
            size=10
        )
    )
    fig.update_layout(
        margin=dict(t=10, b=10, l=10, r=10),
        legend=dict(
            yanchor="top",
            y=0.99,
            xanchor="right",
            x=0.99
        )
    )
    return fig

def make_chippeak_figure(histone, treat, table_data, DATA_PATH):
    # read data
    with open(DATA_PATH.joinpath(f'chippeak_for_dash/{treat}_{histone}_ChIP_seq_Rep1.pic.bin'), 'rb') as f:
        rep1 = pickle.load(f)
    with open(DATA_PATH.joinpath(f'chippeak_for_dash/{treat}_{histone}_ChIP_seq_Rep2.pic.bin'), 'rb') as f:
        rep2 = pickle.load(f)
    
    figs = []
    fig = make_subplots(rows=8, cols=1,
                        shared_xaxes=True,
                        vertical_spacing=0.03,
                        row_heights=[0.2, 0.1, 0.05, 0.2, 0.1, 0.05, 0.2, 0.1])
    for k, i in enumerate(range(table_data.shape[0])):
        tmp_NLR = table_data.iloc[i, :]["NLR id"]
        tmp_pos = (table_data.iloc[i, :].Start + table_data.iloc[i, :].End)/2
        tmp_strand = table_data.iloc[i, :].Strand
        tmp_exon = eval(table_data.iloc[i, :].Exon)
        
        rep1_region = rep1[f"{tmp_NLR}_region"]
        rep1_region_width = rep1[f"{tmp_NLR}_region_width"]
        rep1_depth = rep1[f"{tmp_NLR}_depth"]
        rep2_region = rep2[f"{tmp_NLR}_region"]
        rep2_region_width = rep2[f"{tmp_NLR}_region_width"]
        rep2_depth = rep2[f"{tmp_NLR}_depth"]
        
        fig.add_trace(
            go.Bar(
                x=np.array(rep1_region) - rep1_region[0], 
                y=rep1_depth, 
                width=rep1_region_width, 
                marker=dict(line=dict(width=0)),
                hoverinfo="text",
                hovertext="Rep1",
                marker_color='orange',
                opacity=0.75,
            ),
            row=1+k*3, col=1
        )
        fig.add_trace(
            go.Bar(
                x=np.array(rep2_region) - rep1_region[0],
                y=rep2_depth, 
                width=rep2_region_width, 
                marker=dict(line=dict(width=0)),
                hoverinfo="text",
                hovertext="Rep2",
                marker_color='slateblue',
                opacity=0.75,
            ),
            row=1+k*3, col=1
        )
        for ind, j in enumerate(tmp_exon):
            if ind > 0:
                fig.add_trace(
                    go.Scatter(
                        x=np.array([tmp_exon[ind-1][1]/100, j[0]/100]) - rep1_region[0],
                        y=[-15, -15],
                        mode='lines',
                        hoverinfo='skip',
                        line = dict(
                            color='black',
                            width=3,
                        ),
                    ),
                    row=2+k*3, col=1
                )
            if (tmp_strand == "-") & (ind == 0):
                fig.add_trace(
                    go.Scatter(
                        x=np.array([j[0]/100, j[1]/100, j[1]/100, j[0]/100]) - rep1_region[0] if (j[1]-j[0] < 300) else np.array([j[0]/100, j[0]/100+3, j[1]/100, j[1]/100, j[0]/100+3, j[0]/100]) - rep1_region[0],
                        y=[-15, -20, -10, -15] if (j[1]-j[0] < 300) else [-15, -20, -20, -10, -10, -15],
                        mode='lines',
                        hoverinfo='skip',
                        line = dict(
                            color='black',
                            width=1,
                        ),
                        fill="toself",
                    ),
                    row=2+k*3, col=1
                )
            elif (tmp_strand == "+") & (ind == len(tmp_exon)-1):
                fig.add_trace(
                    go.Scatter(
                        x=np.array([j[0]/100, j[1]/100, j[0]/100, j[0]/100]) - rep1_region[0] if (j[1]-j[0] < 300) else np.array([j[0]/100, j[1]/100-3, j[1]/100, j[1]/100-3, j[0]/100, j[0]/100]) - rep1_region[0],
                        y=[-20, -15, -10, -20] if (j[1]-j[0] < 300) else [-20, -20, -15, -10, -10, -20],
                        mode='lines',
                        hoverinfo='skip',
                        line = dict(
                            color='black',
                            width=1,
                        ),
                        fill="toself",
                    ),
                    row=2+k*3, col=1
                )
            else:
                fig.add_trace(
                    go.Scatter(
                        x=np.array([j[0]/100, j[0]/100, j[1]/100, j[1]/100, j[0]/100]) - rep1_region[0],
                        y=[-10, -20, -20, -10, -10],
                        mode='lines',
                        hoverinfo='skip',
                        line = dict(
                            color='black',
                            width=1,
                        ),
                        fill="toself",
                    ),
                    row=2+k*3, col=1
                )
        fig.add_annotation(
            x=tmp_pos/100 - rep1_region[0], y=-25,
            text=tmp_NLR,
            showarrow=False,
            row=2+k*3, col=1
        )
        if k == 0:
            fig.update_layout(
                yaxis1 = dict(
                    range=[0, max(40, rep1_depth.max(), rep2_depth.max())]
                ),
                yaxis2 = dict(
                    range=[-30, 0]
                ),
                xaxis2_showticklabels=True,
                xaxis2 = dict(
                    tickmode = 'array',
                    tickvals = np.arange(rep1_region.min() - rep1_region[0], rep1_region.max() - rep1_region[0], 50),
                    ticktext = np.arange(rep1_region.min()*100, rep1_region.max()*100, 5000)
                ),
                yaxis2_showticklabels=False,
                yaxis2_ticklen=0,
                yaxis2_showgrid=False,
            )
        elif k == 1:
            fig.update_layout(
                yaxis4 = dict(
                    range=[0, max(40, rep1_depth.max(), rep2_depth.max())]
                ),
                yaxis5 = dict(
                    range=[-30, 0]
                ),
                xaxis5_showticklabels=True,
                xaxis5 = dict(
                    tickmode = 'array',
                    tickvals = np.arange(rep1_region.min() - rep1_region[0], rep1_region.max() - rep1_region[0], 50),
                    ticktext = np.arange(rep1_region.min()*100, rep1_region.max()*100, 5000)
                ),
                yaxis5_showticklabels=False,
                yaxis5_ticklen=0,
                yaxis5_showgrid=False,
            )
        else:
            fig.update_layout(
                yaxis7 = dict(
                    range=[0, max(40, rep1_depth.max(), rep2_depth.max())]
                ),
                yaxis8 = dict(
                    range=[-30, 0]
                ),
                xaxis8_showticklabels=True,
                xaxis8 = dict(
                    tickmode = 'array',
                    tickvals = np.arange(rep1_region.min() - rep1_region[0], rep1_region.max() - rep1_region[0], 50),
                    ticktext = np.arange(rep1_region.min()*100, rep1_region.max()*100, 5000)
                ),
                yaxis8_showticklabels=False,
                yaxis8_ticklen=0,
                yaxis8_showgrid=False,
            )
    fig.update_layout(
        template="ggplot2",
        showlegend=False,
        margin=dict(t=10, b=10, l=10, r=10),
    )
    fig.update_yaxes(
        fixedrange=True,
    )

    return fig