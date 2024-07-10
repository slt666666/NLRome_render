import pandas as pd
import plotly.graph_objs as go


def make_gene_markers_all_chr(table_data):
    fig = go.Figure()
    for tmp_chr in table_data.Chr.unique():
        tmp_table = table_data[table_data.Chr == tmp_chr]
        fig.add_traces(
            go.Scatter(
                x=[0, tmp_table.End.max()],
                y=[tmp_chr, tmp_chr],
                mode='lines',
                hoverinfo='skip',
                line=dict(color="gray")
            )
        )
        for i in range(tmp_table.shape[0]):
            start = tmp_table.iloc[i, :].Start
            end = tmp_table.iloc[i, :].End
            pos = tmp_table.iloc[i, :].Graph_pos
            NLR_id = tmp_table.iloc[i, :]["NLR id"]
            fig.add_traces(
                go.Scatter(
                    x=[start, end],
                    y=[tmp_chr, tmp_chr],
                    mode='lines',
                    hoverinfo='skip',
                    line=dict(
                        color="black",
                        width=8
                    )
                )
            )
            fig.add_traces(
                go.Scatter(
                    x=[pos],
                    y=[tmp_chr],
                    mode='markers',
                    hoverinfo='text',
                    hovertext=f'{NLR_id}<br>{tmp_table.iloc[i, :].Domain}',
                    name=NLR_id,
                    marker=dict(
                        symbol='triangle-right' if tmp_table.iloc[i, :].Strand == "+" else 'triangle-left',
                        color=tmp_table.iloc[i, :].Color,
                        size=14,
                        line=dict(
                            width=2,
                            color='black'
                        )
                    )
                )
            )
    fig.update_layout(
        showlegend=False,
        height=600,
        margin=dict(t=10, b=10, l=10, r=40),
    )
    fig.update_yaxes(
        categoryorder='array',
        categoryarray=table_data.Chr.unique()[::-1],
        fixedrange=True,
        tickfont=dict(size=15)
    )
    return fig

def make_gene_markers_each_chr(table_data, chr_chosen):
    fig = go.Figure()
    chr_table = table_data[table_data.Chr == chr_chosen]
    fig.add_traces(
        go.Scatter(
            x=[0, chr_table.End.max()],
            y=[chr_chosen, chr_chosen],
            mode='lines',
            hoverinfo='skip',
            line=dict(color="gray")
        )
    )
    for i in range(chr_table.shape[0]):
        start = chr_table.iloc[i, :].Start
        end = chr_table.iloc[i, :].End
        pos = chr_table.iloc[i, :].Graph_pos
        fig.add_traces(
            go.Scatter(
                x=[start, end],
                y=[chr_chosen, chr_chosen],
                mode='lines',
                hoverinfo='skip',
                line=dict(
                    color="black",
                    width=8
                )
            )
        )
        fig.add_trace(
            go.Scatter(
                x=[pos],
                y=[chr_chosen],
                hoverinfo='text',
                hovertext=f'{chr_table.iloc[i, :]["NLR id"]}<br>{chr_table.iloc[i, :].Domain}',
                name=chr_table.iloc[i, :]["NLR id"],
                mode='markers',
                marker=dict(
                    symbol='triangle-right' if chr_table.iloc[i, :].Strand == "+" else 'triangle-left',
                    color=chr_table.iloc[i, :].Color,
                    size=20,
                    line=dict(
                        width=3,
                        color='black'
                    )
                )
            ),
        )
    fig.update_layout(
        showlegend=False,
        margin=dict(t=10, b=20, l=10, r=10),
        height=160,
    )
    fig.update_yaxes(
        type='category',
        fixedrange=True,
        tickfont=dict(size=18)
    )
    return fig

def make_gene_markers(table_data, chr_chosen, zoom_range):
    if chr_chosen == "All":
        fig = make_gene_markers_all_chr(table_data)
    else:
        fig = make_gene_markers_each_chr(table_data, chr_chosen)
    if zoom_range != None:
        fig.update_layout(
            xaxis=dict(
                range=zoom_range
            )
        )
    return fig
