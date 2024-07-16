# Import packages
from dash import Dash, html, dash_table, dcc, callback, Output, Input, State, Patch, ctx
import pandas as pd
import numpy as np
from Bio import Phylo
import plotly.express as px
import plotly.graph_objs as go
import dash_bootstrap_components as dbc
import ipywidgets as ipw
import pathlib

from modules.phylogenetic_tree import make_tree_figure
from modules.genome_view import make_gene_markers
from modules.additional_data import make_expression_figure, make_chippeak_figure

#path
# BASE_PATH = pathlib.Path('__file__').parent.resolve()
BASE_PATH = pathlib.Path(__file__).parent.resolve()
# DATA_PATH = BASE_PATH.joinpath("../data").resolve()
DATA_PATH = BASE_PATH.joinpath("data").resolve()

# read data
table_data = pd.read_csv(DATA_PATH.joinpath("test_table_data.csv"))

# get chromosomes
chromosomes = ["All"]
chromosomes.extend(table_data.Chr.unique())

# read expression data
RNA_seq_data = pd.read_csv(DATA_PATH.joinpath("sample_RNAseq_counts.csv"), index_col=0)

# read ortholog data
ortholog_data = pd.read_csv(DATA_PATH.joinpath("ortholog_data.csv"), index_col=0)

# radio items
histone_list = ["H3K4me3", "H3K9ac", "H3K27me3"]
chip_treats = ["ABA", "CS_seedlings", "CS_leaf", "CS_spikelet_Feekes10", "MeJA", "SA", "Not displayed"]

# make phylogenetic tree
tree_fig, tree_id_table, tree_ids_len = make_tree_figure(DATA_PATH.joinpath("CS_tree.xml"))

# merge table_data and tree_id_table
table_data = pd.merge(table_data, tree_id_table, on="NLR id")

# Initialize the app - incorporate a Dash Bootstrap theme
external_stylesheets = [dbc.themes.CERULEAN]
app = Dash(__name__, external_stylesheets=external_stylesheets, suppress_callback_exceptions=True)
server = app.server

# App layout
app.layout = dbc.Container([
    
    dbc.Row([
        dbc.Col([
            html.Div('Wheat NLRome', className="text-primary text-center fs-3"),
            html.Br(),
            html.P('Chromosome'),
            dcc.Dropdown(
                chromosomes,
                chromosomes[0],
                id='dropdown-buttons-chr'
            ),
            html.Br(),
            html.P('Search', className="fs-5"),
            html.Div([
                dcc.Input(id='input-NLR-id', type='text', placeholder="ex. TraesCS1A03G0009700"),
                html.Button(id='submit-NLR-id', n_clicks=0, children='Submit'),
                html.Div(id='under_button'),
            ]),
            html.Br(),
            html.Div('Additional data', className="fs-5"),
            html.P('Expression data'),
            dcc.Dropdown(
                RNA_seq_data.columns.values,
                [RNA_seq_data.columns.values[0], RNA_seq_data.columns.values[1]],
                id='dropdown-content-exp',
                multi=True,
            ),
            html.Br(),
            html.P('ChIP-seq Peak'),
            dcc.Dropdown(
                histone_list,
                histone_list[0],
                id='dropdown-content-chip1',
            ),
            dcc.Dropdown(
                chip_treats,
                chip_treats[-1],
                id='dropdown-content-chip2',
            ),
            html.Br(),
            html.P('Core NLRs'),
            html.Div('Orthologs in Tritucum/Aegilops species is displayed below.')
        ],
        id="left-column",
        width=2
        ),
        
        dbc.Col([
            dcc.Store(id='stored-range'),
            dbc.Row([
                html.B('NLR locations', className="text-center fs-4"),
                html.Div(
                    className="text-center",
                    children="Click and drag on the below plot to Zoom-in and double-click to Zoom-out completely!",
                ),
                html.Hr(),
                dbc.Col([
                    dcc.Loading(
                        [dcc.Graph(figure={}, id='location-graph')],
                        type="circle",
                    ),
                ]),
            ], id="location-card"),

            dbc.Row([
                dbc.Col([
                    html.B('Phylogenetic tree', className="fs-4"),
                    dcc.Graph(figure=tree_fig, id='phylogenetic-tree')
                ],
                width=4,
                id="phylogeny-card"), 
                dbc.Col([
                    dbc.Row(id="additional-data")
                ],
                width=8,
                id="additional-card"), 
            ]),
            
            dbc.Row([
                dbc.Col([
                    html.B("Ortholog NLRs", className="fs-4"),
                ],
                id="ortholog-card"),
            ]),
            
            dbc.Row([
                html.B('NLR list', className="text-center fs-4"),
                html.Div(id="NLR-list-table"),
            ]),
        ],
        width=10,
        id="right-column")
    ])
], fluid=False)


# Show red marker in phylogenetic tree triggered by hover NLRs
@callback(
    Output('phylogenetic-tree', 'figure'),
    Input('location-graph', 'hoverData'),
)
def display_select_NLR_in_tree(hover_info):
    if hover_info is not None:
        edge_id = table_data.loc[(table_data.Chr == hover_info['points'][0]['y']) & (table_data.Graph_pos == hover_info['points'][0]['x']), "clade_id"]
    else:
        edge_id = None
    patched_figure = Patch()
    if edge_id is not None:
        new_size = np.repeat(0, tree_ids_len)
        new_size[edge_id] = 10
        patched_figure["data"][2]["marker"]["size"] = new_size
    return patched_figure

# Show NLR locations in genome
@callback(
    Output('location-graph', 'figure'),
    Output('under_button', 'children'),
    Output('dropdown-buttons-chr', 'value'),
    Input('dropdown-buttons-chr', 'value'),
    Input('submit-NLR-id', 'n_clicks'),
    State('input-NLR-id', 'value'),
)
def update_graph(chr_chosen, n_clicks, NLR_id):
    under_button = None
    zoom_range = None
    if "submit-NLR-id" == ctx.triggered_id:
        tmp_table_data = table_data.loc[table_data["NLR id"] == NLR_id, :]
        if tmp_table_data.shape[0] > 0:
            chr_chosen = tmp_table_data.Chr.values[0]
            zoom_range = [tmp_table_data.Start.min()-50000, tmp_table_data.End.max()+50000]
        else:
            under_button = f"{NLR_id} is not found in NLR list"
    fig = make_gene_markers(table_data, chr_chosen, zoom_range)
    return fig, under_button, chr_chosen

# Show Additional info
@callback(
    Output('additional-data', 'children'),
    Output('stored-range', 'data'),
    Input('location-graph', 'relayoutData'),
    Input('dropdown-buttons-chr', 'value'),
    [Input(f'dropdown-content-{i}', 'value') for i in ["exp", "chip1", "chip2"]],
    Input('submit-NLR-id', 'n_clicks'),
    State('input-NLR-id', 'value'),
    Input('stored-range', 'data'),
    prevent_initial_call=True,
)
def update_additional_info(relayoutData, chr_chosen, exp_chosen, chip1_chosen, chip2_chosen, n_click, NLR_id, store_range):
    div = []
    if chr_chosen == "All":
        div.extend(
            [
                html.B("Additional data", className="fs-4"),
                html.Br(),
                html.P('Please select chromosome/scaffold/contig from the dropdown menu on the left.')
            ]
        )
    else:
        if "submit-NLR-id" == ctx.triggered_id:
            tmp_pos = table_data.loc[table_data["NLR id"] == NLR_id, ["Chr", "Start", "End"]].values[0]
            tmp_table = table_data[(table_data.Chr == tmp_pos[0]) & (table_data.Start >= tmp_pos[1]-50000) & (table_data.End <= tmp_pos[2]+50000)]
        elif (ctx.triggered_id == "dropdown-content-exp") | (ctx.triggered_id == "dropdown-content-chip1") | (ctx.triggered_id == "dropdown-content-chip2"):
            tmp_table = table_data[(table_data.Chr == chr_chosen) & (table_data.Start >= store_range["start"]) & (table_data.End <= store_range["end"])]
        elif ("xaxis.range[0]" not in relayoutData) | (ctx.triggered_id == "dropdown-buttons-chr"):
            tmp_table = table_data[table_data.Chr == chr_chosen]
        else:
            xmin = relayoutData["xaxis.range[0]"]
            xmax = relayoutData["xaxis.range[1]"]
            tmp_table = table_data[(table_data.Chr == chr_chosen) & (table_data.Start >= xmin) & (table_data.End <= xmax)]
        store_range = {"start": tmp_table.Start.min(), "end": tmp_table.End.max()}
        
        add_figs = []
        # Expression bar plots
        if len(exp_chosen) >= 1:
            add_figs.append(
                [
                    html.B("Expression data", className="fs-4"),
                    html.Div("TPM values"),
                    dcc.Graph(figure=make_expression_figure(RNA_seq_data, exp_chosen, tmp_table["NLR id"].values))
                ]
            )
        
        # Chip-seq peak plots
        if (chip2_chosen != "Not displayed") & (chip2_chosen != None):
            showed_NLR_num = tmp_table.shape[0] 
            if showed_NLR_num > 3:
                add_figs.append(
                    [
                        html.B("ChIP-seq Peak", className="fs-4"),
                        html.P('Please Zoom-in more! (show 1~3 NLRs)')
                    ]
                )
            elif showed_NLR_num >= 1:
                figs = make_chippeak_figure(chip1_chosen, chip2_chosen, tmp_table, DATA_PATH)
                add_figs.append(
                    [
                        html.B("ChIP-seq Peak", className="fs-4"),
                        html.Div("Reads per genome coverage"),
                        dbc.Row(
                            [
                                dcc.Graph(figure=figs),
                            ],
                        ),
                    ]
                )
            else:
                add_figs.append(
                    [
                        html.B("ChIP-seq Peak", className="fs-4"),
                        dbc.Row(
                            [
                                html.P('No NLRs')
                            ]
                        ),
                    ]
                )
        div = [dbc.Col(each_fig, width=12//len(add_figs)) for each_fig in add_figs]
    return div, store_range

# Show other information of viewd NLRs
@callback(
    Output('NLR-list-table', 'children'),
    Output('ortholog-card', 'children'),
    Input('location-graph', 'relayoutData'),
    Input('dropdown-buttons-chr', 'value'),
    Input('submit-NLR-id', 'n_clicks'),
    State('input-NLR-id', 'value'),
    prevent_initial_call=True
)
def update_data(relayoutData, chr_chosen, n_click, NLR_id):
    if chr_chosen == "All":
        tmp_df = table_data
    else:
        if "submit-NLR-id" == ctx.triggered_id:
            tmp_pos = table_data.loc[table_data["NLR id"] == NLR_id, ["Chr", "Start", "End"]].values[0]
            tmp_df = table_data[(table_data.Chr == tmp_pos[0]) & (table_data.Start >= tmp_pos[1]-50000) & (table_data.End <= tmp_pos[2]+50000)]
        elif ("xaxis.range[0]" not in relayoutData) | (ctx.triggered_id == "dropdown-buttons-chr"):
            tmp_df = table_data[table_data.Chr == chr_chosen]
        else:
            xmin = relayoutData["xaxis.range[0]"]
            xmax = relayoutData["xaxis.range[1]"]
            tmp_df = table_data[(table_data.Chr == chr_chosen) & (table_data.Start >= xmin) & (table_data.End <= xmax)]
    
    table_div = html.Div(
        [
            dash_table.DataTable(data=tmp_df.loc[:, ["NLR id","Chr","Start","End","Strand","Domain"]].to_dict('records'),
                                 page_size=10,
                                 style_header={
                                    'backgroundColor': 'rgb(30, 30, 30)',
                                    'color': 'white'
                                 },
                                 style_data={
                                    'backgroundColor': 'rgb(50, 50, 50)',
                                    'color': 'white'
                                 },
                                 style_table={'overflowX': 'auto'},
                                 id="table")
        ]
    )
    
    tmp_NLR_ids = tmp_df["NLR id"].values
    if len(tmp_NLR_ids) > 20:
        ortholog_div = [
            html.B("Ortholog NLRs", className="fs-4"),
            html.P('Please select chromosome/scaffold/contig & Zoom-in more! (show ~20 NLRs)')
        ]
    else:
        ortholog_fig = make_ortholog_figure(ortholog_data, tmp_df["NLR id"].values)
        ortholog_div = [
            html.B("Ortholog NLRs", className="fs-4"),
            dcc.Graph(figure=ortholog_fig, id='ortholog-figure')
        ]
    return table_div, ortholog_div


# Run the app
if __name__ == '__main__':
    app.run_server(debug=False)
    # app.run(debug=True, jupyter_mode="external")
    # app.run(debug=True)
