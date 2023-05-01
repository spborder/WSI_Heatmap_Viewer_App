"""
Initializing FUSION web application for both local and AWS/Web deployment

This file should contain codes to generate:
    - layouts for each page
        - initial visualization layout with default dataset?
    - information dictionaries for available datasets

"""

import os
import sys

import json
import geojson

import pandas as pd
import numpy as np
from glob import glob
import lxml.etree as ET
from geojson import Feature, dump 
import uuid
import zipfile
import shutil

import plotly.express as px
import plotly.graph_objects as go

from dash import dcc, ctx, Dash, dash_table
import dash_bootstrap_components as dbc
import dash_cytoscape as cyto
import dash_leaflet as dl

from dash_extensions.enrich import html
from dash_extensions.javascript import arrow_function

from dataclasses import dataclass, field
from typing import Callable, List, Union
from dash.dependencies import handle_callback_args
from dash.dependencies import Input, Output, State


class LayoutHandler:
    def __init__(self,
                 verbose = False):
        
        self.verbose = verbose

        self.validation_layout = []
        self.layout_dict = {}
        self.description_dict = {}

        self.gen_welcome_layout()
        self.gen_uploader_layout()

    def gen_vis_layout(self,cell_types, center_point, map_dict, spot_dict):

        # Main visualization layout, used in initialization and when switching to the viewer

        # Description and instructions card
        vis_description = [
            html.P('FUSION was designed by the members of the CMI Lab at the University of Florida in collaboration with HuBMAP'),
            html.Hr(),
            html.P('We hope that this tool provides users with an immersive visualization method for understanding the roles of specific cell types in combination with different functional tissue units'),
            html.Hr(),
            html.P('As this tool is still under active development, we welcome any and all feedback. Use the "User Survey" link above to provide comments. Thanks!'),
            html.Hr(),
            html.P('Happy fusing!')         
        ]

        # Slide selection
        """
        slide_select = dbc.Card(
            id = 'slide-select-card',
            children = [
                dbc.CardHeader("Select case from dropdown menu"),
                dbc.CardBody([
                    dbc.Row([
                        dbc.Col(
                            html.Div(
                                id = "slide_select-label",
                                children = [
                                    html.P("Available cases: ")
                                ]
                            ),md=4
                        ),
                        dbc.Col(
                            html.Div(
                                dcc.Dropdown(
                                    slides_available,
                                    slides_available[0],
                                    id = 'slide-select'
                                )
                            ), md=8
                        )
                    ])
                ])
            ],style={'marginBottom':'20px'}
        )
        """

        # View of WSI

        self.initial_overlays = [
            dl.Overlay(
                dl.LayerGroup(
                    dl.GeoJSON(data = map_dict['FTUs'][struct]['geojson'], id = map_dict['FTUs'][struct]['id'], options = dict(color = map_dict['FTUs'][struct]['color']),
                        hoverStyle = arrow_function(dict(weight=5, color = map_dict['FTUs'][struct]['hover_color'], dashArray = '')))),
                name = struct, checked = True, id = struct)
        for struct in map_dict['FTUs']
        ] 

        self.initial_overlays+= [
            dl.Overlay(
                dl.LayerGroup(
                    dl.GeoJSON(data = spot_dict['geojson'], id = spot_dict['id'], options = dict(color = spot_dict['color']),
                        hoverStyle = arrow_function(dict(weight=5, color = spot_dict['hover_color'], dashArray = '')))),
                name = 'Spots', checked = False, id = 'Spots')
        ]

        wsi_view = dbc.Card([
            dbc.CardHeader('Whole Slide Image Viewer'),
            dbc.Row([
                html.Div(
                    dl.Map(center=center_point, zoom = 12, minZoom=11,
                           children = [
                                dl.TileLayer(url = map_dict['url'],id = 'slide-tile'),
                                dl.FeatureGroup(id='feature-group',children = [dl.EditControl(id='edit_control')]),
                                dl.LayerGroup(id='mini-label'),
                                html.Div(id='colorbar-div',children=[dl.Colorbar(id='map-colorbar')]),
                                dl.LayersControl(id='layer-control',children = self.initial_overlays)
                           ],
                           style = {'width':'100%','height':'80vh','margin':'auto','display':'inline-block'},
                           id = 'slide-map'
                    )
                )
            ]),
            dbc.Row([html.Div(id='current-hover')])
        ], style = {'marginBottom':'20px'})

        # Cell type proportions and cell state distributions
        roi_pie = dbc.Card([
            dbc.CardBody([
                html.Div(id = 'roi-pie-holder')
            ])
        ])

        # Stylesheet for cyto plot thingy
        cyto_style = [
            {
            'selector':'node',
            'style':{
                'label':'data(label)',
                'width':45,
                'height':45,
                'background-color':'white',
                'background-fit':'cover',
                'background-clip':'none',
                'background-image-opacity':1,
                #'border-opacity':1,
                #'border-width':2,
                'background-image':'data(url)',
                }
            },
            {
            'selector':'edge',
            'style':{
                'line-width':35,
                'line-color':'blue'
            }
            }
        ]

        # Cell card graphic and hierarchy
        cell_card_types = cell_types.copy()
        cell_card = dbc.Card([
            dbc.CardBody([
                dbc.Row([
                    dbc.Col([dbc.Label('Select a Cell Type to view the Cell Hierarchy')],md=4),
                    dbc.Col([
                        dcc.Dropdown(
                            cell_card_types,cell_card_types[0],id='cell-cards-drop'
                        )
                    ],md=8)
                ]),
                dbc.Row([
                    dbc.Col([
                        dbc.Label("Cell Graphic", html_for="cell-graphic"),
                        html.Img(
                            id = 'cell-graphic',
                            src = './assets/cell_graphics/default_cell_graphic.png',
                            height = '250px',
                            width = '500px'
                        )],md=6),
                    
                    dbc.Col([
                        dbc.Label("Cell Hierarchy",html_for="cell-hierarchy"),
                        cyto.Cytoscape(
                            id = 'cell-hierarchy',
                            layout={'name':'preset'},
                            style = {'width':'100%','height':'400px'},
                            minZoom = 0.5,
                            maxZoom = 3,
                            stylesheet=cyto_style,
                            elements = [
                                {'data': {'id': 'one', 'label': 'Node 1'}, 'position': {'x': 75, 'y': 75}},
                                {'data': {'id': 'two', 'label': 'Node 2'}, 'position': {'x': 200, 'y': 200}},
                                {'data': {'source': 'one', 'target': 'two'}}
                            ]),
                        html.Div(id='label-p'),
                        html.Div(id='id-p'),
                        html.Div(id='notes-p')
                        ],md=6)
                ],align='center')
            ]),
            dbc.CardFooter(
                dbc.Row([
                    dbc.Col([
                        html.Div([
                            dcc.Link('Derived from ASCT+B Kidney v1.2',href='https://docs.google.com/spreadsheets/d/1NMfu1bEGNFcTYTFT-jCao_lSbFD8n0ti630iIpRj-hw/edit#gid=949267305')
                        ])
                    ])
                ])
            )
        ])
        
        ftu_list = ['glomerulus','Tubules']
        plot_types = ['TSNE','UMAP']
        labels = ['Cluster','image_id']+cell_types.copy()
        # Cluster viewer tab
        cluster_card = dbc.Card([
            dbc.Row([
                dbc.Col([
                    dbc.Card(
                        id= 'plot-options',
                        children = [
                            dbc.CardHeader('Plot Options'),
                            dbc.CardBody([
                                dbc.Row([
                                    dbc.Col(dbc.Label('Functional Tissue Unit Type',html_for='ftu-select'),md=4),
                                    dbc.Col([
                                        html.Div(
                                            dcc.Dropdown(
                                                ftu_list,
                                                ftu_list[0],
                                                id='ftu-select'
                                            )
                                        )],md=8
                                    )
                                ]),
                                html.B(),
                                dbc.Row([
                                    dbc.Col(dbc.Label('Type of plot',html_for='plot-select'),md=4),
                                    dbc.Col([
                                        html.Div(
                                            dcc.Dropdown(
                                                plot_types,
                                                plot_types[0],
                                                id='plot-select'
                                            )
                                        )],md=8
                                    )
                                ]),
                                html.B(),
                                dbc.Row([
                                    dbc.Col(dbc.Label('Sample Labels',html_for='label-select'),md=4),
                                    dbc.Col([
                                        html.Div(
                                            dcc.Dropdown(
                                                labels,
                                                labels[0],
                                                id='label-select'
                                            )
                                        )],md=8
                                    )
                                ])
                            ])
                        ]
                    )
                ],md=12)
            ]),
            dbc.Row([
                dbc.Col([
                    html.Div(
                        dcc.Graph(id='cluster-graph',figure=go.Figure())
                    )
                ],md=6),
                dbc.Col([
                    dcc.Tabs([
                        dcc.Tab(
                            dbc.Card(
                                id = 'selected-image-card',
                                children = [
                                    dcc.Loading(
                                        id = 'loading-image',
                                        children = [
                                            dcc.Graph(id='selected-image',figure=go.Figure())
                                        ]
                                    )
                                ]
                            ),label='Selected Images'),
                        dcc.Tab(
                            dbc.Card(
                                id = 'selected-data-card',
                                children = [
                                    dcc.Loading(
                                        id='loading-data',
                                        children = [
                                            dbc.Row(
                                                children = [
                                                    dbc.Col(dcc.Graph(id='selected-cell-types',figure=go.Figure())),
                                                    dbc.Col(dcc.Graph(id='selected-cell-states',figure=go.Figure()))
                                                ]
                                            )
                                        ]
                                    )
                                ]
                            ),label='Selected Cell Data')
                    ]),
                    html.Div(id='selected-image-info')
                ],md=6),
            ],align='center')
        ])

        # Tools for selecting regions, transparency, and cells
        cell_types+=['Max Cell Type','Morphometrics Clusters','Area','Mesangial Area','Mesangial Fraction','Arterial Area','Luminal Fraction','Average TBM Thickness','Average Cell Thickness']
        
        # Converting the cell_types list into a dictionary to disable some
        disable_list = ['Morphometrics Clusters','Area','Mesangial Area','Mesangial Fraction','Arterial Area','Luminal Fraction','Average TBM Thickness','Average Cell Thickness']
        cell_types_list = []
        for c in cell_types:
            if c not in disable_list:
                cell_types_list.append({'label':c,'value':c,'disabled':False})
            else:
                cell_types_list.append({'label':c+' (In Progress)','value':c,'disabled':True})

        # Extracting data tab
        data_options = [
            {'label':'Annotations','value':'Annotations','disabled':False},
            {'label':'Cell Type and State','value':'Cell Type and State','disabled':False},
            {'label':'Slide Metadata','value':'Slide Metadata','disabled':False},
            {'label':'Selected FTUs and Metadata','value':'Selected FTUs and Metadata','disabled':True},
            {'label':'Manual ROIs','value':'Manual ROIs','disabled':True}
        ]
        extract_card = dbc.Card([
            dbc.CardBody([
                dbc.Row([
                    dbc.Label('Select data for download',html_for = 'data-select'),
                    dcc.Dropdown(data_options,placeholder = 'Select Data for Download',multi=True,id='data-select')
                ]),
                dbc.Row([dbc.Label('Download options',html_for = 'data-options')]),
                dbc.Row([
                    html.Div(id='data-options')
                ])
            ])
        ])
        

        mini_options = ['All Main Cell Types','Cell States for Current Cell Type','None']
        tools = [
            dbc.Card(
                id='tools-card',
                children=[
                    dbc.CardHeader("Tools"),
                    dbc.CardBody([
                        dbc.Row([
                            dbc.Col([
                                html.H6("Select Cell for Overlaid Heatmap Viewing",className="cell-select"),
                                html.Div(
                                    id = 'cell-select-div',
                                    children=[
                                        dcc.Dropdown(cell_types_list,cell_types_list[0]['value'],id='cell-drop')
                                    ]
                                )
                            ],md=6),
                            dbc.Col([
                                html.H6("Options for Overlaid Minicharts",className='mini-select'),
                                html.Div(
                                    id='mini-select-div',
                                    children=[
                                        dcc.Dropdown(mini_options,mini_options[0],id='mini-drop')
                                    ]
                                )
                            ])
                        ]),
                        html.Hr(),
                        dbc.Form([
                            dbc.Row([
                                dbc.Col([
                                    dbc.Label(
                                        "Adjust Transparency of Heatmap",
                                        html_for="vis-slider"
                                    ),
                                    dcc.Slider(
                                        id='vis-slider',
                                        min=0,
                                        max=100,
                                        step=10,
                                        value=50
                                    )
                                ])
                            ]),
                            dbc.Row([
                                dbc.Tabs([
                                    #dbc.Tab(thumbnail_select, label = "Thumbnail ROI"),
                                    dbc.Tab(roi_pie, label = "Cell Composition"),
                                    dbc.Tab(cell_card,label = "Cell Card"),
                                    dbc.Tab(cluster_card,label = 'Morphological Clustering'),
                                    dbc.Tab(extract_card,label = 'Download Data')
                                ])
                            ])
                        ])
                    ])
                ]
            )
        ]

        # Separately outputting the functional components of the application for later reference when switching pages
        vis_content = [
            dbc.Row(
                id="app-content",
                children=[
                    dbc.Col(wsi_view,md=6),
                    dbc.Col(tools,md=6)
                ],style={"height":"100vh"}
            )
        ]

        self.current_vis_layout = vis_content
        self.validation_layout.append(vis_content)
        self.layout_dict['vis'] = vis_content
        self.description_dict['vis'] = vis_description

    def gen_builder_layout(self, dataset_handler):

        # This builds the layout for the Dataset Builder functionality, 
        # allowing users to select which datasets/slides are incorporated into their 
        # current viewing instance.

        # Description and instructions card
        builder_description = [
            html.P('Happy fusing!')
        ]

        # Table containing information on each datset in dataset_handler.dataset_reference
        include_columns = ["name","organ","histology_type","stain","omics_type","description"]
        combined_dataset_dict = []
        for d_name in dataset_handler.dataset_names:
            specific_dict = dataset_handler.get_dataset(d_name)
            dataset_dict = {}
            for i in include_columns:
                dataset_dict[i] = specific_dict[i]

            # Adding extra info determined from nested keys (ftu, slide_info)
            if type(specific_dict['metadata'])==list:
                dataset_dict['metadata'] = ','.join(specific_dict['metadata'])
            else:
                dataset_dict['metadata'] = specific_dict['metadata']

            dataset_dict['annotation_type'] = specific_dict['ftu']['annotation_type']
            dataset_dict['FTUs'] = ','.join(list(specific_dict['ftu']['names'].keys()))
            dataset_dict['N_Slides'] = len(specific_dict['slide_info'])

            combined_dataset_dict.append(dataset_dict)
        
        dataset_df = pd.DataFrame.from_records(combined_dataset_dict)
        
        # Table with a bunch of filtering and tooltip info
        table_layout = html.Div([
            dash_table.DataTable(
                id = 'dataset-table',
                columns = [{'name':i,'id':i,'deletable':False,'selectable':True} for i in dataset_df],
                data = dataset_df.to_dict('records'),
                editable = False,
                filter_action='native',
                sort_action = 'native',
                sort_mode = 'multi',
                column_selectable = 'single',
                row_selectable = 'multi',
                row_deletable = False,
                selected_columns = [],
                selected_rows = [],
                page_action='native',
                page_current=0,
                page_size=10,
                style_cell = {
                    'overflow':'hidden',
                    'textOverflow':'ellipsis',
                    'maxWidth':0
                },
                tooltip_data = [
                    {
                        column: {'value':str(value),'type':'markdown'}
                        for column, value in row.items()
                    } for row in dataset_df.to_dict('records')
                ],
                tooltip_duration = None
            )
        ])

        
        builder_layout = [
                    html.H3('Select a Dataset to add slides to current session'),
                    html.Hr(),
                    table_layout,
                    html.H3('Select Slides to include in current session'),
                    html.Hr(),
                    dcc.Loading(html.Div(id='selected-dataset-slides')),
                    html.Hr(),
                    html.H3('Current Metadata'),
                    dcc.Loading(html.Div(id='slide-metadata-plots'))
                ]

        self.current_builder_layout = builder_layout
        self.validation_layout.append(builder_layout)
        self.layout_dict['dataset-builder'] = builder_layout
        self.description_dict['dataset-builder'] = builder_description

    def gen_uploader_layout(self):

        # This builds the layout for the Dataset Uploader functionality,
        # allowing users to upload their own data to be incorporated into the 
        # dataset builder or directly to the current viewing instance.

        # Description and instructions card
        uploader_description = [
            html.P('Happy fusing!')
        ]

        uploader_layout =[
                html.H1('Dataset Uploader')
            ]
        self.current_uploader_layout = uploader_layout
        self.validation_layout.append(uploader_layout)
        self.layout_dict['dataset-uploader'] = uploader_layout
        self.description_dict['dataset-uploader'] = uploader_description

    def gen_welcome_layout(self):

        # welcome layout after initialization and information and buttons to go to other areas

        # Description and instructions card
        welcome_description = [
            html.P('Happy fusing!')
        ]
        
        welcome_layout = [
                html.H1('Welcome to FUSION!')
            ]
        
        self.current_welcome_layout = welcome_layout
        self.validation_layout.append(welcome_layout)
        self.layout_dict['welcome'] = welcome_layout
        self.description_dict['welcome'] = welcome_description

    def gen_initial_layout(self,slide_names):

        # welcome layout after initialization and information and buttons to go to other areas

        # Header
        header = dbc.Navbar(
            dbc.Container([
                dbc.Row([
                    dbc.Col(html.Img(id='logo',src=('./assets/Lab_Logo_white.png'),height='75px'),md='auto'),
                    dbc.Col([
                        html.Div([
                            html.H3('FUSION'),
                            html.P('Functional Unit State Identification and Navigation with WSI')
                        ],id='app-title')
                    ],md=True,align='center')
                ],align='center'),
                dbc.Row([
                    dbc.Col([
                        dbc.NavbarToggler(id='navbar-toggler'),
                        dbc.Collapse(
                            dbc.Nav([
                                dbc.NavItem(
                                    dbc.Button(
                                        'User Survey',
                                        id = 'user-survey-button',
                                        outline = True,
                                        color = 'primary',
                                        href = ' https://ufl.qualtrics.com/jfe/form/SV_1A0CcKNLhTnFCHI',
                                        style = {'textTransform':'none'}
                                    )
                                ),
                                dbc.NavItem(
                                    dbc.Button(
                                        "Cell Cards",
                                        id='cell-cards-button',
                                        outline=True,
                                        color="primary",
                                        href="https://cellcards.org/index.php",
                                        style={"textTransform":"none"}
                                    )
                                ),
                                dbc.NavItem(
                                    dbc.Button(
                                        "Lab Website",
                                        id='lab-web-button',
                                        outline=True,
                                        color='primary',
                                        href='https://cmilab.nephrology.medicine.ufl.edu',
                                        style={"textTransform":"none"}
                                    )
                                )
                            ],navbar=True),
                        id="navbar-collapse",
                        navbar=True)
                    ],
                    md=2)
                    ],
                align='center')
                ], fluid=True),
            dark=True,
            color="dark",
            #sticky="top",
            style={'marginBottom':'20px'}
        )

        # Sidebar
        sider = html.Div([
            dbc.Offcanvas([
                html.Img(id='welcome-logo-side',src=('./assets/Lab_Logo.png'),height='280px',width='250px'),
                dbc.Nav([
                    dbc.NavLink('Welcome',href='/welcome',active='exact'),
                    dbc.NavLink('FUSION Visualizer',href='/vis',active='exact'),
                    dbc.NavLink('Dataset Builder',href='/dataset-builder',active='exact'),
                    dbc.NavLink('Dataset Uploader',href='/dataset-uploader',active='exact',disabled=True)
                ],vertical=True,pills=True)], id={'type':'sidebar-offcanvas','index':0},style={'background-color':"#f8f9fa"}
            )
        ])
        
        # Description and instructions card
        initial_description = [
            html.P('Happy fusing!')
        ]

        description = dbc.Card(
            children = [
                #dbc.CardHeader("Description and Instructions"),
                dbc.CardBody([
                    dbc.Button('Open Sidebar',id={'type':'sidebar-button','index':0},className='mb-3',color='primary',n_clicks=0,style={'marginRight':'5px'}),
                    dbc.Button("View/Hide Description",id={'type':'collapse-descrip','index':0},className='mb-3',color='primary',n_clicks=0,style={'marginLeft':'5px'}),
                    dbc.Collapse(
                        dbc.Row(
                            dbc.Col(
                                html.Div(
                                    id = 'descrip',
                                    children = initial_description,
                                    style={'fontSize':10}
                                )
                            )
                        ),id={'type':'collapse-content','index':0},is_open=False
                    )
                ])
            ],style={'marginBottom':'20px'}
        )

        # Slide select row (seeing if keeping it in every layout will let it be updated by actions in other pages)
        # Slide selection
        slide_select = dbc.Card(
            id = 'slide-select-card',
            children = [
                dbc.CardHeader("Select case from dropdown menu"),
                dbc.CardBody([
                    dbc.Row([
                        dbc.Col(
                            html.Div(
                                id = "slide_select-label",
                                children = [
                                    html.P("Available cases: ")
                                ]
                            ),md=4
                        ),
                        dbc.Col(
                            html.Div(
                                dcc.Dropdown(
                                    slide_names,
                                    slide_names[0],
                                    id = 'slide-select'
                                )
                            ), md=8
                        )
                    ])
                ])
            ],style={'marginBottom':'20px','display':'none'}
        )
        
        welcome_layout = html.Div([
            dcc.Location(id='url'),
            header,
            html.B(),
            dbc.Row(dbc.Col(html.Div(sider))),
            html.B(),
            dbc.Row(
                id = 'descrip-and-instruct',
                children = description
            ),
            dbc.Row(
                id = 'slide-select-row',
                children = slide_select
            ),
            dbc.Container([
                html.H1('Welcome to FUSION!'),
                ],fluid=True,id='container-content')
        ])

        self.current_initial_layout = welcome_layout
        self.validation_layout.append(welcome_layout)
        self.layout_dict['initial'] = welcome_layout
        self.description_dict['initial'] = initial_description


class DatasetHandler:
    def __init__(self,
                reference_path: str,
                verbose = False):
        
        self.reference_path = reference_path
        self.verbose = verbose

        # Reading the dataset_reference_file
        self.dataset_reference = json.load(open(self.reference_path))
        self.n_datasets = len(self.dataset_reference["datasets"])
        self.dataset_names = [self.dataset_reference["datasets"][i]["name"] for i in range(self.n_datasets)]
        self.dataset_sizes = [len(self.dataset_reference["datasets"][i]["slide_info"]) for i in range(self.n_datasets)]

    def get_dataset(self,dataset_name):

        # dataset_name = name of dataset to grab information from
        return self.dataset_reference["datasets"][self.dataset_names.index(dataset_name)]

    def get_slide_dataset(self,slide_name):
        
        dataset = None
        for d in range(len(self.dataset_names)):

            d_slides = [i['name'] for i in self.dataset_reference['datasets'][d]['slide_info']]

            if slide_name in d_slides:
                dataset = self.dataset_names[d]
                break
        
        if dataset is not None:
            return dataset
        else:
            raise ValueError


class DownloadHandler:
    def __init__(self,
                dataset_object,
                verbose = False):
        
        self.dataset_object = dataset_object
        self.verbose = verbose

        # Placeholder for lots of re-formatting and stuff

    def what_data(self,download_options):
        # Figuring out what data to return based on download_options used

        if download_options in ['Aperio XML','Histomics JSON','GeoJSON']:

            return 'annotations'
        
        elif download_options in ['CSV Files','Excel File','RDS File']:

            return 'cell'
        
        elif any([i in ['FTU Properties','Tissue Type','Omics Type','Slide Metadata','FTU Counts'] for i in download_options]):

            return 'metadata'
        
        elif 'man' in download_options:

            return 'manual_rois'
        
        else:

            return 'select_ftus'

    def zip_data(self,download_data_list):

        # Writing temporary directory 
        if not os.path.exists('./assets/FUSION_Download/'):
            os.makedirs('./assets/FUSION_Download/')

        output_file = './assets/FUSION_Download.zip'
        # Writing files to the temp directory
        for d in download_data_list:
            filename = d['filename']
            file_ext = filename.split('.')[-1]

            save_path = './assets/FUSION_Download/'+filename

            if file_ext in ['xml','json','geojson']:

                if file_ext == 'xml':
                    
                    # Writing xml data
                    with open(save_path,'w') as f:
                        f.write(d['content'])
                        f.close()
                else:
                    # Writing JSON data
                    with open(save_path,'w') as f:
                        dump(d['content'],f)

            elif file_ext == 'csv':

                d['content'].to_csv(save_path)
            
            elif file_ext == 'xlsx':
                # If it's supposed to be an excel file, include sheet name as a key
                with pd.ExcelWriter(save_path,mode='a') as writer:
                    d['content'].to_excel(writer,sheet_name=d['sheet'],engine='openpyxl')

        # Writing temporary data to a zip file
        with zipfile.ZipFile(output_file,'w', zipfile.ZIP_DEFLATED) as zip:
            for file in os.listdir('./assets/FUSION_Download/'):
                zip.write('./assets/FUSION_Download/'+file)

        try:
            shutil.rmtree('./assets/FUSION_Download/')
        except OSError as e:
            print(f'OSError removing FUSION_Download directory: {e.strerror}')

    def extract_annotations(self, slide, format):
        
        # Extracting annotations from the current slide object
        annotations = slide.geojson_ftus

        # Making output scale for the coordinates
        width_scale = slide.wsi_dims[0]/slide.slide_bounds[2]
        height_scale = slide.wsi_dims[1]/slide.slide_bounds[3]

        if format=='GeoJSON':
            
            save_name = slide.slide_name.replace('.'+slide.slide_ext,'.geojson')

            # Have to re-normalize coordinates and only save one or two properties
            final_ann = {'type':'FeatureCollection','features':[]}
            for f in annotations['features']:
                f_dict = {'type':'Feature','geometry':{'type':'Polygon','coordinates':[]}}

                scaled_coords = np.squeeze(np.array(f['geometry']['coordinates']))
                scaled_coords[:,0] *= height_scale
                scaled_coords[:,1] *= width_scale
                scaled_coords = scaled_coords.astype(int).tolist()
                f_dict['geometry']['coordinates'] = scaled_coords

                f_dict['properties'] = {'label':f['properties']['label'], 'structure':f['properties']['structure']}

                final_ann['features'].append(f_dict)

        elif format == 'Aperio XML':
            
            save_name = slide.slide_name.replace('.'+slide.slide_ext,'.xml')

            # Initializing xml file
            final_ann = ET.Element('Annotations')

            for l,name in enumerate(list(slide.ftus.keys())):
                # Need some kinda random color here
                ann = ET.SubElement(final_ann,'Annotation',attrib={'Type':'4','Visible':'1','ReadOnly':'0',
                                                            'Incremental':'0','LineColorReadOnly':'0',
                                                            'LineColor':'65280','Id':str(l),'NameReadOnly':'0'})
                
                regions = ET.SubElement(ann,'Regions')

                # Getting the features to add to this layer based on name
                layer_features = [i for i in annotations['features'] if i['properties']['structure']==name]

                for r_idx,f in enumerate(layer_features):

                    region = ET.SubElement(regions,'Region',attrib = {'NegativeROA':'0','ImageFocus':'-1',
                                                                    'DisplayId':'1','InputRegionId':'0',
                                                                    'Analyze':'0','Type':'0','Id':str(r_idx),
                                                                    'Text': f['properties']['label']})

                    verts = ET.SubElement(region,'Vertices')
                    scaled_coords = np.squeeze(np.array(f['geometry']['coordinates']))
                    scaled_coords[:,0] *= height_scale
                    scaled_coords[:,1] *= width_scale
                    scaled_coords = scaled_coords.astype(int).tolist()

                    for v in scaled_coords:
                        ET.SubElement(verts,'Vertex',attrib={'X':str(v[1]),'Y':str(v[0]),'Z':'0'})

                    ET.SubElement(verts,'Vertex',attrib={'X':str(scaled_coords[0][1]),'Y':str(scaled_coords[0][0]),'Z':'0'})

            final_ann = ET.tostring(final_ann,encoding='unicode',pretty_print=True)

        elif format == 'Histomics JSON':
            
            save_name = slide.slide_name.replace('.'+slide.slide_ext,'.json')

            # Following histomics JSON formatting
            final_ann = []

            for ftu_name in list(slide.ftus.keys()):

                output_dict = {'name':ftu_name,'attributes':{},'elements':[]}
                ftu_annotations = [i for i in annotations['features'] if i['properties']['structure']==ftu_name]

                for f in ftu_annotations:
                    scaled_coords = np.squeeze(np.array(f['geometry']['coordinates']))
                    scaled_coords[:,0] *= height_scale
                    scaled_coords[:,1] *= width_scale
                    scaled_coords = scaled_coords.astype(int).tolist()

                    struct_id = uuid.uuid4().hex[:24]
                    struct_dict = {
                        'type':'polyline',
                        'points':[i+[0] for i in scaled_coords],
                        'id':struct_id,
                        'closed':True,
                        'user':{
                            'label':f['properties']['label']
                        }
                    }
                    output_dict['elements'].append(struct_dict)

            final_ann.append(output_dict)

        return [{'filename': save_name, 'content':final_ann}]
    
    """
    def extract_metadata(self,slides, include_meta):
    
    """
    def extract_cell(self, slide, file_format):
        # Output here is a dictionary containing Main_Cell_Types and Cell_States for each FTU
        # Main_Cell_Types is a pd.DataFrame with a column for every cell type and an index of the FTU label
        # Cell states is a dictionary of pd.DataFrames with a column for every cell state and an index of the FTU label for each main cell type

        final_cell_info = {}
        for ftu in list(slide.ftus.keys()):

            ftu_annotations = [i for i in slide.geojson_ftus['features'] if i['properties']['structure']==ftu]
            ftu_main_cells = pd.DataFrame.from_records([i['properties']['Main_Cell_Types'] for i in ftu_annotations])
            ftu_main_cells.index = [i['properties']['label'] for i in ftu_annotations]
            
            cell_types = ftu_main_cells.columns
            ftu_cell_states = {}
            for c in cell_types:

                cell_state_df = pd.DataFrame.from_records([i['properties']['Cell_States'][c] for i in ftu_annotations])
                cell_state_df.index = [i['properties']['label'] for i in ftu_annotations]
                ftu_cell_states[c] = cell_state_df
            

            final_cell_info[ftu] = {'Main_Cell_Types':ftu_main_cells,'Cell_States':ftu_cell_states}

        # Formatting for downloads
        base_save_name = slide.slide_name.replace('.'+slide.slide_ext,'')
        download_data = []
        for ftu in list(final_cell_info.keys()):

            # Main cell types should just be one file
            print(f'file_format: {file_format}')
            if file_format=='Excel File':
                main_content = {'filename':base_save_name+'_Main_Cell_Types.xlsx','sheet':ftu,'content':final_cell_info[ftu]['Main_Cell_Types']}
            elif file_format=='CSV Files':
                main_content = {'filename':base_save_name+f'_Main_Cell_Types_{ftu}.csv','content':final_cell_info[ftu]['Main_Cell_Types']}
            else:
                print('Invalid format (RDS will come later)')
                main_content = {}
            
            download_data.append(main_content)

            # Cell state info
            for mc in final_cell_info[ftu]['Cell_States']:
                if file_format=='Excel File':
                    state_content = {'filename':base_save_name+f'_{ftu}_Cell_States.xlsx','sheet':mc.replace('/',''),'content':final_cell_info[ftu]['Cell_States'][mc]}
                elif file_format=='CSV Files':
                    state_content = {'filename':base_save_name+f'_{ftu}_{mc.replace("/","")}_Cell_States.csv','content':final_cell_info[ftu]['Cell_States'][mc]}
                else:
                    print('Invalid format (RDS will come later)')
                    state_content = {}

                download_data.append(state_content)

        return download_data

    """
    def extract_ftu(self, slide, data):
    

    def extract_manual(self, slide, data):
    """

@dataclass
class Callback:
    func: Callable
    outputs: Union[Output,List[Output]]
    inputs: Union[Input, List[Input]]
    states: Union[State, List[State]] = field(default_factory=list)
    kwargs: dict = field(default_factory=lambda: {"prevent_initial_call":False})

class CallbackManager:
    def __init__(self):
        self._callbacks = []

    def callback(self,*args,**kwargs):
        output, inputs, state, prevent_initial_call = handle_callback_args(
            args, kwargs
        )

        def wrapper(func):
            self._callbacks.append(Callback(func,
                                            output,
                                            inputs,
                                            state,
                                            {'prevent_initial_call':prevent_initial_call}))

        return wrapper
    
    def attach_to_app(self,app):
        for callback in self._callbacks:
            app.callback(
                callback.outputs,callback.inputs,callback.states,**callback.kwargs
            )(callback.func)















































