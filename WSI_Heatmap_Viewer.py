"""

Whole Slide Image heatmap generation and viewer for specific cell types using Ground Truth Spatial Transcriptomics

"""

import os
import sys
import pandas as pd
import numpy as np
import json

from glob import glob

from PIL import Image

try:
    import openslide
except:
    import tiffslide as openslide

import lxml.etree as ET

from tqdm import tqdm

import shapely
from shapely.geometry import Polygon, Point, shape
from skimage.draw import polygon
from skimage.transform import resize
import geojson
import random

import plotly.express as px
from plotly.subplots import make_subplots
import plotly.graph_objects as go
from matplotlib import cm

from dash import dcc, ctx, Dash
import dash_bootstrap_components as dbc
import dash_cytoscape as cyto
import dash_leaflet as dl
import dash_leaflet.express as dlx
from dash_extensions.javascript import assign, Namespace, arrow_function
from dash_extensions.enrich import DashProxy, html, Input, Output, MultiplexerTransform, State

from timeit import default_timer as timer


def gen_layout(cell_types,slides_available, center_point, map_dict, spot_dict):

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
        sticky="top",
        style={'marginBottom':'20px'}
    )

    # Description and instructions card
    description = dbc.Card(
        children = [
            #dbc.CardHeader("Description and Instructions"),
            dbc.CardBody([
                dbc.Button("Hide Instructions",id='collapse-descrip',className='mb-3',color='primary',n_clicks=0),
                dbc.Collapse(
                    dbc.Row(
                        dbc.Col(
                            html.Div(
                                id = 'descrip',
                                children = [
                                    html.P('Click within the thumbnail image to drop a square ROI over that point to view in full resolution.'),
                                    html.Hr(),
                                    html.P('Click within the Whole Slide Image Viewer for fine adjustments in current ROI'),
                                    html.Hr(),
                                    html.P('Select a specific cell type to adjust overlaid heatmap visualization'),
                                    html.Hr(),
                                    html.P('Heatmaps can be generated on a whole ROI basis (using overlapping patches), spot basis (showing raw spot cell type proportions), or on a per-Functional Tissue Unit (FTU) basis showing aggregated cell types for overlapping spots')
                                ],style={'fontSize':10}
                            )
                        )
                    ),id='collapse-content',is_open=False
                )
            ])
        ],style={'marginBottom':'20px'}
    )
    
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
    
    # View of WSI
    # Under development, full WSI view take as input a dict of {structure_name: {geojson:geojson, id:id, bounds_color:color, hover_color:color}}
    wsi_view = dbc.Card([
        dbc.CardHeader('Whole Slide Image Viewer'),
        dbc.Row([
            html.Div(
                dl.Map(center = center_point, zoom = 12, minZoom=11, children = [
                    dl.TileLayer(url = map_dict['url'], id = 'slide-tile'),
                    dl.LayerGroup(id='mini-label'),
                    dl.Colorbar(id='map-colorbar'),
                    dl.LayersControl(id = 'layer-control', children = 
                        [
                            dl.Overlay(
                                dl.LayerGroup(
                                    dl.GeoJSON(data = map_dict['FTUs'][struct]['geojson'], id = map_dict['FTUs'][struct]['id'], options = dict(color = map_dict['FTUs'][struct]['color']),
                                        hoverStyle = arrow_function(dict(weight=5, color = map_dict['FTUs'][struct]['hover_color'], dashArray = '')))),
                                name = struct, checked = True, id = 'Initial_'+struct)
                        for struct in map_dict['FTUs']
                        ] + 
                        [
                            dl.Overlay(
                                dl.LayerGroup(
                                    dl.GeoJSON(data = spot_dict['geojson'], id = spot_dict['id'], options = dict(color = spot_dict['color']),
                                        hoverStyle = arrow_function(dict(weight=5, color = spot_dict['hover_color'], dashArray = '')))),
                                name = 'Spots', checked = False, id = 'Initial_Spots')
                        ]
                    )
                ], style={'width': '85%', 'height': '80vh', 'margin': "auto", "display": "block"}, id = 'slide-map')
            )
        ]),
        dbc.Row([html.Div(id='current-hover')])
    ], style = {'marginBottom':'20px'})

    # Cell type proportions and cell state distributions
    roi_pie = dbc.Card([
        dbc.CardBody([
            dbc.Row([
                dbc.Col([
                    dbc.Label("ROI Cell Proportions", html_for="roi-pie"),
                    dcc.Graph(
                        id="roi-pie",
                        figure=go.Figure()
                    ),
                ],md=6),
                dbc.Col([
                    dbc.Label("Selected Cell States",html_for="state-bar"),
                    dcc.Graph(
                        id="state-bar",
                        figure=go.Figure()
                    )
                ],md=6)
            ])
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
    cell_card = dbc.Card([
        dbc.CardBody([
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
    labels = ['Cluster','image_id','Cell Type']
    # Cluster viewer tab
    cluster_card = dbc.Card([
        dbc.Row([
            dbc.Col([
                html.Div(
                    dcc.Graph(id='cluster-graph',figure=go.Figure())
                )
            ],md=4),
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
                                        dcc.Graph(id='selected-cell-data',figure=go.Figure())
                                    ]
                                )
                            ]
                        ),label='Selected Cell Data')
                ]),
                html.Div(id='selected-image-info')
            ],md=4),
            dbc.Col([
                dbc.Card(
                    id= 'plot-options',
                    children = [
                        dbc.CardHeader('Plot Options'),
                        dbc.CardBody([
                            dbc.Label('Functional Tissue Unit Type',html_for='ftu-select'),
                            dbc.Row([
                                dbc.Col(
                                    html.Div(
                                        dcc.Dropdown(
                                            ftu_list,
                                            ftu_list[0],
                                            id='ftu-select'
                                        )
                                    )
                                )
                            ]),
                            html.B(),
                            dbc.Label('Type of plot',html_for='plot-select'),
                            dbc.Row([
                                dbc.Col(
                                    html.Div(
                                        dcc.Dropdown(
                                            plot_types,
                                            plot_types[0],
                                            id='plot-select'
                                        )
                                    )
                                )
                            ]),
                            html.B(),
                            dbc.Label('Sample Labels',html_for='plot-select'),
                            dbc.Row([
                                dbc.Col(
                                    html.Div(
                                        dcc.Dropdown(
                                            labels,
                                            labels[0],
                                            id='label-select'
                                        )
                                    )
                                )
                            ])
                        ])
                    ]
                )
            ],md=4)
        ],align='center')
    ])

    # Tools for selecting regions, transparency, and cells
    cell_types+=['Max Cell Type','Cell Type States']
    mini_options = ['All Main Cell Types','Cell States for Current Cell Types','None']
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
                                    dcc.Dropdown(cell_types,cell_types[0],id='cell-drop')
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
                                dbc.Tab(cluster_card,label = 'Morphological Clustering')
                            ])
                        ])
                    ])
                ])
            ]
        )
    ]

    main_layout = html.Div([
        header,
        html.B(),
        dbc.Container([
            dbc.Row(
                id = 'descrip-and-instruct',
                children = [description]
            ),
            html.B(),
            dbc.Row(
                id = 'slide-select-row',
                children = [slide_select]
            ),
            html.B(),
            dbc.Row(
                id="app-content",
                children=[wsi_view]
            ),
            html.B(),
            dbc.Row(
                id='tab-info',
                children = [dbc.Col(tools,md=12)]
            )
        ],fluid=True)
    ])

    return main_layout


class WholeSlide:
    def __init__(self,
                image_url,
                slide_name,
                slide_info_dict,
                ftu_path,
                spot_path):
        
        self.image_url = image_url
        self.ftu_path = ftu_path
        self.spot_path = spot_path
        self.slide_info_dict = slide_info_dict
        self.slide_bounds = self.slide_info_dict['bounds']
        self.slide_name = slide_name
        self.slide_path = self.slide_info_dict['slide_path']

        # Efficiency test (grouping together ftus and spots into boxes with 50 structures in each)
        self.group_n = 50

        # Processing ftus and spots
        self.ftus = self.process_ftus()
        self.spots = self.process_spots()

    def process_spots(self):

        with open(self.spot_path) as f:
            geojson_polys = geojson.load(f)

        self.geojson_spots = geojson_polys

        # Parsing through geojson spots
        spot_polys = {}

        spot_groups = []
        group_count = 0
        group_bounds = []

        spot_polys['polygons'] = []
        spot_polys['barcodes'] = []
        spot_polys['main_counts'] = []
        spot_polys['cell_states'] = []

        for g in geojson_polys['features']:

            group_count+=1
            current_spot = shape(g['geometry'])
            current_bounds = list(current_spot.bounds)

            if not group_bounds == []:
                new_mins = [np.minimum(current_bounds[j],group_bounds[j]) for j in range(0,2)]
                new_maxs = [np.maximum(current_bounds[j],group_bounds[j]) for j in range(2,4)]

                group_bounds = new_mins+new_maxs
            else:
                group_bounds = current_bounds

            spot_polys['polygons'].append(current_spot)
            spot_polys['barcodes'].append(g['properties']['label'])
            spot_polys['main_counts'].append(g['properties']['Main_Cell_Types'])
            spot_polys['cell_states'].append(g['properties']['Cell_States'])

            if group_count==self.group_n-1:
                
                spot_groups.append({
                    'box':shapely.geometry.box(*group_bounds),
                    'polygons':spot_polys['polygons'],
                    'barcodes':spot_polys['barcodes'],
                    'main_counts':spot_polys['main_counts'],
                    'cell_states':spot_polys['cell_states']
                })

                spot_polys['polygons'] = []
                spot_polys['barcodes'] = []
                spot_polys['main_counts'] = []
                spot_polys['cell_states'] = []
                group_count = 0
                group_bounds = []

        spot_groups.append({
            'box':shapely.geometry.box(*group_bounds),
            'polygons':spot_polys['polygons'],
            'barcodes':spot_polys['barcodes'],
            'main_counts':spot_polys['main_counts'],
            'cell_states':spot_polys['cell_states']
        })

        return spot_groups
        
    def process_ftus(self):
        # Reading files in geojson format
        with open(self.ftu_path) as f:
            geojson_polys = geojson.load(f)

        self.geojson_ftus = geojson_polys

        # Parsing through info stored in geojson file
        ftu_polys = {}
        poly_ids = []
        for f in geojson_polys['features']:
            try:
                poly_ids.append(f['properties']['label'])
            except:
                continue

        poly_names = np.unique([i.split('_')[0] for i in poly_ids])
        self.ann_ids = {}

        # The difference here between the spots and the ftus is that the ftu groups will be a dictionary for each ftu
        # This opens it up to being controlled downstream (e.g. if turning off one FTU)
        ftu_groups = {}

        for p in poly_names:
            
            # Just initializing a dictionary format with empty lists
            self.ann_ids[p] = []

            poly_idx = [i for i in range(len(poly_ids)) if p in poly_ids[i]]
            p_features = [geojson_polys[i] for i in poly_idx]

            ftu_polys[p] = {}
            ftu_polys[p]['polygons'] = []
            ftu_polys[p]['barcodes'] = []
            ftu_polys[p]['main_counts'] = []
            ftu_polys[p]['cell_states'] = []

            ftu_groups[p] = []
            group_count = 0
            group_bounds = []

            # Iterating through each polygon in a given FTU
            for i in p_features:

                group_count+=1
                current_ftu = shape(i['geometry'])
                current_bounds = list(current_ftu.bounds)

                # Updating the groups bounding box
                if not group_bounds == []:
                    new_mins = [np.minimum(current_bounds[j],group_bounds[j]) for j in range(0,2)]
                    new_maxs = [np.maximum(current_bounds[j],group_bounds[j]) for j in range(2,4)]

                    group_bounds = new_mins+new_maxs
                else:
                    group_bounds = current_bounds

                # Adding info to the dictionary
                ftu_polys[p]['polygons'].append(current_ftu)
                ftu_polys[p]['barcodes'].append(i['properties']['label'])
                ftu_polys[p]['main_counts'].append(i['properties']['Main_Cell_Types'])
                ftu_polys[p]['cell_states'].append(i['properties']['Cell_States'])

                # Adding group info to the group list for a given ftu
                if group_count==self.group_n-1:
                    ftu_groups[p].append({
                        'box':shapely.geometry.box(*group_bounds),
                        'polygons':ftu_polys[p]['polygons'],
                        'barcodes':ftu_polys[p]['barcodes'],
                        'main_counts':ftu_polys[p]['main_counts'],
                        'cell_states':ftu_polys[p]['cell_states']
                        })

                    # resetting back to empty/0
                    ftu_polys[p]['polygons'] = []
                    ftu_polys[p]['barcodes'] = []
                    ftu_polys[p]['main_counts'] = []
                    ftu_polys[p]['cell_states'] = []
                    group_count = 0
                    group_bounds = []

            # Getting the last group which will have length < self.group_n
            ftu_groups[p].append({
                'box':shapely.geometry.box(*group_bounds),
                'polygons':ftu_polys[p]['polygons'],
                'barcodes':ftu_polys[p]['barcodes'],
                'main_counts':ftu_polys[p]['main_counts'],
                'cell_states':ftu_polys[p]['cell_states']
            })

        return ftu_groups

    def find_intersecting_spots(self,box_poly):
        
        # Find intersecting bounding boxes for a given box poly (each box will contain up to self.group_n spots)
        intersecting_groups = [i for i in range(0,len(self.spots)) if self.spots[i]['box'].intersects(box_poly)]

        # Searching only in intersecting groups for spots that intersect
        intersect_spots = []
        intersect_barcodes = []
        intersect_counts = []
        intersect_states = []
        for group_idx in intersecting_groups:

            intersect_idxes = [i for i in range(0,len(self.spots[group_idx]['polygons'])) if self.spots[group_idx]['polygons'][i].intersects(box_poly)]
            intersect_barcodes.extend([self.spots[group_idx]['barcodes'][i] for i in intersect_idxes])
            intersect_spots.extend([self.spots[group_idx]['polygons'][i] for i in intersect_idxes])
            intersect_counts.extend([self.spots[group_idx]['main_counts'][i] for i in intersect_idxes])
            intersect_states.extend([self.spots[group_idx]['cell_states'][i] for i in intersect_idxes])

        return {'polys':intersect_spots, 'barcodes':intersect_barcodes, 'main_counts':intersect_counts, 'states':intersect_states}    

    def find_intersecting_ftu(self,box_poly,ftu=None):

        intersect_barcodes = []
        intersect_ftus = []
        intersect_counts = []
        intersect_states = []

        if ftu is not None:
            ftu_list = [ftu]
        else:
            ftu_list = self.ann_ids.keys()

        for ann in ftu_list:
            
            # Which groups of FTUs intersect with the query box_poly
            group_intersect = [i for i in range(0,len(self.ftus[ann])) if self.ftus[ann][i]['box'].intersects(box_poly)]

            for g in group_intersect:
                group_polygons = self.ftus[ann][g]['polygons']
                group_barcodes = self.ftus[ann][g]['barcodes']
                group_counts = self.ftus[ann][g]['main_counts']
                group_states = self.ftus[ann][g]['cell_states']

                # Within group intersections
                intersect_idxes = [i for i in range(0,len(group_polygons)) if box_poly.intersects(group_polygons[i])]

                intersect_barcodes.extend([group_barcodes[i] for i in intersect_idxes])
                intersect_ftus.extend([group_polygons[i] for i in intersect_idxes])
                intersect_counts.extend([group_counts[i] for i in intersect_idxes])
                intersect_states.extend([group_states[i] for i in intersect_idxes])

        return {'polys':intersect_ftus, 'barcodes':intersect_barcodes, 'main_counts':intersect_counts, 'states':intersect_states}


class SlideHeatVis:
    def __init__(self,
                app,
                layout,
                wsi,
                cell_graphics_key,
                asct_b_table,
                cluster_metadata,
                slide_info_dict,
                run_type = None):
                
        # Setting some app-related things
        self.app = app
        self.app.title = "FUSION"
        self.app.layout = layout
        self.app._favicon = './assets/favicon.ico'

        self.run_type = run_type

        # clustering related properties (and also cell types, cell states, image_ids, etc.)
        self.metadata = cluster_metadata

        self.slide_list = list(slide_info_dict.keys())
        self.slide_paths = [slide_info_dict[i]['slide_path'] for i in self.slide_list]
        self.slide_info_dict = slide_info_dict
        self.wsi = wsi

        self.cell_graphics_key = json.load(open(cell_graphics_key))
        # Inverting the graphics key to get {'full_name':'abbreviation'}
        self.cell_names_key = {}
        for ct in self.cell_graphics_key:
            self.cell_names_key[self.cell_graphics_key[ct]['full']] = ct

        # Number of main cell types to include in pie-charts
        self.plot_cell_types_n = 5

        # ASCT+B table for cell hierarchy generation
        self.table_df = asct_b_table    

        # FTU settings
        self.ftus = list(self.wsi.ann_ids.keys())
        self.ftu_colors = {
            'Glomeruli':'#390191',
            'Tubules':'#e71d1d',
            'Arterioles':'#b6d7a8',
            'Spots':'#dffa00'
        }

        self.current_ftu_layers = ['Glomeruli','Tubules','Arterioles']
        self.pie_ftu = self.current_ftu_layers[-1]

        self.current_chart_coords = [0,0]

        # Initializing some parameters
        self.current_cell = 'PT'
        self.current_barcodes = []

        # Cell Hierarchy related properties
        self.node_cols = {
            'Anatomical Structure':{'abbrev':'AS','x_start':50,'y_start':75,
                                    'base_url':'https://www.ebi.ac.uk/ols/ontologies/uberon/terms?iri=http%3A%2F%2Fpurl.obolibrary.org%2Fobo%2FUBERON_'},
            'Cell Types':{'abbrev':'CT','x_start':250,'y_start':0,
                          'base_url':'https://www.ebi.ac.uk/ols/ontologies/uberon/terms?iri=http%3A%2F%2Fpurl.obolibrary.org%2Fobo%2FCL_'},
            'Genes':{'abbrev':'BGene','x_start':450,'y_start':75,
                     'base_url':'https://www.genenames.org/data/gene-symbol-report/#!/hgnc_id/'}
        }

        # Colormap settings (customize later)
        self.color_map = cm.get_cmap('jet')
        self.cell_vis_val = 0.5
        self.ftu_style_handle = assign("""function(feature,context){
            const {color_key,current_cell,fillOpacity,ftu_colors} = context.props.hideout;
            var cell_value = feature.properties.Main_Cell_Types[current_cell];
            if (cell_value==0){
                cell_value = 0.0;
            }

            if (cell_value==1){
                cell_value = 1.0;
            }
            const fillColor = color_key[cell_value];
            var style = {};
            style.fillColor = fillColor;
            style.fillOpacity = fillOpacity;

            return style;
            }
            """
        )

        self.app.callback(
            [Output('roi-pie','figure'),Output('state-bar','figure')],
            [Input('slide-map','zoom'),Input('slide-map','viewport')],
            State('slide-map','bounds')
        )(self.update_roi_pie)

        self.app.callback(
            [Output('cell-graphic','src'),Output('cell-hierarchy','elements'),
             Output('layer-control','children')],
             [Input('cell-drop','value'),Input('vis-slider','value')]
        )(self.update_cell)

        self.app.callback(
            [Output('state-bar','figure'),Output('roi-pie','figure')],
            Input('roi-pie','clickData'),
            prevent_initial_call=True
        )(self.update_state_bar)

        self.app.callback(
            Output('current-hover','children'),
            [Input('glom-bounds','hover_feature'),
             Input('spot-bounds','hover_feature'),
             Input('tub-bounds','hover_feature'),
             Input('art-bounds','hover_feature')],
             prevent_initial_call=True
        )(self.get_hover)

        self.app.callback(
            Output('mini-label','children'),
            [Input('glom-bounds','click_feature'),
            Input('spot-bounds','click_feature'),
            Input('tub-bounds','click_feature'),
            Input('art-bounds','click_feature')],
            prevent_initial_call=True
        )(self.get_click)

        self.app.callback(
            [Output('collapse-content','is_open'),
            Output('collapse-descrip','children')],
            [Input('collapse-descrip','n_clicks'),
            Input('collapse-descrip','children')],
            [State('collapse-content','is_open')],
            prevent_initial_call=True
        )(self.view_instructions)

        self.app.callback(
            [Output('slide-tile','url'), Output('layer-control','children'), Output('slide-map','center'),
             Output('roi-pie','figure'),Output('state-bar','figure')],
            Input('slide-select','value'),
            prevent_initial_call=True
        )(self.ingest_wsi)
        
        self.app.callback(
            [Output('label-p','children'),
            Output('id-p','children'),
            Output('notes-p','children')],
            Input('cell-hierarchy','tapNodeData'),
            prevent_initial_call=True
        )(self.get_cyto_data)

        self.app.callback(
            [Input('ftu-select','value'),
            Input('plot-select','value'),
            Input('label-select','value')],
            Output('cluster-graph','figure')
        )(self.update_graph)

        self.app.callback(
            [Input('cluster-graph','clickData'),
            Input('cluster-graph','selectedData')],
            Output('selected-image','figure'),
            prevent_initial_call=True
        )(self.update_selected)
       
        # Comment out this line when running on the web
        if self.run_type == 'local':
            self.app.run_server(debug=True,use_reloader=True,port=8000)

        elif self.run_type == 'AWS':
            self.app.run_server(host = '0.0.0.0',debug=False,use_reloader=False,port=8000)

    def view_instructions(self,n,text,is_open):
        if text == 'View Instructions':
            new_text = 'Hide Instructions'
        else:
            new_text = 'View Instructions'
        if n:
            return [not is_open,new_text]
        return [is_open,new_text]
    
    def update_roi_pie(self,zoom,viewport,bounds):

        # Making a box-poly from the bounds
        if len(bounds)==2:
            bounds_box = shapely.geometry.box(bounds[0][1],bounds[0][0],bounds[1][1],bounds[1][0])
        else:
            bounds_box = shapely.geometry.box(*bounds)

        # Getting a dictionary containing all the intersecting spots with this current ROI
        intersecting_ftus = {}
        if 'Spots' in self.current_ftu_layers:
            intersecting_spots = self.wsi.find_intersecting_spots(bounds_box)
            intersecting_ftus['Spots'] = intersecting_spots

        for ftu in self.current_ftu_layers:
            if not ftu=='Spots':
                intersecting_ftus[ftu] = self.wsi.find_intersecting_ftu(bounds_box,ftu)

        self.current_ftus = intersecting_ftus
        # Now we have main cell types, cell states, by ftu

        included_ftus = list(intersecting_ftus.keys())
        included_ftus = [i for i in included_ftus if len(intersecting_ftus[i]['polys'])>0]
        # Making subplots for each 
        spec_list = [[{'type':'domain','rowspan':len(included_ftus)-1}]]
        if len(included_ftus)>1:
            spec_list[0].extend([{'type':'domain'}])
        if len(included_ftus)>2:
            for inc_f in range(0,len(included_ftus)-2):
                spec_list.append([None,{'type':'domain'}])
        
        if len(included_ftus)>1:
            combined_pie = make_subplots(
                rows = len(included_ftus)-1, cols = 2,
                subplot_titles = included_ftus,
                specs = spec_list
            )
        else:
            combined_pie = go.Figure()

        # Iterating through intersecting_ftus and getting combined main cell proportions
        figs_to_include = []
        for f in included_ftus:
            
            counts_data = pd.DataFrame(intersecting_ftus[f]['main_counts']).sum(axis=0).to_frame()
            counts_data.columns = [f]

            # Normalizing to sum to 1
            counts_data[f] = counts_data[f]/counts_data[f].sum()
            # Only getting the top-5
            counts_data = counts_data.sort_values(by=f,ascending=False).iloc[0:self.plot_cell_types_n,:]
            counts_data = counts_data.reset_index()
            f_pie = px.pie(counts_data,values=f,names='index')
            figs_to_include.append(f_pie)

        if 'data' in figs_to_include[0]:
            for trace in range(len(figs_to_include[0]['data'])):
                combined_pie.append_trace(figs_to_include[0]['data'][trace],row=1,col=1)
        else:
            combined_pie.append_trace(figs_to_include[0],row=1,col=1)
        
        if len(included_ftus)>1:
            for i,figure in enumerate(figs_to_include[1:]):
                if 'data' in figure:
                    for trace in range(len(figure['data'])):
                        combined_pie.append_trace(figure['data'][trace],row=i+1,col=2)
                else:
                    combined_pie.append_trace(figure,row=i+1,col=2)

        # Picking cell + ftu for cell state proportions plot
        top_cell = counts_data['index'].tolist()[0]
        pct_states = pd.DataFrame([i[top_cell] for i in intersecting_ftus[f]['states']]).sum(axis=0).to_frame()
        pct_states = pct_states.reset_index()
        pct_states.columns = ['Cell State','Proportion']
        pct_states['Proportion'] = pct_states['Proportion']/pct_states['Proportion'].sum()

        state_bar = go.Figure(px.bar(pct_states,x='Cell State', y = 'Proportion', title = f'Cell State Proportions for {self.cell_graphics_key[top_cell]["full"]} in {f}'))

        return combined_pie, state_bar

    def make_new_pie_subplots(self,big_ftu):

        # Making the new focus-ftu the first in the list (probably another way to do this)
        included_ftus = list(self.current_ftus.keys())
        included_ftus = [i for i in included_ftus if len(self.current_ftus[i]['polys'])>0 and not i==big_ftu]
        included_ftus = [big_ftu]+included_ftus
        self.current_ftu_layers = [i for i in self.current_ftu_layers if not i == big_ftu]
        self.current_ftu_layers = [big_ftu]+self.current_ftu_layers

        # Re-generating pie-chart subplots with new focus-ftu
        spec_list = [[{'type':'domain','rowspan':len(included_ftus)-1}]]
        if len(included_ftus)>1:
            spec_list[0].extend([{'type':'domain'}])
        if len(included_ftus)>2:
            for inc_f in range(0,len(included_ftus)-2):
                spec_list.append([None,{'type':'domain'}])
        
        if len(included_ftus)>1:
            combined_pie = make_subplots(
                rows = len(included_ftus)-1, cols = 2,
                subplot_titles = included_ftus,
                specs = spec_list
            )
        else:
            combined_pie = go.Figure()

        # Iterating through intersecting_ftus and getting combined main cell proportions
        figs_to_include = []
        for f in included_ftus:
            
            counts_data = pd.DataFrame(self.current_ftus[f]['main_counts']).sum(axis=0).to_frame()
            counts_data.columns = [f]

            # Normalizing to sum to 1
            counts_data[f] = counts_data[f]/counts_data[f].sum()
            # Only getting the top-5
            counts_data = counts_data.sort_values(by=f,ascending=False).iloc[0:self.plot_cell_types_n,:]
            counts_data = counts_data.reset_index()
            f_pie = px.pie(counts_data,values=f,names='index')
            figs_to_include.append(f_pie)

        
        if 'data' in figs_to_include[0]:
            for trace in range(len(figs_to_include[0]['data'])):
                combined_pie.append_trace(figs_to_include[0]['data'][trace],row=1,col=1)
        else:
            combined_pie.append_trace(figs_to_include[0],row=1,col=1)
        
        if len(included_ftus)>1:
            for i,figure in enumerate(figs_to_include[1:]):
                if 'data' in figure:
                    for trace in range(len(figure['data'])):
                        combined_pie.append_trace(figure['data'][trace],row=i+1,col=2)
                else:
                    combined_pie.append_trace(figure,row=i+1,col=2)

        return combined_pie

    def update_state_bar(self,cell_click):
        
        if not cell_click is None:
            self.pie_cell = cell_click['points'][0]['label']

            if not self.current_ftu_layers[cell_click['points'][0]['curveNumber']] == self.pie_ftu:
                # Re-generating pie-charts with the new FTU as the left big pie-chart
                self.pie_ftu = self.current_ftu_layers[cell_click['points'][0]['curveNumber']]
                new_pie_chart = self.make_new_pie_subplots(self.pie_ftu)
                self.current_pie_chart = new_pie_chart

            else:
                new_pie_chart = self.current_pie_chart

        else:
            self.pie_cell = self.current_cell
            self.pie_ftu = self.current_ftu_layers[-1]

            new_pie_chart = self.make_new_pie_subplots(self.pie_ftu)
            self.current_pie_chart = new_pie_chart

        pct_states = pd.DataFrame([i[self.pie_cell] for i in self.current_ftus[self.pie_ftu]['states']]).sum(axis=0).to_frame()
        pct_states = pct_states.reset_index()
        pct_states.columns = ['Cell State', 'Proportion']
        pct_states['Proportion'] = pct_states['Proportion']/pct_states['Proportion'].sum()

        state_bar = go.Figure(px.bar(pct_states,x='Cell State', y = 'Proportion', title = f'Cell State Proportions for {self.cell_graphics_key[self.pie_cell]["full"]} in {self.pie_ftu}'))

        return state_bar, new_pie_chart
    
    def update_hex_color_key(self):

        # Iterate through all structures (and spots) in current wsi,
        # concatenate all of their proportions of specific cell types together
        # scale it with self.color_map (make sure to multiply by 255 after)
        # convert uint8 RGB colors to hex
        # create look-up table for original value --> hex color
        # add that as get_color() function in style dict (fillColor) along with fillOpacity
        raw_values_list = []
        id_list = []
        # iterating through current ftus
        for f in self.wsi.ftus:
            for g in self.wsi.ftus[f]:
                # get main counts for this ftu
                ftu_counts = pd.DataFrame(g['main_counts'])[self.current_cell].tolist()
                raw_values_list.extend(ftu_counts)

        for g in self.wsi.spots:
            spot_counts = pd.DataFrame(g['main_counts'])[self.current_cell].tolist()
            raw_values_list.extend(spot_counts)

        raw_values_list = np.unique(raw_values_list)
        # Converting to RGB
        rgb_values = np.uint8(255*self.color_map(np.uint8(255*raw_values_list)))[:,0:3]
        hex_list = []
        for row in range(rgb_values.shape[0]):
            hex_list.append('#'+"%02x%02x%02x" % (rgb_values[row,0],rgb_values[row,1],rgb_values[row,2]))

        self.hex_color_key = {i:j for i,j in zip(raw_values_list,hex_list)}

    def update_cell(self,cell_val,vis_val):
        
        # Updating current cell prop
        self.current_cell = self.cell_names_key[cell_val]
        self.cell_vis_val = vis_val/100

        self.update_hex_color_key()

        vis_val = vis_val/100

        # Changing fill and fill-opacity properties for structures and adding that as a property

        # Modifying ftu and spot geojson to add fill color and opacity
        for f in self.wsi.geojson_ftus['features']:
            
            cell_pct = f['properties']['Main_Cell_Types'][self.current_cell]
            hex_color = self.hex_color_key[cell_pct]
            f['properties']['fillColor'] = hex_color
            f['properties']['fillOpacity'] = vis_val

        map_dict = {
            'url':self.wsi.image_url,
            'FTUs':{
                'Glomeruli': {
                    'geojson':{'type':'FeatureCollection', 'features': [i for i in self.wsi.geojson_ftus['features'] if i['properties']['structure']=='Glomeruli']},
                    'id': 'glom-bounds',
                    'color': '#390191',
                    'hover_color':'#666'
                },
                'Tubules': {
                    'geojson':{'type':'FeatureCollection', 'features': [i for i in self.wsi.geojson_ftus['features'] if i['properties']['structure']=='Tubules']},
                    'id':'tub-bounds',
                    'color': '#e71d1d',
                    'hover_color': '#ff0b0a'
                },
                'Arterioles': {
                    'geojson':{'type':'FeatureCollection', 'features': [i for i in self.wsi.geojson_ftus['features'] if i['properties']['structure']=='Arterioles']},
                    'id':'art-bounds',
                    'color': '#b6d7a8',
                    'hover_color': '#50f207'
                }
            }
        }

        spot_dict = {
            'geojson':self.wsi.geojson_spots,
            'id': 'spot-bounds',
            'color': '#dffa00',
            'hover_color':'#9caf00'
        }

        new_children = [
            dl.Overlay(
                dl.LayerGroup(
                    dl.GeoJSON(data = map_dict['FTUs'][struct]['geojson'], id = map_dict['FTUs'][struct]['id'], options = dict(style=self.ftu_style_handle),
                            hideout = dict(color_key = self.hex_color_key, current_cell = self.current_cell, fillOpacity = vis_val, ftu_colors=self.ftu_colors),
                            hoverStyle = arrow_function(dict(weight=5, color = map_dict['FTUs'][struct]['hover_color'], dashArray = '')))),
                    name = struct, checked = True, id = self.wsi.slide_info_dict['key_name']+'_'+struct)
            for struct in map_dict['FTUs']
            ]
        
        new_children += [
            dl.Overlay(
                dl.LayerGroup(
                    dl.GeoJSON(data = spot_dict['geojson'], id = spot_dict['id'], options = dict(style = self.ftu_style_handle),
                            hideout = dict(color_key = self.hex_color_key, current_cell = self.current_cell, fillOpacity = vis_val, ftu_colors= self.ftu_colors),
                            hoverStyle = arrow_function(dict(weight=5, color = spot_dict['hover_color'], dashArray = '')))),
                    name = 'Spots', checked = False, id = self.wsi.slide_info_dict['key_name']+'_Spots')
            ]
    
        # Loading the cell-graphic and hierarchy image
        cell_graphic = self.cell_graphics_key[self.current_cell]['graphic']
        cell_hierarchy = self.gen_cyto()

        return cell_graphic, cell_hierarchy, new_children

    def get_hover(self,glom_hover,spot_hover,tub_hover,art_hover):

        hover_text = ''
        if 'glom-bounds.hover_feature' in ctx.triggered_prop_ids:
            if not glom_hover is None:
                hover_text = f'Glomerulus, {self.current_cell}: {round(glom_hover["properties"]["Main_Cell_Types"][self.current_cell],3)}'
        
        if 'spot-bounds.hover_feature' in ctx.triggered_prop_ids:
            if not spot_hover is None:
                hover_text = f'Spot, {self.current_cell}: {round(spot_hover["properties"]["Main_Cell_Types"][self.current_cell],3)}'

        if 'tub-bounds.hover_feature' in ctx.triggered_prop_ids:
            if not tub_hover is None:
                hover_text = f'Tubule, {self.current_cell}: {round(tub_hover["properties"]["Main_Cell_Types"][self.current_cell],3)}'

        if 'art-bounds.hover_feature' in ctx.triggered_prop_ids:
            if not art_hover is None:
                hover_text = f'Arteriole, {self.current_cell}: {round(art_hover["properties"]["Main_Cell_Types"][self.current_cell],3)}'

        return hover_text
    
    def get_click(self,glom_click,spot_click,tub_click,art_click):

        if 'glom-bounds.click_feature' in ctx.triggered_prop_ids:
            if not glom_click is None:
                chart_coords = np.mean(np.squeeze(glom_click['geometry']['coordinates']),axis=0)
                chart_dict_data = glom_click['properties']['Main_Cell_Types']
                chart_labels = list(chart_dict_data.keys())
                chart_data = [chart_dict_data[j] for j in chart_labels]

                if not all([i==j for i,j in zip(chart_coords,self.current_chart_coords)]):
                    self.current_chart_coords = chart_coords

                    mini_pie_chart = dl.Minichart(
                        data = chart_data, 
                        labels = chart_labels, 
                        lat=chart_coords[1],
                        lon=chart_coords[0],
                        height=100,
                        width=100,
                        labelMinSize=2,
                        type='pie',id=f'glom_pie_click{random.randint(0,1000)}')

                    return [mini_pie_chart]

        if 'tub-bounds.click_feature' in ctx.triggered_prop_ids:
            if not tub_click is None:
                chart_coords = np.mean(np.squeeze(tub_click['geometry']['coordinates']),axis=0)
                chart_dict_data = tub_click['properties']['Main_Cell_Types']
                chart_labels = list(chart_dict_data.keys())
                chart_data = [chart_dict_data[j] for j in chart_labels]

                if not all([i==j for i,j in zip(chart_coords,self.current_chart_coords)]):
                    self.current_chart_coords = chart_coords

                    mini_pie_chart = dl.Minichart(
                        data = chart_data, 
                        labels = chart_labels, 
                        lat=chart_coords[1],
                        lon=chart_coords[0],
                        height=100,
                        width=100,
                        labelMinSize=2,
                        type='pie',id=f'tub_pie_click{random.randint(0,1000)}')

                    return [mini_pie_chart]

        if 'spot-bounds.click_feature' in ctx.triggered_prop_ids:
            if not spot_click is None:
                chart_coords = np.mean(np.squeeze(spot_click['geometry']['coordinates']),axis=0)
                chart_dict_data = spot_click['properties']['Main_Cell_Types']
                chart_labels = list(chart_dict_data.keys())
                chart_data = [chart_dict_data[j] for j in chart_labels]

                if not all([i==j for i,j in zip(chart_coords,self.current_chart_coords)]):
                    self.current_chart_coords = chart_coords

                    mini_pie_chart = dl.Minichart(
                        data = chart_data, 
                        labels = chart_labels, 
                        lat=chart_coords[1],
                        lon=chart_coords[0],
                        height=100,
                        width=100,
                        labelMinSize=2,
                        type='pie',id=f'spot_pie_click{random.randint(0,1000)}')

                    return [mini_pie_chart]

        if 'art-bounds.click_feature' in ctx.triggered_prop_ids:
            if not art_click is None:
                chart_coords = np.mean(np.squeeze(art_click['geometry']['coordinates']),axis=0)
                chart_dict_data = art_click['properties']['Main_Cell_Types']
                chart_labels = list(chart_dict_data.keys())
                chart_data = [chart_dict_data[j] for j in chart_labels]

                if not all([i==j for i,j in zip(chart_coords,self.current_chart_coords)]):
                    self.current_chart_coords = chart_coords

                    mini_pie_chart = dl.Minichart(
                        data = chart_data, 
                        labels = chart_labels, 
                        lat=chart_coords[1],
                        lon=chart_coords[0],
                        height=100,
                        width=100,
                        labelMinSize=2,
                        type='pie',id=f'art_pie_click{random.randint(0,1000)}')

                    return [mini_pie_chart]
             
    def gen_cyto(self):

        cyto_elements = []

        # Getting cell sub-types under that main cell
        cell_subtypes = self.cell_graphics_key[self.current_cell]['subtypes']

        # Getting all the rows that contain these sub-types
        table_data = self.table_df.dropna(subset = ['CT/1/ABBR'])
        cell_data = table_data[table_data['CT/1/ABBR'].isin(cell_subtypes)]

        # cell type
        cyto_elements.append(
            {'data':{'id':'Main_Cell',
                     'label':self.current_cell,
                     'url':'./assets/cell.png'},
            'classes': 'CT',
            'position':{'x':self.node_cols['Cell Types']['x_start'],'y':self.node_cols['Cell Types']['y_start']},
                     }
        )

        # Getting the anatomical structures for this cell type
        an_structs = cell_data.filter(regex=self.node_cols['Anatomical Structure']['abbrev']).dropna(axis=1)

        an_start_y = self.node_cols['Anatomical Structure']['y_start']
        col_vals = an_structs.columns.values.tolist()
        col_vals = [i for i in col_vals if 'LABEL' in i]

        for idx,col in enumerate(col_vals):
            cyto_elements.append(
                {'data':{'id':col,
                         'label':an_structs[col].tolist()[0],
                         'url':'./assets/kidney.png'},
                'classes':'AS',
                'position':{'x':self.node_cols['Anatomical Structure']['x_start'],'y':an_start_y}
                         }
            )
            
            if idx>0:
                cyto_elements.append(
                    {'data':{'source':col_vals[idx-1],'target':col}}
                )
            an_start_y+=75
        
        last_struct = col
        cyto_elements.append(
            {'data':{'source':last_struct,'target':'Main_Cell'}}
        )
        
        cell_start_y = self.node_cols['Cell Types']['y_start']
        gene_start_y = self.node_cols['Genes']['y_start']
        for idx_1,c in enumerate(cell_subtypes):

            matching_rows = table_data[table_data['CT/1/ABBR'].str.match(c)]

            if not matching_rows.empty:
                cell_start_y+=75

                cyto_elements.append(
                    {'data':{'id':f'ST_{idx_1}',
                             'label':c,
                             'url':'./assets/cell.png'},
                    'classes':'CT',
                    'position':{'x':self.node_cols['Cell Types']['x_start'],'y':cell_start_y}}
                )
                cyto_elements.append(
                    {'data':{'source':'Main_Cell','target':f'ST_{idx_1}'}}
                )

                # Getting genes
                genes = matching_rows.filter(regex=self.node_cols['Genes']['abbrev']).dropna(axis=1)
                col_vals = genes.columns.values.tolist()
                col_vals = [i for i in col_vals if 'LABEL' in i]

                for idx,col in enumerate(col_vals):
                    cyto_elements.append(
                        {'data':{'id':col,
                                 'label':genes[col].tolist()[0],
                                 'url':'./assets/gene.png'},
                        'classes':'G',
                        'position':{'x':self.node_cols['Genes']['x_start'],'y':gene_start_y}}
                    )

                    cyto_elements.append(
                        {'data':{'source':col,'target':f'ST_{idx_1}'}}
                    )
                    gene_start_y+=75

        return cyto_elements

    def get_cyto_data(self,clicked):

        if not clicked is None:
            if 'ST' in clicked['id']:
                table_data = self.table_df.dropna(subset=['CT/1/ABBR'])
                table_data = table_data[table_data['CT/1/ABBR'].str.match(clicked['label'])]

                label = clicked['label']
                try:
                    id = table_data['CT/1/ID'].tolist()[0]
                    # Modifying base url to make this link to UBERON
                    base_url = self.node_cols['Cell Types']['base_url']
                    new_url = base_url+id.replace('CL:','')

                except IndexError:
                    print(table_data['CT/1/ID'].tolist())
                    id = ''
                
                try:
                    notes = table_data['CT/1/NOTES'].tolist()[0]
                except:
                    print(table_data['CT/1/NOTES'])
                    notes = ''

            elif 'Main_Cell' not in clicked['id']:
                
                table_data = self.table_df.dropna(subset=[clicked['id']])
                table_data = table_data[table_data[clicked['id']].str.match(clicked['label'])]

                base_label = '/'.join(clicked['id'].split('/')[0:-1])
                label = table_data[base_label+'/LABEL'].tolist()[0]

                id = table_data[base_label+'/ID'].tolist()[0]
                
                if self.node_cols['Anatomical Structure']['abbrev'] in clicked['id']:
                    base_url = self.node_cols['Anatomical Structure']['base_url']

                    new_url = base_url+id.replace('UBERON:','')
                else:
                    base_url = self.node_cols['Genes']['base_url']

                    new_url = base_url+id.replace('HGNC:','')

                try:
                    notes = table_data[base_label+'/NOTES'].tolist()[0]
                except KeyError:
                    notes = ''

            else:
                label = ''
                id = ''
                notes = ''
                new_url = ''
        else:
            label = ''
            id = ''
            notes = ''
            new_url = ''

        return f'Label: {label}', dcc.Link(f'ID: {id}', href = new_url), f'Notes: {notes}'
    
    def ingest_wsi(self,slide_name):

        new_key_name = self.slide_info_dict[slide_name]['key_name']
        old_key_name = self.wsi.slide_info_dict['key_name']
        new_url = self.wsi.image_url.replace(old_key_name,new_key_name)

        old_ftu_path = self.wsi.ftu_path
        old_spot_path = self.wsi.spot_path

        new_ftu_path = old_ftu_path.replace(self.wsi.slide_name.replace('.svs',''),slide_name.replace('.svs',''))
        new_spot_path = old_spot_path.replace(self.wsi.slide_name.replace('.svs',''),slide_name.replace('.svs',''))
        new_slide = WholeSlide(new_url,slide_name,self.slide_info_dict[slide_name],new_ftu_path,new_spot_path)

        self.wsi = new_slide

        self.update_hex_color_key()

        # Getting map_dict and spot_dict for overlays
        with open(new_ftu_path) as f:
            geojson_polys = geojson.load(f)

        with open(new_spot_path) as f:
            spot_geojson_polys = geojson.load(f)

        map_dict = {
            'url':new_slide.image_url,
            'FTUs':{
                'Glomeruli': {
                    'geojson':{'type':'FeatureCollection', 'features': [i for i in geojson_polys['features'] if i['properties']['structure']=='Glomeruli']},
                    'id': 'glom-bounds',
                    'color': '#390191',
                    'hover_color':'#666'
                },
                'Tubules': {
                    'geojson':{'type':'FeatureCollection', 'features': [i for i in geojson_polys['features'] if i['properties']['structure']=='Tubules']},
                    'id':'tub-bounds',
                    'color': '#e71d1d',
                    'hover_color': '#ff0b0a'
                },
                'Arterioles': {
                    'geojson':{'type':'FeatureCollection', 'features': [i for i in geojson_polys['features'] if i['properties']['structure']=='Arterioles']},
                    'id':'art-bounds',
                    'color': '#b6d7a8',
                    'hover_color': '#50f207'
                }
            }
        }

        spot_dict = {
            'geojson':spot_geojson_polys,
            'id': 'spot-bounds',
            'color': '#dffa00',
            'hover_color':'#9caf00'
        }

        new_children = [
            dl.Overlay(
                dl.LayerGroup(
                    dl.GeoJSON(data = map_dict['FTUs'][struct]['geojson'], id = map_dict['FTUs'][struct]['id'], options = dict(style = self.ftu_style_handle),
                                hideout = dict(color_key = self.hex_color_key, current_cell = self.current_cell, fillOpacity = self.cell_vis_val),
                                hoverStyle = arrow_function(dict(weight=5, color = map_dict['FTUs'][struct]['hover_color'], dashArray = '')))),
                    name = struct, checked = True, id = self.wsi.slide_info_dict['key_name']+'_'+struct)
            for struct in map_dict['FTUs']
            ]
        
        new_children += [
            dl.Overlay(
                dl.LayerGroup(
                    dl.GeoJSON(data = spot_dict['geojson'], id = spot_dict['id'], options = dict(style = self.ftu_style_handle),
                            hideout = dict(color_key = self.hex_color_key, current_cell = self.current_cell, fillOpacity = self.cell_vis_val),
                            hoverStyle = arrow_function(dict(weight=5, color = spot_dict['hover_color'], dashArray = '')))),
                    name = 'Spots', checked = False, id = self.wsi.slide_info_dict['key_name']+'_Spots')
            ]
    
        
        new_url = self.wsi.image_url

        # Getting new pie chart and state bar chart
        new_pie_chart, new_state_bar = self.update_roi_pie([],[],self.wsi.slide_bounds)

        center_point = [(self.wsi.slide_bounds[1]+self.wsi.slide_bounds[3])/2,(self.wsi.slide_bounds[0]+self.wsi.slide_bounds[2])/2]


        return new_url, new_children, center_point, new_pie_chart, new_state_bar

    def update_graph(self,ftu,plot,label):
        
        self.current_ftu = ftu
        # Filtering by selected FTU
        #current_data = self.metadata[self.metadata['ftu_type'].str.match(ftu)]
        current_data = []
        for f in self.metadata:
            if 'ftu_type' in f:
                if f['ftu_type'] == ftu:
                    current_data.append(f)

        if plot=='TSNE':
            #plot_data_x = current_data['x_tsne'].tolist()
            #plot_data_y = current_data['y_tsne'].tolist()
            plot_data_x = [i['x_tsne'] for i in current_data]
            plot_data_y = [i['y_tsne'] for i in current_data]

        elif plot=='UMAP':
            #plot_data_x = current_data['x_umap'].tolist()
            #plot_data_y = current_data['y_umap'].tolist()
            plot_data_x = [i['x_umap'] for i in current_data]
            plot_data_y = [i['y_umap'] for i in current_data]

        #custom_data = list(current_data.index)
        #label_data = current_data[label].tolist()

        custom_data = [i['ftu_name'] for i in current_data]
        # If the label is image_id or cluster
        try:
            label_data = [i[label] for i in current_data]
        except:
            # If the label is a main cell type or cell states of a main cell type
            try:
                label_data = [i['Main_Cell_Types'][label] for i in current_data]
            except:
                # Need to add something here for using cell states as a label
                label_data[i['Cell_States'] for i in current_data]


        graph_df = pd.DataFrame({'x':plot_data_x,'y':plot_data_y,'ID':custom_data,'Label':label_data})

        cluster_graph = go.Figure(px.scatter(graph_df,x='x',y='y',custom_data=['ID'],color='Label'))
        cluster_graph.update_layout(
            margin=dict(l=0,r=0,t=0,b=0)
        )

        return cluster_graph

    def grab_image(self,sample_info):

        slide_name = sample_info['image_id']
        if type(slide_name)==str:
            slide_name = [slide_name.replace('V10S15-103_','')]
        else:
            slide_name = [i.replace('V10S15-103_','') for i in slide_name.tolist()]

        #print(slide_name)
        img_list = []
        for idx,s in enumerate(slide_name):
            # openslide needs min_x, min_y, width, height
            try:
                min_x = int(sample_info['Min_x_coord'])
                min_y = int(sample_info['Min_y_coord'])
                width = int(sample_info['Max_x_coord'])-min_x
                height = int(sample_info['Max_y_coord'])-min_y
            except:
                min_x = int(sample_info['Min_x_coord'].tolist()[idx])
                min_y = int(sample_info['Min_y_coord'].tolist()[idx])
                width = int(sample_info['Max_x_coord'].tolist()[idx])-min_x
                height = int(sample_info['Max_y_coord'].tolist()[idx])-min_y      

            slide_path = [i for i in self.slide_paths if s in i]

            slide_path = slide_path[0]

            wsi = openslide.OpenSlide(slide_path)
            slide_region = wsi.read_region((min_x,min_y),0,(width,height))
            openslide.OpenSlide.close(wsi)

            img_list.append(resize(np.uint8(np.array(slide_region))[:,:,0:3],output_shape=(256,256,3)))

        return img_list        

    def update_selected(self,hover,selected):

        if 'cluster-graph.selectedData' in list(ctx.triggered_prop_ids.keys()):
            sample_ids = [i['customdata'][0] for i in selected['points']]
            #sample_info = self.metadata.loc[sample_ids]
            sample_info = [i for i in self.metadata if self.metadata['ftu_name'] in sample_ids]
        else:
            if hover is not None:
                sample_id = hover['points'][0]['customdata']
                #sample_info = self.metadata.loc[sample_id]
                sample_info = [i for i in self.metadata if self.metadata['ftu_name']==sample_id]
            else:
                #sample_info = self.metadata.iloc[0,:]
                sample_info = self.metadata[0]

        current_image = self.grab_image(sample_info)
        if len(current_image)==1:
            selected_image = go.Figure(px.imshow(current_image[0]))
        else:
            selected_image = go.Figure(px.imshow(np.stack(current_image,axis=0),animation_frame=0,binary_string=True,labels=dict(animation_frame=self.current_ftu)))
        
        selected_image.update_layout(
            margin=dict(l=0,r=0,t=0,b=0)
        )

        return selected_image
    


#if __name__ == '__main__':
def app(*args):

    run_type = 'local'

    slide_info_dict = {
        'XY01_IU-21-015F.svs':{
            'key_name':'XY01IU21015F',
            'bounds':[-121.4887696318595,0.0,-121.29151201587271,0.19456360965996605],
            'wsi_dims':[22012,21566]
        },
        'XY02_IU-21-016F.svs':{
            'key_name':'XY02IU21016F',
            'bounds':[-121.48876728121257,0.0,-121.29107290809065,0.18546976820936942],
            'wsi_dims':[22061,20558]
        },
        'XY03_IU-21-019F.svs':{
            'key_name':'XY03IU21019F',
            'bounds':[-121.48876962947176,0.0,-121.29142240206029,0.19455461087467554],
            'wsi_dims':[22022,21565]
        },
        'XY04_IU-21-020F.svs':{
            'key_name':'XY04IU21020F',
            'bounds':[-121.48876962708414,0.0,-121.29123421301975,0.19454563737274247],
            'wsi_dims':[22043,21564]
        }
    }

    try:
        run_type = os.environ['RUNTYPE']
    except:
        print(f'Using {run_type} run type')

    """
    print(f'Using {run_type} run type')
    print(f'Current working directory is: {os.getcwd()}')
    print(f'Contents of current working directory is: {os.listdir(os.getcwd())}')
    """
    if run_type == 'local':
        # For local testing
        base_dir = '/mnt/c/Users/Sam/Desktop/HIVE/'

        available_slides = glob(base_dir+'FFPE/*.svs')
        slide_names = [i.split('/')[-1] for i in available_slides]
        slide_name = slide_names[0]

        slide_url = 'http://localhost:5000/rgb/'+slide_info_dict[slide_name]["key_name"]+'/{z}/{x}/{y}.png?r=B04&r_range=[46,168]&g=B03&g_range=[46,168]&b=B02&b_range=[46,168]&'
        ftu_path = base_dir+'SpotNet_NonEssential_Files/CellAnnotations_GeoJSON/'+slide_name.replace('.svs','_scaled.geojson')
        spot_path = base_dir+'SpotNet_NonEssential_Files/CellAnnotations_GeoJSON/'+slide_name.replace('.svs','_Spots_scaled.geojson')
        cell_graphics_path = 'graphic_reference.json'
        asct_b_path = 'Kidney_v1.2 - Kidney_v1.2.csv'

        #metadata_path = './assets/cluster_metadata/'
        metadata_paths = [base_dir+'SpotNet_NonEssential_Files/CellAnnotations_GeoJSON/'+s.replace('.svs','_scaled.geojson') for s in list(slide_info_dict.keys())]

    elif run_type == 'web' or run_type=='AWS':
        # For test deployment
        if run_type == 'web':
            base_dir = os.getcwd()+'/mysite/'
            slide_info_path = base_dir+'/slide_info/'
        else:
            base_dir = os.getcwd()
            slide_info_path = base_dir+'assets/slide_info/'

        available_slides = glob(slide_info_path+'*.svs')
        slide_names = [i.split('/')[-1] for i in available_slides]
        slide_name = slide_names[0]

        slide_url = 'http://0.0.0.0:5000/rgb/'+slide_info_dict[slide_name]['key_name']+'/{z}/{x}/{y}.png?r=B04&r_range=[46,168]&g=B03&g_range=[46,168]&b=B02&b_range=[46,168]&'
        spot_path = slide_info_path+slide_name.replace('.svs','_Spots_scaled.geojson')
        ftu_path = slide_info_path+slide_name.replace('.svs','_scaled.geojson')
        cell_graphics_path = base_dir+'graphic_reference.json'
        asct_b_path = base_dir+'Kidney_v1.2 - Kidney_v1.2.csv'

        #metadata_path = base_dir+'assets/cluster_metadata/'
        metadata_paths = [slide_info_path+s.replace('.svs','_scaled.geojson') for s in list(slide_info_dict.keys())]
    
    # Adding slide paths to the slide_info_dict
    for slide,path in zip(slide_names,available_slides):
        slide_info_dict[slide]['slide_path'] = path

    # Reading dictionary containing paths for specific cell types
    cell_graphics_key = cell_graphics_path
    cell_graphics_json = json.load(open(cell_graphics_key))
    cell_names = []
    for ct in cell_graphics_json:
        cell_names.append(cell_graphics_json[ct]['full'])

    # Reading in clustering metadata (old way, separate json files for each structure)
    #glom_metadata = json.load(open(metadata_path+'FFPE_SpTx_Glomeruli.json'))
    #tub_metadata = json.load(open(metadata_path+'FFPE_SpTx_Tubules.json'))

    #metadata = pd.DataFrame.from_dict(glom_metadata,orient='index')
    #metadata = pd.concat([metadata,pd.DataFrame.from_dict(tub_metadata,orient='index')],axis=0,ignore_index=True)

    # Compiling info for clustering metadata
    metadata = []
    for m in metadata_paths:
        # Reading geojson file
        with open(m) as f:
            current_geojson = geojson.load(f)
        
        # Iterating through features and adding properties to metadata
        # This will not include any of the geometries data which would be the majority
        for f in current_geojson['features']:
            metadata.append(f['properties'])

    # Adding ASCT+B table to files
    asct_b_table = pd.read_csv(asct_b_path,skiprows=list(range(10)))

    #wsi = Slide(slide_path,spot_path,counts_path,ftu_path,ann_ids,counts_def_path)
    wsi = WholeSlide(slide_url,slide_name,slide_info_dict[slide_name],ftu_path,spot_path)

    external_stylesheets = [dbc.themes.LUX]

    # Calculating center point for initial layout
    current_slide_bounds = slide_info_dict[slide_name]['bounds']
    center_point = [(current_slide_bounds[1]+current_slide_bounds[3])/2,(current_slide_bounds[0]+current_slide_bounds[2])/2]
    
    # Getting map_dict and spot_dict for overlays
    with open(ftu_path) as f:
        geojson_polys = geojson.load(f)

    with open(spot_path) as f:
        spot_geojson_polys = geojson.load(f)



    map_dict = {
        'url':wsi.image_url,
        'FTUs':{
            'Glomeruli': {
                'geojson':{'type':'FeatureCollection', 'features': [i for i in geojson_polys['features'] if i['properties']['structure']=='Glomeruli']},
                'id': 'glom-bounds',
                'color': '#390191',
                'hover_color':'#666'
            },
            'Tubules': {
                'geojson':{'type':'FeatureCollection', 'features': [i for i in geojson_polys['features'] if i['properties']['structure']=='Tubules']},
                'id':'tub-bounds',
                'color': '#e71d1d',
                'hover_color': '#ff0b0a'
            },
            'Arterioles': {
                'geojson':{'type':'FeatureCollection', 'features': [i for i in geojson_polys['features'] if i['properties']['structure']=='Arterioles']},
                'id':'art-bounds',
                'color': '#b6d7a8',
                'hover_color': '#50f207'
            }
        }
    }

    spot_dict = {
        'geojson':spot_geojson_polys,
        'id': 'spot-bounds',
        'color': '#dffa00',
        'hover_color':'#9caf00'
    }

    main_layout = gen_layout(cell_names,slide_names,center_point,map_dict,spot_dict)

    main_app = DashProxy(__name__,external_stylesheets=external_stylesheets,transforms = [MultiplexerTransform()])
    vis_app = SlideHeatVis(main_app,main_layout,wsi,cell_graphics_key,asct_b_table,metadata, slide_info_dict,run_type)

    if run_type=='web':
        return vis_app.app

# Comment this portion out for web running
if __name__=='__main__':
    app()
