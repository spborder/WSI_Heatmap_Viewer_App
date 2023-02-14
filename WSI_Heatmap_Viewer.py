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
import geojson

import plotly.express as px
from plotly.subplots import make_subplots
import plotly.graph_objects as go
from matplotlib import cm

from dash import dcc, ctx, Dash
import dash_bootstrap_components as dbc
import dash_cytoscape as cyto
from dash_extensions.enrich import DashProxy, html, Input, Output, MultiplexerTransform, State

from timeit import default_timer as timer


def gen_layout(cell_types,slides_available,thumb):

    # Header
    header = dbc.Navbar(
        dbc.Container([
            dbc.Row([
                dbc.Col(html.Img(id='logo',src=('./assets/Lab_Logo_white.png'),height='75px'),md='auto'),
                dbc.Col([
                    html.Div([
                        html.H3('Whole Slide Image Cell Distribution'),
                        html.P('Spatial Transcriptomics')
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
                                    href='https://github.com/SarderLab',
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
        style={'margin-bottom':'20px'}
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
        ],style={'margin-bottom':'20px'}
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
        ],style={'margin-bottom':'20px'}
    )
    

    # View of WSI
    wsi_view = [
        dbc.Card(
            id='wsi-card',
            children=[
                dbc.CardHeader("Whole Slide Image Viewer"),
                dbc.CardBody([
                    html.Div(
                        id='cell-heat-loader',
                        children=[
                            dcc.Loading(
                                id="cell-loading",
                                type="cube",
                                children=[
                                    dcc.Graph(
                                        id="wsi-view",
                                        figure=go.Figure())

                                ]
                            ),
                            html.Div(
                                id = 'current-hover'
                            )
                        ]
                    ),
                    dbc.CardFooter(
                        dcc.RadioItems(
                            id= 'vis-types',
                            options=[
                                {'label':html.Div('Heatmap',style={'margin-right':'20px','padding-left':'15px'}),'value':'Heatmap'},
                                {'label':html.Div('Spots',style={'margin-right':'20px','padding-left':'15px'}),'value':'Spots'},
                                {'label':html.Div('Functional Tissue Units',style={'margin-right':'20px','padding-left':'15px'}),'value':'FTUs'}
                            ],
                            value = 'Heatmap',
                            inline = True
                        )
                    )
                ])
            ],style={'margin-bottom':'20px'}
        )
    ]

    # Contents of tabs
    thumbnail_select = dbc.Card(
        children = [
            dbc.CardHeader("Thumbnail View"),
            dbc.CardBody([
                dcc.Graph(
                    id="thumb-img",
                    figure=go.Figure(px.imshow(thumb))
                )]
            ),
            html.Div(),
            dbc.CardFooter([
                dbc.Row([
                    dbc.Col([
                        dbc.Label("ROI Size", html_for="roi-slider"),
                        dcc.Slider(
                            id="roi-slider",
                            min=5,
                            max=50,
                            step=5,
                            value=10
                        )
                    ],md=6),
                    dbc.Col([
                        dbc.Label("Thumbnail Transparency",html_for='thumb-slider'),
                        dcc.Slider(
                            id='thumb-slider',
                            min=0,
                            max=100,
                            step=10,
                            value=50
                        )
                    ],md=6)
                ])
            ])
    ],style={'margin-bottom':'20px'})

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
    

    # Tools for selecting regions, transparency, and cells
    tools = [
        dbc.Card(
            id='tools-card',
            children=[
                dbc.CardHeader("Tools"),
                dbc.CardBody([
                    html.H6("Select Cell for Heatmap View",className="cell-select"),
                    html.Div(
                        id = "cell-select-div",
                        children=[
                            dcc.Dropdown(cell_types,cell_types[0],id='cell-drop')
                        ]
                    ),
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
                            ],md=8),
                            dbc.Col([
                                dbc.Label(
                                    "Thumbnail Heatmap",
                                    html_for='thumb-heat'
                                ),
                                dcc.Checklist(
                                    id = 'thumb-heat',
                                    options = ["View Thumbnail Heatmap"],
                                    value = ["View Thumbnail Heatmap"]
                                )
                            ],md=4)
                        ]),
                        dbc.Row([
                            dbc.Tabs([
                                #dbc.Tab(thumbnail_select, label = "Thumbnail ROI"),
                                dbc.Tab(roi_pie, label = "Cell Composition"),
                                dbc.Tab(cell_card,label = "Cell Card")
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
                children = description
            ),
            html.B(),
            dbc.Row(
                id = 'slide-select-row',
                children = slide_select
            ),
            html.B(),
            dbc.Row(
                id="app-content",
                children=[dbc.Col(wsi_view,md=6),dbc.Col(thumbnail_select,md=6)]
            ),
            html.B(),
            dbc.Row(
                id='tab-info',
                children = [dbc.Col(tools,md=12)]
            )
        ],fluid=True)
    ])

    return main_layout


class Slide:
    def __init__(self,
                slide_path,
                spot_path,
                counts_path,
                ftu_path,
                ann_ids,
                counts_def_df):

        self.slide_path = slide_path
        self.spot_path = spot_path
        self.counts_path = counts_path
        self.counts_def_df = counts_def_df
        self.ftu_path = ftu_path
        self.ann_ids = ann_ids

        self.slide_name = self.slide_path.split(os.sep)[-1]
        slide_extension = self.slide_name.split('.')[-1]
        self.slide_name = self.slide_name.replace('.'+slide_extension,'')

        self.thumbnail_path = f'./assets/thumbnail_masks/{self.slide_name}/'
        if not os.path.exists(self.thumbnail_path):
            self.thumbnail_path = f'/home/samborder/mysite/assets/thumbnail_masks/{self.slide_name}/'


        self.current_cell = ''
        self.current_cell_thumb = Image.fromarray(np.uint8(np.zeros((256,256))))

        # Grouping together a set number of shapely objects and generating bounding boxes to reduce intersection calls
        self.group_n = 50

        if not self.counts_def_df is None:
            if type(self.counts_def_df)==str:
                self.counts_def_df = pd.read_csv(self.counts_def_df)

        self.counts = pd.read_csv(self.counts_path,index_col=0,engine='python')
        self.cell_types = list(self.counts.index)

        if not self.counts_def_df is None:
            self.counts_data = self.compress_counts()
            self.cell_types = self.counts_data['main_cell_types'].columns.tolist()

        self.slide = self.read_wsi()
        self.dimensions = self.slide.dimensions
        self.spots = self.read_spots()

        # Reading FTUs from ftu_path if provided
        if not self.ftu_path is None:
            self.ftus = self.read_ftus()

        self.thumb = self.get_thumbnail()
        print(f'Got thumbnail')
        if not os.path.exists(self.slide_path.replace('.svs','_thumb.png')):
            self.thumb.save(self.slide_path.replace('.svs','_thumb.png'))

    def compress_counts(self):
        
        sub_types_list = self.counts_def_df['Sub_Types'].tolist()
        cell_states_list = self.counts_def_df['Cell_States'].tolist()

        main_types_list = self.counts_def_df['Main_Types'].tolist()
        
        # Normalize to sum to 1
        norm_counts = self.counts/self.counts.sum(axis=0)

        slide_compressed_counts = {}
        slide_compressed_counts['main_cell_types'] = pd.DataFrame()
        
        for m in range(0,len(main_types_list)):
            main_name = main_types_list[m]
            sub_list = sub_types_list[m].split('.')
            state_list = cell_states_list[m].split('.')

            slide_compressed_counts[main_name] = {}

            sub_df = norm_counts.loc[norm_counts.index.isin(sub_list)]
            sub_sum_df = sub_df.sum(axis=0)
            if slide_compressed_counts['main_cell_types'].empty:
                slide_compressed_counts['main_cell_types'] = pd.DataFrame(data = sub_sum_df.values,columns = [main_name],index=list(sub_sum_df.index))
            else:
                new_df = pd.DataFrame(data = sub_sum_df.values, columns = [main_name],index=list(sub_sum_df.index))
                slide_compressed_counts['main_cell_types'] = pd.concat([slide_compressed_counts['main_cell_types'],new_df],axis=1)

            pct_count_df = (sub_df/sub_sum_df).fillna(0)
            slide_compressed_counts[main_name]['pct_subtypes'] = pct_count_df

            state_pct_df = pct_count_df.copy()

            # Setting the index of state percentages
            #if self.use_all_states:
            state_pct_df.index = [state_list[i] for i in range(len(state_list)) if sub_list[i] in list(state_pct_df.index)]

            # Combining states with the same name
            state_pct_df = state_pct_df.groupby(level=0).sum()

            # Re-normalizing (if certain cell states are removed) (replacing inf/nan with zero (for structures that have a sum of zero across included cell states))
            state_pct_df = (state_pct_df/state_pct_df.sum(axis=0)).fillna(0)

            # Sorting index so each one is consistent
            state_pct_df = state_pct_df.sort_index()

            slide_compressed_counts[main_name]['pct_states'] = state_pct_df

        # Adding in rows for missing main cell types
        if slide_compressed_counts['main_cell_types'].shape[1]<len(main_types_list):
            add_column_names = [i for i in main_types_list if i not in slide_compressed_counts['main_cell_types'].columns]
            for add in add_column_names:
                added_df = pd.DataFrame({add:[0]*slide_compressed_counts['main_cell_types'].shape[0]},index=list(slide_compressed_counts['main_cell_types'].index))
                slide_compressed_counts['main_cell_types'] = pd.concat([slide_compressed_counts['main_cell_types'],added_df],axis=1)

        slide_compressed_counts['main_cell_types'] = slide_compressed_counts['main_cell_types'].sort_index(axis=1)
    
        return slide_compressed_counts

    def read_wsi(self):
        return openslide.OpenSlide(self.slide_path)

    def get_thumbnail(self):

        if not os.path.exists(self.slide_path.replace('.svs','_thumb.png')):
            # Returns thumbnail PIL Image for a given slide with size (512,512)
            return self.slide.get_thumbnail((256,256))
        else:
            return Image.open(self.slide_path.replace('.svs','_thumb.png'))

    def read_xml(self,filepath):
        return ET.parse(filepath)

    def read_regions(self,region):
        Vertices = region.findall('./Vertices/Vertex')

        coords = []
        for Vertex in Vertices:
            coords.append((np.float32(Vertex.attrib['X']),np.float32(Vertex.attrib['Y'])))

        if len(coords)>2:
            reg_poly = Polygon(coords)
            barcode = region.attrib['Text']

            return reg_poly,barcode
        else:
            return None, None

    def read_spots(self):

        spot_polys = {}
        spot_polys['polygons'] = []
        spot_polys['barcodes'] = []

        spot_groups = []
        group_count = 0
        group_box = []

        spots = self.read_xml(self.spot_path).getroot().findall('Annotation[@Id="1"]/Regions/Region')
        for spot in tqdm(spots,desc='Parsing spots'):

            poly,barcode = self.read_regions(spot)

            if not poly==None:
                #spot_polys['polygons'].append(poly)
                #spot_polys['barcodes'].append(barcode)
                group_count+=1

                # Get polygon bounds, compare to group current bounds
                poly_bounds = list(poly.bounds)
                if not group_box == []:
                    # minimums
                    new_mins = [np.minimum(group_box[i],poly_bounds[i]) for i in range(0,2)]
                    new_maxs = [np.maximum(group_box[i],poly_bounds[i]) for i in range(2,4)]

                    group_box = new_mins+new_maxs
                else:
                    group_box = poly_bounds

                spot_polys['polygons'].append(poly)
                spot_polys['barcodes'].append(barcode)

                # Adding collected data to spot_groups list
                if group_count==self.group_n-1:
                    spot_groups.append({'box':shapely.geometry.box(*group_box),'polygons':spot_polys['polygons'],'barcodes':spot_polys['barcodes']})
                    spot_polys['polygons'] = []
                    spot_polys['barcodes'] = []
                    group_count = 0
                    group_box = []

        # Adding the remaining info to the spot_groups list (last one will have length<group count)
        spot_groups.append({'box':shapely.geometry.box(*group_box),'polygons':spot_polys['polygons'],'barcodes':spot_polys['barcodes']})

        return spot_groups

    def get_image(self,poly_box):
        
        # Returns PIL RGBA image from patch coordinates
        image = self.slide.read_region((int(poly_box[0]),int(poly_box[1])),0,(int(poly_box[2]-poly_box[0]),int(poly_box[3]-poly_box[1])))

        return image

    def find_intersecting_spots(self,box_poly):
        
        # Find intersecting bounding boxes for a given box poly (each box will contain up to self.group_n spots)
        intersecting_groups = [i for i in range(0,len(self.spots)) if self.spots[i]['box'].intersects(box_poly)]

        # Searching only in intersecting groups for spots that intersect
        intersect_spots = []
        intersect_barcodes = []
        for group_idx in intersecting_groups:

            intersect_idxes = [i for i in range(0,len(self.spots[group_idx]['polygons'])) if self.spots[group_idx]['polygons'][i].intersects(box_poly)]
            intersect_barcodes.extend([self.spots[group_idx]['barcodes'][i] for i in intersect_idxes])
            intersect_spots.extend([self.spots[group_idx]['polygons'][i] for i in intersect_idxes])


        #intersect_idxes = [i for i in range(0,len(self.spots['polygons'])) if self.spots['polygons'][i].intersects(box_poly)]

        #intersect_barcodes = [self.spots['barcodes'][i] for i in intersect_idxes]
        #intersect_spots = [self.spots['polygons'][i] for i in intersect_idxes]

        return intersect_spots, intersect_barcodes

    def read_ftus(self):

        if '.xml' in self.ftu_path:

            ftu_polys = {}
            ftu_tree = self.read_xml(self.ftu_path).getroot()
            
            for ann,id in self.ann_ids.items():
                
                ftu_polys[ann] = {}
                ftu_polys[ann]['polygons'] = []
                ftu_polys[ann]['barcodes'] = []

                ftus = ftu_tree.findall('Annotation[@Id="'+str(id)+'"]/Regions/Region')
                for idx,struct in tqdm(enumerate(ftus),desc = f'Parsing FTUs: {ann}'):

                    poly, barcode = self.read_regions(struct)

                    if not poly==None:
                        
                        ftu_polys[ann]['polygons'].append(poly)
                        ftu_polys[ann]['barcodes'].append(f'{ann}_{idx}')

        elif '.geojson' in self.ftu_path:

            # Reading files in geojson format
            with open(self.ftu_path) as f:
                geojson_polys = geojson.load(f)

            # Parsing through info stored in geojson file
            ftu_polys = {}
            poly_ids = []
            for f in geojson_polys:
                poly_ids.append(f['properties']['label'])

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

                #ftu_polys[p]['polygons'] = [shape(i['geometry']) for i in p_features]
                #ftu_polys[p]['barcodes'] = [i['properties']['label'] for i in p_features]
                #ftu_polys[p]['main_counts'] = [i['properties']['Main_Cell_Types'] for i in p_features]
                #ftu_polys[p]['cell_states'] = [i['properties']['Cell_States'] for i in p_features]

        return ftu_groups

    def find_intersecting_ftu(self,box_poly):

        intersect_barcodes = []
        intersect_ftus = []
        intersect_counts = []
        intersect_states = []
        for ann in list(self.ann_ids.keys()):

            # Which groups of FTUs intersect with the query box_poly
            group_intersect = [i for i in range(0,len(self.ftus[ann])) if self.ftus[ann][i]['box'].intersects(box_poly)]

            for g in group_intersect:
                group_polygons = self.ftus[ann][g]['polygons']
                group_barcodes = self.ftus[ann][g]['barcodes']
                group_counts = self.ftus[ann][g]['main_counts']
                group_states = self.ftus[ann][g]['cell_states']

                # Within group intersections
                intersect_idxes = [i for i in range(0,len(group_polygons)) if group_polygons[i].intersects(box_poly)]

                intersect_barcodes.extend([group_barcodes[i] for i in intersect_idxes])
                intersect_ftus.extend([group_polygons[i] for i in intersect_idxes])
                intersect_counts.extend([group_counts[i] for i in intersect_idxes])
                intersect_states.extend([group_states[i] for i in intersect_idxes])
                    

            #intersect_idxes= [i for i in range(0,len(self.ftus[ann]['polygons'])) if self.ftus[ann]['polygons'][i].intersects(box_poly)]
            #intersect_barcodes.extend([self.ftus[ann]['barcodes'][i] for i in intersect_idxes])
            #intersect_ftus.extend([self.ftus[ann]['polygons'][i] for i in intersect_idxes])
            #intersect_counts.extend([self.ftus[ann]['main_counts'][i] for i in intersect_idxes])
            #intersect_states.extend([self.ftus[ann]['cell_states'][i] for i in intersect_idxes])


        return {'polys':intersect_ftus, 'barcodes':intersect_barcodes, 'main_counts':intersect_counts, 'states':intersect_states}

    def thumbnail_overlay(self,cell_val,vis_val):
        
        if not cell_val == self.current_cell or self.current_cell == '':
            self.current_cell = cell_val

            # reading cell thumbnail
            cell_thumb = Image.open(self.thumbnail_path+f'{cell_val.replace("/","")}_thumbnail_vis.png')
            
            vis = np.array(cell_thumb)[:,:,0:3]
            zero_mask = np.where(np.sum(vis,axis=2)==0,0,vis_val)
            vis_mask_4d = np.concatenate((vis,zero_mask[:,:,None]),axis=-1)
            vis_mask_4d = Image.fromarray(np.uint8(vis_mask_4d)).convert('RGBA')

            self.current_cell_thumb = vis_mask_4d
            vis_overlay = self.thumb.copy()
            vis_overlay.paste(self.current_cell_thumb,mask=self.current_cell_thumb)
        
        elif cell_val == self.current_cell:

            vis = np.array(self.current_cell_thumb)[:,:,0:3]
            zero_mask = np.where(np.sum(vis,axis=2)==0,0,vis_val)
            vis_mask_4d = np.concatenate((vis,zero_mask[:,:,None]),axis=-1)
            vis_mask_4d = Image.fromarray(np.uint8(vis_mask_4d)).convert('RGBA')

            self.current_cell_thumb = vis_mask_4d
            vis_overlay = self.thumb.copy()
            vis_overlay.paste(self.current_cell_thumb,mask=self.current_cell_thumb)


        return vis_overlay



class SlideHeatVis:
    def __init__(self,
                app,
                layout,
                wsi,
                cell_graphics_key,
                asct_b_table,
                run_type = None):
                
        self.app = app
        self.app.title = "WSI Heatmap Viewer"
        self.app.layout = layout
        self.app._favicon = './assets/favicon.ico'

        self.run_type = run_type

        self.wsi = wsi
        # size here is in the form [width,height]
        self.wsi_size = self.wsi.dimensions
        print(f'wsi_size: {self.wsi_size}')

        self.cell_graphics_key = json.load(open(cell_graphics_key))
        # Inverting the graphics key to get {'full_name':'abbreviation'}
        self.cell_names_key = {}
        for ct in self.cell_graphics_key:
            self.cell_names_key[self.cell_graphics_key[ct]['full']] = ct


        # ASCT+B table for cell hierarchy generation
        self.table_df = asct_b_table    

        self.original_thumb = self.wsi.thumb.copy()
        self.roi_size = 2

        # FTU settings
        self.ftus = list(self.wsi.ann_ids.keys())
        self.ftu_boundary_thickness = 10
        # Random color generation to start (maybe add the ability to customize these later)
        self.ftu_color = {
            ann_name:[np.random.randint(0,255), np.random.randint(0,255), np.random.randint(0,255)]
            for ann_name in self.ftus
        }
        
        # Initializing some parameters
        self.patch_size = self.roi_size*30
        self.batch_size = 16
        self.current_cell = 'PT'
        self.current_barcodes = []

        # Cell Hierarchy related properties
        self.node_cols = {
            'Anatomical Structure':{'abbrev':'AS','x_start':25,'y_start':75},
            'Cell Types':{'abbrev':'CT','x_start':225,'y_start':0},
            'Genes':{'abbrev':'BGene','x_start':425,'y_start':75}
        }

        # Colormap settings (customize later)
        self.color_map = cm.get_cmap('jet')

        # Initial region overlay on thumbnail image
        initial_coords = [0,0]
        self.square_roi = self.gen_square(initial_coords)

        self.app.callback(
            [Output('thumb-img','figure'),Output('wsi-view','figure'),
            Output('roi-pie','figure'),Output('state-bar','figure'),
            Output('cell-graphic','src'),Output('cell-hierarchy','elements')],
            [Input('thumb-img','clickData'), Input('wsi-view','clickData'),
            Input('cell-drop','value'),Input('vis-slider','value'),
            Input('roi-slider','value'),Input('vis-types','value'),
            Input('thumb-heat','value'),Input('slide-select','value')]
        )(self.put_square)

        self.app.callback(
            Output('state-bar','figure'),
            Input('roi-pie','clickData')
        )(self.update_state_bar)

        self.app.callback(
            Output('current-hover','children'),
            Input('wsi-view','hoverData')
        )(self.get_hover)

        self.app.callback(
            [Output('collapse-content','is_open'),
            Output('collapse-descrip','children')],
            [Input('collapse-descrip','n_clicks'),
            Input('collapse-descrip','children')],
            [State('collapse-content','is_open')]
        )(self.view_instructions)

        self.app.callback(
            [Output('label-p','children'),
            Output('id-p','children'),
            Output('notes-p','children')],
            Input('cell-hierarchy','tapNodeData')
        )(self.get_cyto_data)
       
        # Comment out this line when running on the web
        if self.run_type == 'local':
            self.app.run_server(debug=True,use_reloader=True,port=8000)

        elif self.run_type == 'AWS':
            self.app.run_server(debug=True,use_reloader=False,port=8000)

    def view_instructions(self,n,text,is_open):
        if text == 'View Instructions':
            new_text = 'Hide Instructions'
        else:
            new_text = 'View Instructions'
        if n:
            return [not is_open,new_text]
        return [is_open,new_text]
    
    def gen_square(self,center_coords):

        # Given a center coordinate (x,y), generate a square of set size for pasting on thumbnail image
        square_mask = np.zeros((256,256))
        square_center = [i if i-(self.roi_size/2)>=0 else int(self.roi_size/2) for i in center_coords]

        square_poly = Polygon([
            (square_center[0]-(self.roi_size/2),square_center[1]+(self.roi_size/2)),
            (square_center[0]-(self.roi_size/2),square_center[1]-(self.roi_size/2)),
            (square_center[0]+(self.roi_size/2),square_center[1]-(self.roi_size/2)),
            (square_center[0]+(self.roi_size/2),square_center[1]+(self.roi_size/2)),
            (square_center[0]-(self.roi_size/2),square_center[1]+(self.roi_size/2))
            ])

        x_coords = [int(i[0]) for i in list(square_poly.exterior.coords)]
        y_coords = [int(i[1]) for i in list(square_poly.exterior.coords)]
        cc,rr = polygon(y_coords,x_coords,(square_mask.shape[1],square_mask.shape[0]))
        square_mask[cc,rr] = 1
        square_mask_3d = np.stack((square_mask,square_mask,square_mask),axis=2)
        square_mask_3d[:,:,0]*=255
        square_mask_3d[:,:,1]*=0
        square_mask_3d[:,:,2]*=0

        return square_mask_3d, square_poly
 
    def update_square_location(self, click_coords,thumb_view,cell_val,vis_val):
        
        new_square, square_center = self.gen_square(click_coords)

        zero_mask = np.where(np.sum(new_square.copy(),axis=2)==0,0,255)
        square_mask_4d = np.concatenate((new_square,zero_mask[:,:,None]),axis=-1)
        rgba_mask = Image.fromarray(np.uint8(square_mask_4d),'RGBA')

        if thumb_view:
            self.current_thumb = self.wsi.thumbnail_overlay(self.cell_names_key[cell_val],vis_val)
        else:
            self.current_thumb = self.original_thumb.copy()

        annotated_thumb = self.current_thumb.copy()
        annotated_thumb.paste(rgba_mask,mask=rgba_mask)

        return annotated_thumb, square_center

    def gen_patch_centers(self,image):
        # Creating patch coordinates on wsi view
        square_dim = int(self.batch_size**0.5)
        img_max_y, img_max_x = np.shape(image)[0], np.shape(image)[1]

        min_x, max_x = int(self.patch_size/2), img_max_x-int(self.patch_size/2)
        min_y, max_y = int(self.patch_size/2), img_max_y-int(self.patch_size/2)

        center_x_base = np.linspace(min_x,max_x,num=square_dim).astype(int).tolist()
        center_y_base = np.linspace(min_y,max_y,num=square_dim).astype(int).tolist()

        center_x_base = [i if i>=0 else 0 for i in center_x_base]
        center_y_base = [i if i>=0 else 0 for i in center_y_base]

        center_x = []
        center_y = []
        for i in range(square_dim):
            center_x.extend([center_x_base[i]]*square_dim)
            center_y.extend(center_y_base)
        
        batch_centers = list(zip(center_x,center_y))

        return batch_centers

    def make_cell_mask(self,poly_list,cell_vis_mask,wsi_vis_bounds,weight_list,shape_type,ftu_names = None):

        if shape_type == 'ftus':
            outline_mask_2D = np.zeros(np.shape(cell_vis_mask))
            outline_mask_3D = np.stack((outline_mask_2D,outline_mask_2D,outline_mask_2D),axis=-1)
        
        for p,w in tqdm(zip(poly_list,weight_list),total = len(poly_list)):
            x_coords = [int(i[0]-wsi_vis_bounds[0]) for i in list(p.exterior.coords)]
            y_coords = [int(i[1]-wsi_vis_bounds[1]) for i in list(p.exterior.coords)]
            
            intermed_mask = np.empty(np.shape(cell_vis_mask))
            intermed_mask[:] = np.nan

            if shape_type == 'patch':
                min_x = min(x_coords)
                max_x = max(x_coords)
                min_y = min(y_coords)
                max_y = max(y_coords)

                intermed_mask[min_y:max_y,min_x:max_x] = w

            elif shape_type == 'spots':
                cc,rr = polygon(y_coords,x_coords,(np.shape(cell_vis_mask)[0],np.shape(cell_vis_mask)[1]))
                intermed_mask[cc,rr] = w

            elif shape_type == 'ftus':
                # making shape outline
                outline = p.buffer(self.ftu_boundary_thickness)
                outline_x = [int(i[0]-wsi_vis_bounds[0]) for i in list(outline.exterior.coords)]
                outline_y = [int(i[1]-wsi_vis_bounds[1]) for i in list(outline.exterior.coords)]

                cc_o,rr_o = polygon(outline_y,outline_x,(np.shape(cell_vis_mask)[0],np.shape(cell_vis_mask)[1]))

                # making innerpolygon
                cc,rr = polygon(y_coords,x_coords,(np.shape(cell_vis_mask)[0],np.shape(cell_vis_mask)[1]))
                intermed_mask[cc,rr] = w

                outline_mask_3D[cc_o,rr_o,:] = self.ftu_color[ftu_names[poly_list.index(p)].split('_')[0]]
                outline_mask_3D[cc,rr,:] = [0,0,0]

            intermed_mask[intermed_mask==0] = np.nan
            cell_vis_mask = np.nanmean(np.stack((cell_vis_mask,intermed_mask),axis=-1),axis=-1)

        if shape_type == 'ftus':
            return cell_vis_mask, outline_mask_3D
        else:
            return cell_vis_mask

    def gen_cell_vis(self,cell_val,vis_val):
        
        # Getting intersecting spot polygons and barcodes
        intersect_spots, intersect_barcodes = self.wsi.find_intersecting_spots(self.current_projected_poly)

        self.current_spots = intersect_spots
        self.current_barcodes = intersect_barcodes
        intersect_counts = self.wsi.counts_data['main_cell_types'][self.wsi.counts_data['main_cell_types'].index.isin(intersect_barcodes)]

        self.current_counts = intersect_counts
        self.current_intersect_weight = [(spt.intersection(self.current_projected_poly).area)/(spt.area) for spt in self.current_spots]

        # Getting current cell value
        if not intersect_counts.empty:
            #intersecting_cell_vals = intersect_counts.loc[cell_val]
            intersecting_cell_vals = intersect_counts[cell_val]

        # Getting dimensions of proj_poly
        proj_bounds = list(self.current_projected_poly.bounds)
        cell_vis_overlay = np.empty(((int(proj_bounds[3]-proj_bounds[1]),int(proj_bounds[2]-proj_bounds[0]))))
        cell_vis_overlay[:] = np.nan

        if len(intersect_spots)>0:
            # Making a new patch for each point in self.patch_centers

            if self.current_vis_type == 'Heatmap':
                patch_list = []
                weight_list = []
                for center in self.patch_centers:
                    patch_poly = shapely.geometry.box(
                        center[0]-(self.patch_size/2)+proj_bounds[0],
                        center[1]-(self.patch_size/2)+proj_bounds[1],
                        center[0]+(self.patch_size/2)+proj_bounds[0],
                        center[1]+(self.patch_size/2)+proj_bounds[1],
                        center[0]-(self.patch_size/2)+proj_bounds[0]
                        )
                    intersect_weight = []
                    for spt in intersect_spots:
                        intersect_weight.append((spt.intersection(patch_poly).area)/(patch_poly.area))

                    mult_cell_vals = (intersecting_cell_vals*intersect_weight).sum()

                    patch_list.append(patch_poly)
                    weight_list.append(mult_cell_vals)

                cell_mask = self.make_cell_mask(patch_list,cell_vis_overlay,proj_bounds,weight_list,'patch')
            
            elif self.current_vis_type == 'Spots':

                cell_mask = self.make_cell_mask(self.current_spots,cell_vis_overlay,proj_bounds,intersecting_cell_vals.tolist(),'spots')

            elif self.current_vis_type == 'FTUs':
                counts = [i[cell_val] for i in self.current_ftus['main_counts']]
                cell_mask, outlines = self.make_cell_mask(self.current_ftus['polys'],cell_vis_overlay,proj_bounds,counts,'ftus',self.current_ftus['barcodes'])
                
        else:
            cell_mask = np.zeros((np.shape(cell_vis_overlay)))
        
        self.current_cell_distribution = cell_mask

        cell_heatmap = cell_mask.copy()
        cell_heatmap[np.where(np.isnan(cell_heatmap))] = 0

        color_cell_heatmap = np.uint8(255*self.color_map(np.uint8(255*cell_heatmap))[:,:,0:3])
        zero_mask = np.where(cell_heatmap==0,0,vis_val)
        cell_mask_4d = np.concatenate((color_cell_heatmap,zero_mask[:,:,None]),axis=-1)
        heatmap = Image.fromarray(np.uint8(cell_mask_4d)).convert('RGBA')

        #heatmap = self.current_vis.copy()

        if self.current_vis_type == 'FTUs':

            # overlapping ftu boundaries
            zero_mask_outlines = np.where(outlines.sum(axis=-1)==0,0,255)
            outline_mask_4D = np.concatenate((outlines,zero_mask_outlines[:,:,None]),axis=-1)
            outline_mask_4D = Image.fromarray(np.uint8(outline_mask_4D)).convert('RGBA')

            self.current_ftu_outlines = outline_mask_4D

        return heatmap

    def gen_projected_poly(self):

        projected_coords = list(self.current_square_poly.bounds)
        proj_x = [int(projected_coords[0]*(self.wsi_size[0]/256)),int(projected_coords[2]*(self.wsi_size[0]/256))]
        proj_y = [int(projected_coords[1]*(self.wsi_size[1]/256)),int(projected_coords[3]*(self.wsi_size[1]/256))]

        projected_poly = Polygon([
            (proj_x[0],proj_y[1]),
            (proj_x[0],proj_y[0]),
            (proj_x[1],proj_y[0]),
            (proj_x[1],proj_y[1]),
            (proj_x[0],proj_y[1])
        ])

        return projected_poly

    def gen_roi_pie(self):

        # Getting current counts and using those to generate a relative pie-chart
        if isinstance(self.current_counts,pd.DataFrame):
            combined_counts = self.current_counts.sum(axis=0).to_frame()
            combined_counts.columns = ['Current ROI']
            combined_counts = combined_counts.reset_index()
        else:
            combined_counts = self.current_counts.to_frame()
            combined_counts.columns = ['Current ROI']
            combined_counts = combined_counts.reset_index()
        
        combined_counts['Current ROI'] = combined_counts['Current ROI']/combined_counts['Current ROI'].sum()
        
        # Currently just doing top 5
        #combined_counts = combined_counts[combined_counts['Current ROI']>0]
        combined_counts = combined_counts.sort_values(by=['Current ROI'],ascending=False).iloc[0:5,:]

        cell_types_pie = px.pie(combined_counts,values='Current ROI',names = 'index')

        # By default just getting the top cell type
        top_cell = combined_counts['index'].tolist()[0]

        # Getting pct_states
        pct_states = self.wsi.counts_data[top_cell]['pct_states']
        # Getting spots from columns
        intersect_states = pct_states[pct_states.columns.intersection(self.current_barcodes)]
        aggregated_states = intersect_states.sum(axis=1).to_frame()
        aggregated_states = aggregated_states.reset_index()
        aggregated_states.columns = ['Cell State','Proportion']
        aggregated_states['Proportion'] = aggregated_states['Proportion']/aggregated_states['Proportion'].sum()

        state_bar = go.Figure(px.bar(aggregated_states,x='Cell State', y = 'Proportion',title=f'Cell State Proportions for {self.cell_graphics_key[top_cell]["full"]}'))

        return cell_types_pie, state_bar

    def put_square(self,thumb_click_point,wsi_click_point,cell_val,vis_val,roi_size,vis_type,thumb_view,slide_val):
        
        vis_val = int((vis_val/100)*255)
        self.patch_size = roi_size*30
        self.current_vis_type = vis_type
        self.current_cell = self.cell_names_key[cell_val]

        if ctx.triggered_id in ['thumb-img','wsi-view','slide-select',None]:
            
            if ctx.triggered_id=='slide-select':
                self.wsi=self.ingest_wsi(slide_val)

            if ctx.triggered_id == 'thumb-img':
                click_point = [thumb_click_point['points'][0]['x'],thumb_click_point['points'][0]['y']]
            elif ctx.triggered_id == 'wsi-view':
                # Getting wsi-view click coordinates
                wsi_click_coords = [wsi_click_point['points'][0]['x'],wsi_click_point['points'][0]['y']]
                # Getting wsi-view origin points in thumbnail-view space
                current_origin = [i-(self.roi_size/2) for i in self.current_click]
                current_width,current_height = self.current_wsi_image.size
                wsi_view_dims = [current_width,current_height]
                click_point = [k+(i*(self.roi_size/j)) for i,j,k in zip(wsi_click_coords,wsi_view_dims,current_origin)]

            elif ctx.triggered_id is None or ctx.triggered_id=='slide-select':
                click_point = [int(self.roi_size/2),int(self.roi_size/2)]

            new_square_image, square_poly = self.update_square_location(click_point,thumb_view,cell_val,vis_val)

            self.current_click = click_point

            self.current_square = new_square_image
            self.current_square_poly = square_poly

            self.current_projected_poly = self.gen_projected_poly()

            self.current_ftus = self.wsi.find_intersecting_ftu(self.current_projected_poly)

            # Timing current image extraction
            start = timer()
            self.current_wsi_image = self.wsi.get_image(list(self.current_projected_poly.bounds))
            end = timer()

            print(f'time to extract slide region: {end-start} seconds')

            self.patch_centers = self.gen_patch_centers(self.current_wsi_image)
            # Getting cell overlay
            # Timing for cell visualization
            start = timer()
            self.current_overlay = self.gen_cell_vis(self.current_cell,vis_val)
            end = timer()

            print(f'time to generate cell visualization with {len(self.current_ftus["polys"])} FTUs and {len(self.current_spots)} spots with ROI size: {self.current_projected_poly.area}: {end-start} (seconds)')
            vis_overlay = self.current_wsi_image.copy()
            vis_overlay.paste(self.current_overlay,mask=self.current_overlay)

            if self.current_vis_type=='FTUs':
                vis_overlay.paste(self.current_ftu_outlines,mask=self.current_ftu_outlines)
            
            self.current_wsi_view = vis_overlay

            wsi_fig = go.Figure(px.imshow(vis_overlay))
            square_fig = go.Figure(px.imshow(new_square_image))

        if ctx.triggered_id == 'cell-drop' or ctx.triggered_id == 'vis-types':
            
            start = timer()
            new_overlay = self.gen_cell_vis(self.current_cell,vis_val)
            end = timer()

            print(f'time to generate cell visualization with {len(self.current_ftus["polys"])} FTUs and {len(self.current_spots)} spots with ROI size: {self.current_projected_poly.area}: {end-start} (seconds)')
            
            print(thumb_view)
            if thumb_view == ['View Thumbnail Heatmap']:
                self.current_square, square_poly = self.update_square_location(self.current_click,thumb_view,cell_val,vis_val)

            self.current_overlay = new_overlay
            vis_overlay = self.current_wsi_image.copy()
            vis_overlay.paste(new_overlay,mask = new_overlay)

            if self.current_vis_type=='FTUs':
                vis_overlay.paste(self.current_ftu_outlines,mask=self.current_ftu_outlines)
            
            self.current_wsi_view = vis_overlay

            wsi_fig = go.Figure(px.imshow(self.current_wsi_view))
            square_fig = go.Figure(px.imshow(self.current_square))

        if ctx.triggered_id == 'vis-slider':

            vis = np.array(self.current_overlay)[:,:,0:3]
            zero_mask = np.where(np.sum(vis,axis=2)==0,0,vis_val)
            vis_mask_4d = np.concatenate((vis,zero_mask[:,:,None]),axis=-1)
            vis_mask_4d = Image.fromarray(np.uint8(vis_mask_4d)).convert('RGBA')

            self.current_overlay = vis_mask_4d
            vis_overlay = self.current_wsi_image.copy()
            vis_overlay.paste(self.current_overlay,mask=self.current_overlay)

            if thumb_view == 'View Thumbnail Heatmap':
                self.current_square, square_poly = self.update_square_location(self.current_click,thumb_view,cell_val,vis_val)

            if self.current_vis_type=='FTUs':
                vis_overlay.paste(self.current_ftu_outlines,mask=self.current_ftu_outlines)

            self.current_wsi_view = vis_overlay

            wsi_fig = go.Figure(px.imshow(self.current_wsi_view))
            square_fig = go.Figure(px.imshow(self.current_square))

        if ctx.triggered_id == 'roi-slider':
            # Resizing ROI
            self.roi_size = roi_size
            new_square_image, square_poly = self.update_square_location(self.current_click,thumb_view,cell_val,vis_val)

            self.current_square = new_square_image
            self.current_square_poly = square_poly

            self.current_projected_poly = self.gen_projected_poly()
            self.current_ftus = self.wsi.find_intersecting_ftu(self.current_projected_poly)

            self.current_wsi_image = self.wsi.get_image(list(self.current_projected_poly.bounds))

            self.patch_centers = self.gen_patch_centers(self.current_wsi_image)

            # Getting cell overlay
            self.current_overlay = self.gen_cell_vis(self.current_cell,vis_val)
            vis_overlay = self.current_wsi_image.copy()
            vis_overlay.paste(self.current_overlay,mask=self.current_overlay)

            if self.current_vis_type=='FTUs':
                vis_overlay.paste(self.current_ftu_outlines,mask=self.current_ftu_outlines)

            self.current_wsi_view = vis_overlay
            wsi_fig = go.Figure(px.imshow(self.current_wsi_view))
            square_fig = go.Figure(px.imshow(self.current_square))

        if ctx.triggered_id == 'thumb-heat':
            
            if thumb_view=='View Thumbnail Heatmap':

                self.current_square, square_poly = self.update_square_location(self.current_click,thumb_view,cell_val,vis_val)

            else:

                self.current_square, square_poly = self.update_square_location(self.current_click,thumb_view,cell_val,vis_val)

            wsi_fig = go.Figure(px.imshow(self.current_wsi_view))
            square_fig = go.Figure(px.imshow(self.current_square))

        wsi_fig.update_layout(
            margin=dict(l=10,r=10,t=10,b=10),
        )

        square_fig.update_layout(
            margin=dict(l=10,r=10,t=10,b=10)
        )
        
        cell_types_pie, state_barchart = self.gen_roi_pie()

        # Loading the cell-graphic and hierarchy image
        cell_graphic = self.cell_graphics_key[self.current_cell]['graphic']
        #cell_hierarchy = self.cell_graphics_key[self.current_cell]['hierarchy']
        cell_hierarchy = self.gen_cyto()

        if cell_graphic == "":
            cell_graphic = './assets/cell_graphics/default_cell_graphic.png'
        if cell_hierarchy == "":
            cell_hierarchy = './assets/cell_graphics/default_cell_hierarchy.png'

        return square_fig, wsi_fig, cell_types_pie, state_barchart, cell_graphic, cell_hierarchy

    def update_state_bar(self,cell_click):
        
        if not cell_click is None:
            self.pie_cell = cell_click['points'][0]['label']
        else:
            self.pie_cell = self.current_cell

        # Getting pct_states
        pct_states = self.wsi.counts_data[self.pie_cell]['pct_states']
        # Getting spots from columns
        intersect_states = pct_states[pct_states.columns.intersection(self.current_barcodes)]
        aggregated_states = intersect_states.sum(axis=1).to_frame()
        aggregated_states = aggregated_states.reset_index()
        aggregated_states.columns = ['Cell State','Proportion']
        aggregated_states['Proportion'] = aggregated_states['Proportion']/aggregated_states['Proportion'].sum()

        state_bar = go.Figure(px.bar(aggregated_states,x='Cell State', y = 'Proportion',title=f'Cell State Proportions for {self.cell_graphics_key[self.pie_cell]["full"]}'))

        return state_bar
    
    def get_hover(self,wsi_hover):

        hover_text = ''
        if not wsi_hover is None:
            hover_point = [wsi_hover['points'][0]['x'],wsi_hover['points'][0]['y']]
        else:
            hover_point = [0,0]

        # Getting hover data relative to wsi coordinates
        square_bounds = list(self.current_projected_poly.bounds)
        if self.current_vis_type in ['Spots','FTUs']:
            wsi_point = Point(hover_point[0]+square_bounds[0],hover_point[1]+square_bounds[1])

            if self.current_vis_type == 'Spots':
                spot_idx = np.argwhere([i.intersects(wsi_point) for i in self.current_spots])[0]
                barcode = self.current_barcodes[spot_idx[0]]
                counts = self.current_counts.loc[barcode][self.current_cell]
                hover_text += f'{self.cell_graphics_key[self.current_cell]["full"]}: {str(round(counts,3))}'
            
            if self.current_vis_type == 'FTUs':
                ftu_idx = np.argwhere([i.intersects(wsi_point) for i in self.current_ftus['polys']])[0]
                cell_val = self.current_ftus['main_counts'][ftu_idx[0]][self.current_cell]
                hover_text += f'{self.cell_graphics_key[self.current_cell]["full"]}: {str(round(cell_val,3))}'
        else:
            # For heatmap visualization
            hover_text += f'{self.cell_graphics_key[self.current_cell]["full"]}: {str(round(self.current_cell_distribution[hover_point[1],hover_point[0]]))}'

        return hover_text
        
    def gen_cyto(self):

        cyto_elements = []

        # Getting cell sub-types under that main cell
        cell_subtypes = self.cell_graphics_key[self.current_cell]['subtypes']

        # Getting all the rows that contain these sub-types
        table_data = self.table_df.dropna(subset = ['CT/1/ABBR'])
        cell_data = table_data[table_data['CT/1/ABBR'].isin(cell_subtypes)]

        # cell type
        cyto_elements.append(
            {'data':{'id':'Main_Cell','label':self.current_cell},'position':{'x':self.node_cols['Cell Types']['x_start'],'y':self.node_cols['Cell Types']['y_start']}}
        )

        # Getting the anatomical structures for this cell type
        an_structs = cell_data.filter(regex=self.node_cols['Anatomical Structure']['abbrev']).dropna(axis=1)

        an_start_y = self.node_cols['Anatomical Structure']['y_start']
        col_vals = an_structs.columns.values.tolist()
        col_vals = [i for i in col_vals if 'LABEL' in i]

        for idx,col in enumerate(col_vals):
            cyto_elements.append(
                {'data':{'id':col,'label':an_structs[col].tolist()[0]},'position':{'x':self.node_cols['Anatomical Structure']['x_start'],'y':an_start_y}}
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
                    {'data':{'id':f'ST_{idx_1}','label':c},'position':{'x':self.node_cols['Cell Types']['x_start'],'y':cell_start_y}}
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
                        {'data':{'id':col,'label':genes[col].tolist()[0]},'position':{'x':self.node_cols['Genes']['x_start'],'y':gene_start_y}}
                    )

                    cyto_elements.append(
                        {'data':{'source':col,'target':f'ST_{idx_1}'}}
                    )
                    gene_start_y+=75

        return cyto_elements

    def get_cyto_data(self,clicked):

        if 'ST' in clicked['id']:
            table_data = self.table_df.dropna(subset=['CT/1/ABBR'])
            table_data = table_data[table_data['CT/1/ABBR'].str.match(clicked['label'])]
            print(table_data)

            label = clicked['label']
            try:
                id = table_data['CT/1/ID'].tolist()[0]
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
            try:
                notes = table_data[base_label+'/NOTES'].tolist()[0]
            except KeyError:
                notes = ''

        else:
            label = ''
            id = ''
            notes = ''

        return f'Label: {label}', f'ID: {id}', f'Notes: {notes}'
    
    def ingest_wsi(self,slide_name):

        old_paths = [self.wsi.slide_path, self.wsi.spot_path,self.wsi.counts_path,
                    self.wsi.ftu_path,self.wsi.ann_ids,self.wsi.counts_def_df]
        
        wsi_ext = slide_name.split('.')[-1]
        new_paths = [i.replace(self.wsi.slide_name,slide_name.replace('.'+wsi_ext,'')) if type(i)==str else i for i in old_paths]

        new_slide = Slide(*new_paths)

        return new_slide


#if __name__ == '__main__':
def app(*args):

    run_type = 'local'

    try:
        run_type = os.environ['RUNTYPE']
    except:
        print(f'Using {run_type} run type')

    
    if run_type == 'local':
        # For local testing
        base_dir = '/mnt/c/Users/Sam/Desktop/HIVE/'

        available_slides = glob(base_dir+'FFPE/*.svs')
        slide_names = [i.split('/')[-1] for i in available_slides]
        slide_name = slide_names[0]

        slide_path = base_dir+'FFPE/'+slide_name
        spot_path = base_dir+'FFPE/Spot_Coordinates_large/'+slide_name.replace('.svs','_Large.xml')
        counts_path = base_dir+'counts_data/FFPE/CellTypeFractions_SpotLevel/V10S15-103_'+slide_name.replace('.svs','_cellfract.csv')
        counts_def_path = base_dir+'counts_data/Cell_SubTypes_Grouped.csv'
        ftu_path = base_dir+'SpotNet_NonEssential_Files/CellAnnotations_GeoJSON/'+slide_name.replace('.svs','.geojson')
        cell_graphics_path = 'graphic_reference.json'
        asct_b_path = 'Kidney_v1.2 - Kidney_v1.2.csv'

    elif run_type == 'web':
        # For test deployment
        base_dir = os.getcwd()+'/mysite/'

        available_slides = glob(base_dir+'/slide_info/*.svs')
        slide_names = [i.split('/')[-1] for i in available_slides]
        slide_name = slide_names[0]

        slide_path = base_dir+'/slide_info/'+slide_name
        spot_path = slide_path.replace('.svs','_Large.xml')
        counts_path = base_dir+'/slide_info/V10S15-103_'+slide_name.replace('.svs','_cellfract.csv')
        counts_def_path = slide_path.replace(slide_name,'Cell_SubTypes_Grouped.csv')
        ftu_path = slide_path.replace('.svs','.geojson')
        cell_graphics_path = base_dir+'graphic_reference.json'
        asct_b_path = base_dir+'Kidney_v1.2 - Kidney_v1.2.csv'

    ann_ids = None

    # Reading dictionary containing paths for specific cell types
    cell_graphics_key = cell_graphics_path
    cell_graphics_json = json.load(open(cell_graphics_key))
    cell_names = []
    for ct in cell_graphics_json:
        cell_names.append(cell_graphics_json[ct]['full'])


    # Adding ASCT+B table to files
    asct_b_table = pd.read_csv(asct_b_path,skiprows=list(range(10)))

    wsi = Slide(slide_path,spot_path,counts_path,ftu_path,ann_ids,counts_def_path)

    external_stylesheets = [dbc.themes.LUX]

    main_layout = gen_layout(cell_names,slide_names,wsi.thumb)

    main_app = DashProxy(__name__,external_stylesheets=external_stylesheets,transforms = [MultiplexerTransform()])
    vis_app = SlideHeatVis(main_app,main_layout,wsi,cell_graphics_key,asct_b_table,run_type)

    if run_type=='web':
        return vis_app.app

# Comment this portion out for web running
if __name__=='__main__':
    app()
