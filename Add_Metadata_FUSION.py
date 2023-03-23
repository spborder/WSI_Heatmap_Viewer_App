"""

Script for adding data to exiting GeoJSON/Histomics annotation files 

Input can be either:
    (1) GeoJSON formatted
    (2) JSON formatted
    (3) CSV/XLSX tables (with bounding box coordinates as columns data)
    (4) Maybe eventually RDS for some more sequencing-related analyses

"""


import os
import sys
import numpy as np

from glob import glob

import pandas as pd
from shapely.geometry import Polygon, shape
import shapely
import json
import geojson
from geojson import Feature, dump

import uuid

from tqdm import tqdm


class CurrentAnnotations:
    def __init__(self,
                current_ann_path,
                wsi_dims_data = None):

        self.current_ann_path = current_ann_path

        if 'geojson' in self.current_ann_path:
            self.ann_type = 'geojson'
        else:
            self.ann_type = 'histomics'

        self.wsi_dims_data = wsi_dims_data

        if not self.wsi_dims_data is None:
            self.make_output_scale()

        self.current_json_poly, self.current_ftus = self.load_current_annotations()

    def make_output_scale(self):
        # longitude = left-right
        # latitude = up-down
        # output_crs is in format [min longitude, min latitude, max longitude, max latitude]
        output_width = self.wsi_dims_data['bounds'][2]-self.wsi_dims_data['bounds'][0]        
        output_height = self.wsi_dims_data['bounds'][3]-self.wsi_dims_data['bounds'][1]

        # wsi dims is in format [width, height]
        self.width_scale = output_width/self.wsi_dims_data['wsi_dims'][0]
        self.height_scale = output_height/self.wsi_dims_data['wsi_dims'][1]

    def load_current_annotations(self):

        if self.ann_type == 'geojson':
            with open(self.current_ann_path) as f:
                json_polys = geojson.load(f)

            # Loading list of polygons
            poly_list = []
            for p in json_polys['features']:
                poly_list.append(shape(p['geometry']))
            
        else:
            with open(self.current_ann_path) as f:
                json_polys = json.load(f)
            
            print(f'Length of json_polys: {len(json_polys)}')
            # Loading list of polygons
            poly_list = []
            self.n_polys_in_structure = []
            for p in json_polys:
                if 'elements' in p:
                    boundary_coords = p['elements']
                    self.n_polys_in_structure.append(len(boundary_coords))
                    for b in boundary_coords:
                        # Removing "z" coordinate, which is all zeros in our case
                        shapely_poly_coords = [i[:-1] for i in b['points']]
                        poly_list.append(Polygon(shapely_poly_coords))

        return json_polys, poly_list

    def add_metadata(self,add_meta):
        """
        In this case, add_meta is a list of dictionaries with the polygon info and the metadata to add.

        add_meta: [
          {'polygon': {
              'type': 'polygon' or 'box',
              'coords': coordinates
          }
          'meta_labels':{
              'add-label-1': label(str) or value(num)
          }},
          etc.    
        ]
        """
        # Iterating through each member in add_meta list
        for m in tqdm(add_meta,desc='Adding Metadata'):

            # Creating shapely polygon for this object
            if m['polygon']['type'] == 'polygon':
                if self.wsi_dims_data is None:
                    m_poly = Polygon(m['polygon']['coords'])
                else:
                    # Normalizing coordinates in bounds of WSI (used for Terracotta deployment, geographic coordinate reference scheme (CRS) (epsg:32611))
                    norm_coords = []
                    for c in m['polygon']['coords']:
                        norm_coords.append(
                            (
                                self.wsi_dims_data['bounds'][0]+(c[0]*self.width_scale),
                                self.wsi_dims_data['bounds'][1]+(c[1]*self.height_scale)
                            )
                        )
                    m_poly = Polygon(norm_coords)

            elif m['polygon']['type'] == 'box':
                if self.wsi_dims_data is None:

                    edited_coords = []
                    for k in m['polygon']['coords']:
                        if len(k)==1:
                            edited_coords.append(k)
                        elif len(k)==2:
                            edited_coords.extend(k)

                    m_poly = shapely.geometry.box(*edited_coords)

                else:
                    # Normalizing coordinates in bounds of WSI
                    norm_coords = []
                    for c in m['polygon']['coords']:
                        norm_coords.extend(
                            [
                                self.wsi_dims_data['bounds'][0]+(c[0]*self.width_scale),
                                self.wsi_dims_data['bounds'][1]+(c[1]*self.height_scale)
                            ]
                        )
                    m_poly = shapely.geometry.box(*norm_coords)

            # Finding intersection in current_ftus
            intersect_list = [i for i in range(0,len(self.current_ftus)) if m_poly.intersects(self.current_ftus[i])]

            # If there is more than one intersection than only assign metadata to largest percentage area intersection?
            if len(intersect_list)>1:
                # Finding intersection between current metadata bounding box and bounding boxes of intersecting FTUs
                intersect_areas = [m_poly.intersection(shapely.geometry.box(*list(self.current_ftus[i].bounds))).area/m_poly.area for i in intersect_list]
                # Finding the ftu with the largest percentage overlap
                if any([i>=0.8 for i in intersect_areas]):
                    intersect_ftu = intersect_list[np.argmax(intersect_areas)]
                else:
                    intersect_ftu = None

            elif len(intersect_list)==1:
                if m_poly.intersection(shapely.geometry.box(*list(self.current_ftus[intersect_list[0]].bounds))).area/m_poly.area >= 0.8:
                    intersect_ftu = intersect_list[0]
                else:
                    intersect_ftu = None

            # If there is an overlapping ftu with this metadata object
            if not intersect_ftu is None:
                # Now, just add the meta_labels data to that index in self.current_json_poly
                if self.ann_type == 'geojson':
                    if not 'properties' in self.current_json_poly['features'][intersect_ftu]:
                        self.current_json_poly['features'][intersect_ftu]['properties'] = m['meta_labels']
                    else:
                        for l in m['meta_labels']:
                            if not l in self.current_json_poly['features'][intersect_ftu]['properties']:
                                self.current_json_poly['features'][intersect_ftu]['properties'][l] = m['meta_labels'][l]
                
                elif self.ann_type == 'histomics':
                    struct_idx, ftu_idx = self.find_aligning_index(intersect_ftu)

                    for l in m['meta_labels']:
                        if l not in self.current_json_poly[struct_idx]['elements'][ftu_idx]['user']:
                            self.current_json_poly[struct_idx]['elements'][ftu_idx]['user'][l] = m['meta_labels'][l]
            else:
                print('No intersecting FTUs')
                print(m['meta_labels'])

    def save_metadata(self,save_path):

        with open(save_path,'w') as f:
            dump(self.current_json_poly,f)

    def find_aligning_index(self,idx):
        # Special case for Histomics annotations which are stored as a list for each structure
        # This function just translates the intersect index in the total list (combined) into the index for the specific structure
        test_val = idx
        struct_idx = -1
        while test_val>=0 and struct_idx<len(self.n_polys_in_structure):
            idx = test_val
            struct_idx+=1
            test_val -= self.n_polys_in_structure[struct_idx]
        
        return struct_idx,idx


def main(meta_args):

    # Add a '*' to reference more than one metadata
    current_meta_path = meta_args['metadata_path']
    
    if '*' in current_meta_path:
        current_metas = glob(current_meta_path)
    else:
        current_metas = current_meta_path

    # slide_id used for separating data for each slide to add to the slide-specific annotations file
    add_meta_path = meta_args['add_meta_path']
    slide_id = meta_args['slide_id']
    extra_labels = meta_args['extra_labels']
    if extra_labels == []:
        extra_labels = None
    ignore_columns = meta_args['ignore_columns']
    poly_type = meta_args['poly_type']

    slide_info_dict = meta_args['slide_info_dict']

    if slide_info_dict == []:
        slide_info_dict = None


    for slide_id_idx,meta in enumerate(add_meta_path):

        slide_id_id = slide_id[slide_id_idx] 
        # Have to first parse the metadata into slide-dictionaries 
        if 'json' in meta:
            if not 'geojson' in meta:
                with open(meta) as f:
                    all_add_meta_dict = json.load(f)
            else:
                with open(meta) as f:
                    all_add_meta_dict = geojson.load(f)

            # Creating the add_meta_list
            if not slide_id_id is None:
                # TODO: make this slide_id_id thing also apply to geojson files, not sure how conversion to dataframes works with those because the metadata would be a bit more nested
                meta_frame = pd.DataFrame.from_dict(all_add_meta_dict,orient='index')

                for slide in meta_frame[slide_id_id].unique():
                    slide_metadata = meta_frame[meta_frame[slide_id_id].str.match(slide)]
                    add_meta_data_list = slide_metadata.to_dict('records')
                    print(f'Number of objects in: {slide}: {len(add_meta_data_list)}')
                    if poly_type == 'box':
                        add_metadata = [
                            {
                                'polygon':{
                                    'type':'box',
                                    'coords': [[i['Min_x_coord'], i['Min_y_coord']],[i['Max_x_coord'],i['Max_y_coord']]]
                                },
                                'meta_labels':i
                            }
                            for i in add_meta_data_list
                        ]
                    else:
                        print('need to add functionality here')

                    if not slide_info_dict is None:
                        current_annotations = CurrentAnnotations(current_ann_path = [i for i in current_metas if slide in i][0],wsi_dims_data=slide_info_dict[slide+'.svs'])
                    else:
                        current_annotations = CurrentAnnotations(current_ann_path = [i for i in current_metas if slide in i][0])

                    current_annotations.add_metadata(add_metadata)
                    current_annotations.save_metadata([i for i in current_metas if slide in i][0])


        elif 'csv' in meta:

            # This is for feature files from Nick
            if type(slide_id_id)==list:
                
                for idx,slide in enumerate(slide_id_id):
                    print(slide)
                    if not slide_info_dict is None:
                        current_annotations = CurrentAnnotations(current_ann_path = current_metas[idx],wsi_dims_data = slide_info_dict[slide+'.svs'])
                    else:
                        current_annotations = CurrentAnnotations(current_ann_path = current_metas[idx])

                    if extra_labels is None:
                        meta_df = pd.read_csv(meta.replace('*',slide))
                        add_meta_data_list = meta_df.to_dict('records')

                        # Do the same kinda thing as with JSON since there'll probably be columns for bbox coords
                        add_metadata = [
                            {
                                'polygon':{
                                    'type':'box',
                                    'coords': [[i['y1'],i['x1']],[i['y2'],i['x2']]]
                                },
                                'meta_labels': {j:k for j,k in i.items() if j not in ignore_columns}
                            }
                            for i in add_meta_data_list
                        ]

                        current_annotations.add_metadata(add_metadata)
                    else:
                        for e in extra_labels:
                            meta_df = pd.read_csv(meta.replace('*',slide+e))
                            add_meta_data_list = meta_df.to_dict('records')

                            # Do the same kinda thing as with JSON since there'll probably be columns for bbox coords
                            add_metadata = [
                                {
                                    'polygon':{
                                        'type':'box',
                                        'coords': [i['y1'],i['x1'],i['y2'],i['x2']]
                                    },
                                    'meta_labels': {j:k for j,k in i.items() if j not in ignore_columns}
                                }
                                for i in add_meta_data_list
                            ]

                            current_annotations.add_metadata(add_metadata)            

                    current_annotations.save_metadata(current_metas[idx])

if __name__=='__main__':

    # Reading JSON file with parameters:
    meta_args = sys.argv[-1]
    meta_args = json.load(open(meta_args))

    main(meta_args)