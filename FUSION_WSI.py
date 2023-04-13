"""

External file to hold the WholeSlide class used in FUSION


"""

import geojson
import json

import numpy as np

import shapely
from shapely.geometry import Polygon, Point, shape
from skimage.draw import polygon


class WholeSlide:
    def __init__(self,
                image_url,
                slide_name,
                slide_info_dict,
                ftu_path,
                spot_path,
                verbose = False):
        
        self.image_url = image_url
        self.ftu_path = ftu_path
        self.spot_path = spot_path
        self.slide_info_dict = slide_info_dict
        self.slide_bounds = self.slide_info_dict['bounds']
        self.slide_name = slide_name

        self.slide_ext = self.slide_name.split('.')[-1]

        if 'slide_path' in self.slide_info_dict:
            self.slide_path = self.slide_info_dict['slide_path']

        self.morphometrics_list = ['Area','Mesangial Area','Mesangial Fraction','Luminal Fraction','Arterial Area','Average TBM Thickness','Average Cell Thickness']

        # More print statements:
        """
        print(f'Slide Name in WholeSlide: {self.slide_name}')
        print(f'FTU Path in WholeSlide: {self.ftu_path}')
        print(f'Image URL in WholeSlide: {self.image_url}')
        """

        # Efficiency test (grouping together ftus and spots into boxes with 50 structures in each)
        self.group_n = 50

        # Processing ftus and spots
        if self.ftu_path is not None:
            self.ftus = self.process_ftus()
        else:
            self.ftus = {}
            
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
            ftu_polys[p]['Cluster'] = []

            # Adding morphometrics
            for m in self.morphometrics_list:
                ftu_polys[p][m] = []

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

                if 'Cluster' in i['properties']:
                    ftu_polys[p]['Cluster'].append(float(i['properties']['Cluster']))

                for m in self.morphometrics_list:
                    if m in i['properties']:
                        ftu_polys[p][m].append(i['properties'][m])
                    else:
                        ftu_polys[p][m].append(0)

                # Adding group info to the group list for a given ftu
                if group_count==self.group_n-1:

                    group_dict = {
                        'box':shapely.geometry.box(*group_bounds),
                        'polygons':ftu_polys[p]['polygons'],
                        'barcodes':ftu_polys[p]['barcodes'],
                        'main_counts':ftu_polys[p]['main_counts'],
                        'cell_states':ftu_polys[p]['cell_states'],
                        'Cluster':ftu_polys[p]['Cluster']
                    }

                    for m in self.morphometrics_list:
                        group_dict[m] = ftu_polys[p][m]

                    ftu_groups[p].append(group_dict)

                    # resetting back to empty/0
                    ftu_polys[p]['polygons'] = []
                    ftu_polys[p]['barcodes'] = []
                    ftu_polys[p]['main_counts'] = []
                    ftu_polys[p]['cell_states'] = []
                    ftu_polys[p]['Cluster'] = []

                    group_count = 0
                    group_bounds = []

            # Adding last group
            group_dict = {
                'box':shapely.geometry.box(*group_bounds),
                'polygons':ftu_polys[p]['polygons'],
                'barcodes':ftu_polys[p]['barcodes'],
                'main_counts':ftu_polys[p]['main_counts'],
                'cell_states':ftu_polys[p]['cell_states'],
                'Cluster':ftu_polys[p]['Cluster']
            }

            for m in self.morphometrics_list:
                group_dict[m] = ftu_polys[p][m]

            ftu_groups[p].append(group_dict)

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












