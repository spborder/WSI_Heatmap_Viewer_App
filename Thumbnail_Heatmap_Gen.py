"""

Whole slide image cell heatmap generation projected on thumbnail image

Procedure:
    - Read spots and counts
    - Read WSI
    - Read thumbnail
    - Get spot center locations
    - Project spot center points to thumbnail coordinates
    - Buffer spot points by spot_buffer
    - Get bounding box of all spots
    - Generate list of patch center coordinates equally distributed over tissue bounding box
    - Generate equal sized patches for each patch coordinate
    - Find intersecting spots with each patch
    - Combine counts info from intersecting spots (mean)
    - Assign mean value for each cell type to each patch
    - Combine all patches (ignore zero'd locations)
    - Apply colormap
    - Save heatmap thumbnail for each cell type separately

"""

import os
import sys
import numpy as np
import pandas as pd

import openslide

from PIL import Image

import lxml.etree as ET
from tqdm import tqdm

import shapely
from shapely.geometry import Polygon, Point, shape
from skimage.draw import polygon

from matplotlib import cm


class Thumbnail_Heatmap_Generator:
    def __init__(self,
                wsi_path,
                spot_path,
                counts_path,
                counts_def_path,
                save_path):

        self.wsi_path = wsi_path
        self.spot_path = spot_path
        self.counts_path = counts_path
        self.counts_def_path = counts_def_path
        self.save_path = save_path

        if not os.path.exists(self.save_path):
            os.makedirs(self.save_path)

        self.group_n = 50
        self.spot_patch_buffer = 10
        self.thumbnail_size = 256
        self.patch_n = 900
    
        self.patch_size = 20
        self.color_map = cm.get_cmap('jet')

        self.counts_def_df = pd.read_csv(self.counts_def_path)
        self.counts = pd.read_csv(self.counts_path,index_col=0,engine='python')
        
        print('Compressing counts data')
        self.counts_data = self.compress_counts()
        self.cell_names = self.counts_data['main_cell_types'].columns.values.tolist()

        self.wsi = self.read_wsi()

        self.wsi_width, self.wsi_height = self.wsi.dimensions


        print('Getting thumbnail')
        self.thumbnail = self.get_thumbnail()
        print('Reading spots')
        self.spots = self.read_spots()


        self.generate_thumbnail_map()
        

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
        return openslide.OpenSlide(self.wsi_path)
    
    def get_thumbnail(self):

        if not os.path.exists(self.save_path+'Thumbnail.png'):

            thumbnail_img = self.wsi.get_thumbnail((self.thumbnail_size,self.thumbnail_size))
            thumbnail_img.save(self.save_path+'Thumbnail.png')
            return thumbnail_img
        else:
            return Image.open(self.save_path+'Thumbnail.png')

    def read_xml(self,filepath):
        return ET.parse(filepath)
    
    def read_regions(self,region):

        Vertices = region.findall('./Vertices/Vertex')

        coords = []
        for Vertex in Vertices:
            x_coord = np.float32(Vertex.attrib['X'])
            y_coord = np.float32(Vertex.attrib['Y'])

            scaled_x_coord = int(self.thumbnail_size*(x_coord/self.wsi_width))
            scaled_y_coord = int(self.thumbnail_size*(y_coord/self.wsi_height))
            coords.append((scaled_x_coord,scaled_y_coord))

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
                group_count+=1

                # Get polygon bounds, compare to group current bounds ( converting to thumbnail dimensions )
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

        return intersect_spots, intersect_barcodes

    def gen_patch_centers(self,bbox):

        square_dim = int(self.patch_n**0.5)
        box_max_y, box_max_x = bbox[3], bbox[2]

        min_x, max_x = int(self.patch_size/2), box_max_x-int(self.patch_size/2)
        min_y, max_y = int(self.patch_size/2), box_max_y-int(self.patch_size/2)

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

    def generate_thumbnail_map(self):

        # Initializing all thumbnail heatmaps
        cell_mask = np.empty((self.thumbnail_size,self.thumbnail_size))
        cell_mask[:] = np.nan
        thumbnail_heatmap = {}
        for c in self.cell_names:
            thumbnail_heatmap[c] = cell_mask.copy()

        # Finding total spot bounding box
        tissue_bbox = []
        for s_g in self.spots:
            if not tissue_bbox == []:
                new_mins = [np.minimum(tissue_bbox[i],list(s_g['box'].bounds)[i]) for i in range(0,2)]
                new_maxs = [np.maximum(tissue_bbox[i],list(s_g['box'].bounds)[i]) for i in range(2,4)]

                tissue_bbox = new_mins+new_maxs
            else:
                tissue_bbox = list(s_g['box'].bounds)

        # Finding equally distributed spots over bounding box area
        patch_center_points = self.gen_patch_centers(tissue_bbox)
        # Iterating through patch centers and generating a box

        for center in tqdm(patch_center_points,desc='Iterating through patches'):
            patch_poly = shapely.geometry.box(
                center[0]-int(self.patch_size/2),
                center[1]-int(self.patch_size/2),
                center[0]+int(self.patch_size/2),
                center[1]+int(self.patch_size/2)
            )
            # Finding intersecting spots with this patch
            intersecting_spots, intersecting_barcodes = self.find_intersecting_spots(patch_poly)
            
            if len(intersecting_spots)>0:
                # Getting counts data for all intersecting spots
                intersect_counts = self.counts_data['main_cell_types'][self.counts_data['main_cell_types'].index.isin(intersecting_barcodes)]
                patch_bounds = list(patch_poly.bounds)
                patch_bounds = [int(i) for i in patch_bounds]

                all_sums = intersect_counts.sum(axis=0)
                for c in self.cell_names:
                    cell_val = all_sums.loc[c]

                    intermed_mask = np.empty((self.thumbnail_size,self.thumbnail_size))
                    intermed_mask[:] = np.nan
                    intermed_mask[patch_bounds[1]:patch_bounds[3],patch_bounds[0]:patch_bounds[2]] = cell_val

                    intermed_mask[intermed_mask==0] = np.nan
                    thumbnail_heatmap[c] = np.nanmean(np.stack((thumbnail_heatmap[c],intermed_mask),axis=-1),axis=-1)

        for c in self.cell_names:
            cell_heatmap = thumbnail_heatmap[c].copy()
            cell_heatmap[np.where(np.isnan(cell_heatmap))] = 0

            color_cell_heatmap = np.uint8(255*self.color_map(np.uint8(255*cell_heatmap))[:,:,0:3])
            zero_mask = np.where(cell_heatmap==0,0,255)
            cell_mask_4d = np.concatenate((color_cell_heatmap,zero_mask[:,:,None]),axis=-1)
            cell_vis = Image.fromarray(np.uint8(cell_mask_4d)).convert('RGBA')

            cell_vis.save(self.save_path+f'{c.replace("/","")}_thumbnail_vis.png')

def main():

    slide_name = 'XY04_IU-21-020F.svs'
    
    base_dir = '/mnt/c/Users/Sam/Desktop/HIVE/'
    slide_path = base_dir+'FFPE/'+slide_name
    spot_path = base_dir+'FFPE/Spot_Coordinates_large/'+slide_name.replace('.svs','_Large.xml')
    counts_path = base_dir+'counts_data/FFPE/CellTypeFractions_SpotLevel/V10S15-103_'+slide_name.replace('.svs','_cellfract.csv')
    counts_def_path = base_dir+'counts_data/Cell_SubTypes_Grouped.csv'
    save_path = os.getcwd()+'/Thumbnail_Masks/'+slide_name.replace('.svs','/')

    thumbnail_generator = Thumbnail_Heatmap_Generator(slide_path,spot_path,counts_path,counts_def_path,save_path)



if __name__=='__main__':

    main()

