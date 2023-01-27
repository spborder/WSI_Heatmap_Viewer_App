"""

Pre-processing for SVS WSI image files for easy ingestion into WSI_Heatmap_Viewer.py application

preprocessing steps:
    - Extract and save thumbnail
        - pre-compute wsi-level cell distribution heatmaps and save overlays
    - Extract full resolution tiles of set size (10000,10000) or something like that
        - Can specify only tissue regions or get all
        - Record bbox coordinates for each tile and refer to those when the thumbnail is clicked

"""

import os
import sys
import numpy as np

from skimage.filters import gaussian, threshold_otsu
from PIL import Image

import openslide

slide_name = 'XY01_IU-21-015F.svs'
base_dir = '/mnt/c/Users/Sam/Desktop/HIVE/'
slide_path = base_dir+'FFPE/'+slide_name

save_path = base_dir+'FFPE/Slide_Tiles/'+slide_name.replace('.svs','')+'/'

tile_size = (1500,1500)


class TissueTileGenerator:
    def __init__(self,
                wsi,
                tile_size,
                save_path):

        self.wsi = wsi
        self.tile_size = tile_size
        self.save_path = save_path

        print('Generating tissue mask')
        self.tissue_mask = self.segment_tissue()

        if not os.path.exists(self.save_path):
            os.makedirs(self.save_path)

        if not os.path.exists(self.save_path+'tissue_mask.png'):
            Image.fromarray(np.uint8(255*self.tissue_mask)).save(self.save_path+'tissue_mask.png')

        self.all_tiles = self.get_all_tiles()

        self.tile_number = len(self.all_tiles)
        self.tile_n = 0

    def segment_tissue(self):

        # Mostly borrowed from HistomicsTK
        # Return mask of tissue regions in thumbnail image used to define where tissue tiles will be extracted
        tissue_mask = np.zeros((512,512))

        thumbnail_image = self.wsi.get_thumbnail((512,512))
        grayscale_thumb = 255 - np.mean(thumbnail_image,axis=-1)

        gauss_blurred = gaussian(grayscale_thumb,sigma=0.1, output=None, mode='nearest', preserve_range=True)
        thresh_val = threshold_otsu(gauss_blurred[gauss_blurred>0])

        grayscale_thumb[grayscale_thumb<thresh_val] = 0

        tissue_mask = 0 + (grayscale_thumb>0)

        return tissue_mask
        
    def get_all_tiles(self):
        poly_list = []


        return poly_list

    def __iter__(self):


        return self
    
    def __next__(self):
        if self.tile_n<self.tile_number:

            next_tile = self.all_tiles[self.tile_n]
            self.tile_n+=1
            return next_tile
        
        else:
            raise StopIteration
        



def main(slide_path,save_path,tile_size):

    wsi = openslide.OpenSlide(slide_path)

    tile_gen = TissueTileGenerator(wsi,tile_size,save_path)

    # Iterating through available tiles


if __name__=='__main__':


    main(slide_path,save_path,tile_size)









