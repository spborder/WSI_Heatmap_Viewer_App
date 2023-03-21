"""

Reading RDS file in Python and extracting spot information


Can't get rpy2 to work locally, re-trying with just R implementation of these two functions

"""

import os
import sys

import numpy as np
import pandas as pd

from glob import glob
from tqdm import tqdm


#os.environ['R_HOME'] = '/mnt/c/Program Files/R/R-4.2.3'
#os.environ['PATH'] = '/mnt/c/Program Files/R/R-4.2.3/bin/x64'+':'+os.environ['PATH']
#print(os.environ['PATH'])

# Version 3.1.0 for Python 3.6
import rpy2
import rpy2.robjects.packages as rpackages

from rpy2.robjects.vectors import StrVector

class VisiumExtractor:
    def __init__(self):
        
        # Getting the R utility package and selecting CRAN mirror 
        self.utils = rpackages.importr('utils')
        self.utils.chooseCRANmirror(ind=1)

        # Checking if Seurat is installed in base R
        if not rpackages.isinstalled('Seurat'):
            self.utils.install_packages(('Seurat'))

        # Initializing R function to extract spot cell-type info
        rpy2.robjects.r('''
            # Function takes the RDS file name and the output file name as input
            ExtractCellType <- function(input_file, save_file) {
                
                read_vis_file <- readRDS(input_file)

                cell_type_fract <- GetAssayData(read_vis_file@assays[["predictions"]])

                # Normalizing on a per-column basis so cell-type fractions add to 1
                cell_type_fract <- cell_type_fract[1:nrow(cell_type_fract)-1,]
                cell_type_norm <- cell_type_fract/colSums(cell_type_fract)
                cell_type_norm[is.na(cell_type_norm)] = 0

                # Writing normalized cell type fractions to csv file
                write.csv(cell_type_norm,save_file)
                
            }
        
        ''')

        # Initializing R function to extract spot coordinates info
        rpy2.robjects.r('''
            # Function to extract spot coordinates from RDS file
            ExtractSpotCoords <- function(input_file,save_file) {
                
                # Reading visium rds file
                read_vis_file <- readRDS(input_file)

                # Getting spot barcodes and coordinates
                spot_coords <- read_vis_file@images[["slice1"]]@coordinates
                
                write.csv(spot_coords,save_file)
                
            }
        
        ''')

        self.get_spot_coords = rpy2.robjects.r['ExtractSpotCoords']
        self.get_spot_cell_types = rpy2.robjects.r['ExtractCellType']

    def process_directory(self,input_dir, output_dir):

        # Extract spot coordinates and cell types from a folder of RDS files
        rds_files = glob(input_dir+'*.RDS')

        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

        for f in tqdm(rds_files,desc='Processing Folder of RDS Files'):
            
            file_name = f.split(os.sep)[-1]
            output_spot_file = output_dir+file_name.replace('.RDS','_spot_coords.csv')
            output_cell_file = output_dir+file_name.replace('.RDS','_cell_fract.csv')

            self.get_spot_coords(f,output_spot_file)
            self.get_spot_cell_types(f,output_cell_file)

    def process_single(self,input_file,output_dir):

        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

        file_name = input_file.split(os.sep)[-1]
        output_spot_file = output_dir+file_name.replace('.RDS','_spot_coords.csv')
        output_cell_file = output_dir+file_name.replace('.RDS','_cell_fract.csv')

        print(f'Processing single file: {input_file}')
        self.get_spot_coords(input_file,output_spot_file)
        self.get_spot_cell_types(input_file,output_cell_file)


def main():


    # Testing if this works
    base_dir = '/mnt/c/Users/Sam/Desktop/HIVE/counts_data/'
    rds_single_file = base_dir+'Visium_Data/V10S15-103_XY01_IU-21-015F.RDS'
    
    output_dir = base_dir+'Test_Extract/'
    processor = VisiumExtractor()
    processor.process_single(rds_single_file,output_dir)


if __name__=='__main__':
    main()


