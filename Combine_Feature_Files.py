"""


Combining the feature files from Nick



"""

import pandas as pd
import os

from glob import glob


def main():

    # Directory with all the feature files separated by FTU
    base_dir = 'C:\\Users\\Sam\\Desktop\\HIVE\\FFPE\\Feature Files\\'
    ftus = ['arteries','gloms','s_gloms','tubs']

    slides = ['XY01_IU-21-015F','XY02_IU-21-016F','XY03_IU-21-019F','XY04_IU-21-020F']
    
    for s in slides:
        slide_df = pd.DataFrame()
        for f in ftus:
            ftu_feats = pd.read_csv(base_dir+s+'_'+f+'.csv',index_col=0)
            if slide_df.empty:
                slide_df = ftu_feats
            else:
                slide_df = pd.concat([slide_df,ftu_feats],axis=0,ignore_index=True).fillna(0)

        slide_df.to_csv(base_dir+s+'_AllFeatures.csv',index=False)



if __name__=='__main__':
    main()


















