# -*- coding: utf-8 -*-
"""
Created on Wed Mar 12 09:44:07 2025

@author: Jackson Vaughn
"""

import os
import numpy as np
import pandas as pd
import seaborn as sns
from PIL import Image
import matplotlib.pyplot as plt
from matplotlib import colormaps
import matplotlib.colors as mcolors
from sc_functions import read_dfl_boolean_section, load_bathymetry, plot_bathymetry

#####################################################################################################################################
#User input variables

#Bathymetry data for wherever you need mapped
bathymetry = r'C:\Users\kw24171\OneDrive - University of Bristol\Desktop\Seawater_Data\Galap_bathymetry\ne_10m_bathymetry_all'
#Set to False if you don't want to add a map :(
addmap = True
#A max min lat lons for your map
lon_min, lon_max = -92.5, -88.5
lat_min, lat_max = -2, 1.4

#Path to the Atlantis core descriptions excel - you CANNOT have the excel open when you run this
excel_path = r"C:\Users\kw24171\OneDrive - University of Bristol\Desktop\Sediment_Cores\AT50_09_Core_descriptions.xlsx"

#Path to main FleXRay folder
folder_path = r'C:\Users\kw24171\OneDrive - University of Bristol\Desktop\Sediment_Cores\SC_Data\FleXRay'


#Where you want to save the files
output_path = r'C:\Users\kw24171\OneDrive - University of Bristol\Desktop\Sediment_Cores\SC_Images'

#A list of elements. If you want to plot the raw counts put ('element') e.g. ('Ba')
#If you want to plot the element normalized to another element, put ('element','normalizing_element') e.g. ('Ba'/'Ca')
elements_norms = [('Ba'), ('Ba','Al'), ('Mn','Al'), ('Fe','Al'), ('Al'), ('S')]




######################################################################################################################################
#Code
if addmap == True:
    import cartopy
    import cartopy.crs as ccrs
    from cartopy.mpl.geoaxes import GeoAxes
    from mpl_toolkits.axes_grid1.inset_locator import inset_axes
#must also change size for multicores

df_excel = pd.read_excel(excel_path, sheet_name='Preserved cores', header=3)

sample_ID = 'AL5157_PC01'
aaaa = df_excel['Core name'].values
depth_sw = df_excel.loc[df_excel['Core name'].values == sample_ID, 'Water depth (m)']

xpos_col = 'X Position (mm)'

#Extract elements and normalizations into seperate lists - appends False if no norm
elements_norms = [e if isinstance(e, tuple) else (e,) for e in elements_norms]
elements = [e[0] for e in elements_norms]
norms = [e[1] if len(e) > 1 else False for e in elements_norms]

# print(f"  ├── {subfolder}") # unhash if lost :(
for folder in os.listdir(folder_path):
    folder_path1 = os.path.join(folder_path, folder)  #Full path of subfolder
    if os.path.isdir(folder_path1):  #Check if it's a directory
        #Get sample ID from folder. It is a bit tedious since there is a 0 missing from the mcs in excel -
        #you will also need to rename some of your folders where you have a typo: you put AL5172_PC01 instead of AL5172-PC01 for example
        try:
            sample_ID = folder.split('_')[1].split('-')[1].replace('0','')
            sample_ID = folder.split('_')[1].split('-')[0]+'_'+sample_ID
        except:
            print('Fix your folder name typos')
        for subfolder in os.listdir(folder_path1):
            subfolder_path = os.path.join(folder_path1, subfolder)
            if os.path.isdir(subfolder_path):  #Check if it's a directory
            
                for ifile in os.listdir(subfolder_path): #Get necessary file paths
                    if 'settings.dfl' in ifile:
                        br_settings = os.path.join(subfolder_path,ifile)
                    elif 'batch result.txt' in ifile:
                        br_file = os.path.join(subfolder_path,ifile)
                    elif '.tif' in ifile and 'LowRes' not in ifile and 'Card' not in ifile:
                        img_path = os.path.join(subfolder_path,ifile)
                        
                #Get core image, crop and rotate
                img = Image.open(img_path)
                rotated_img = img.rotate(-90, expand=True)  #expand=True prevents cropping after rotation
                width, height = rotated_img.size
                top = 0
                bottom = height #Keep full height
                if 'PC' in sample_ID:
                    crop_amount = int(0.17 * width) #Crop ~17% off both sides, for push cores
                    left = crop_amount
                    right = width - crop_amount
                    cropped_img = rotated_img.crop((left, top, right, bottom))
                else:
                    crop_amount = int(0.05 * width) #Crop ~5% off both sides, for muc cores - could do even less really
                    left = crop_amount
                    right = width - crop_amount
                    cropped_img = rotated_img.crop((left, top, right, bottom))
                img_array = np.asarray(cropped_img) #Convert to NumPy array to display
                        
                #Get depth, lat, lon from excel    
                depth_sw = df_excel.loc[df_excel['Core name'].values == sample_ID, 'Water depth (m)'].iloc[0]
                lat = df_excel.loc[df_excel['Core name'].values == sample_ID, 'Latitude (N)'].iloc[0]
                lon = df_excel.loc[df_excel['Core name'].values == sample_ID, 'Longitude (E)'].iloc[0]
                
                
                #Unpack XRF data to dataframe
                with open(br_file, "r", encoding="utf-8") as file:
                    lines = file.readlines()
                #Find the header line (the first line containing tab-separated values)
                for i, line in enumerate(lines):
                    if "filename" in line:
                        header_index = i
                        break
                df = pd.read_csv(br_file, sep="\t", skiprows=header_index, encoding="utf-8") #Extract to dataframe
                df = df.loc[np.array(df['validity']==1), :] #Remove invalid measurements
                df = df.reset_index(drop=True)

                min_depth = df.loc[0,xpos_col]
                df[xpos_col] = df[xpos_col]-min_depth #Get depth from surface of sediment, not top of core container


                cmap = colormaps['tab10'] #Check documentation for colorblind colormaps, there are many. Find qualitative one. Can also make your own with just a list of colors e.g. ['red', 'blue', etc.] - check documentation for available colors
                color_list = [mcolors.to_hex(cmap(i)) for i in range(cmap.N)]
                colors_element_dict = dict(zip(elements,color_list))
                
                element_data = read_dfl_boolean_section(br_settings) #Create a dictionary according to the settings file, where valid elements are True and invalid elements are False
                
                
                #Plotting figure
                fig, ax = plt.subplots(1, len(elements)+1, figsize=(len(elements)*3,7))

                for count, element in enumerate(elements): 
                    if element_data.get(element) == True: #Checking the settings file
                        c = colors_element_dict.get(element) #Get color
                        if norms[count] == 'Rc1': #Ba excess according to I.L. Hendy. 0.0075 average crustal Ba/Al ratio according to Dymond et al., 1992
                            # col = df[element]-((0.0075)*df['Al'])
                            # col = (df[element]/(df['Al']+1))
                            # col = df[element]-df['Al']
                            col = df['Ba'].rolling(window=8).corr(df['Al']+1)#cov
                            #Standardize between zero and 1 and subtract col?
                            col = (df['Ba']/(df['Al']+0.1)) * (1-abs(col))
                            ax[count].set_xlabel(f'{element}'+'$_{RCv=12}$', color=c, fontsize=14)
                        elif norms[count] == 'Rc2': 
                            col = df['Ba'].rolling(window=12).corr(df['Al'])#cov
                            col = df['Ba']/(df['Al']+0.1) *abs(col)
                            ax[count].set_xlabel(f'{element}'+'$_{RC=15 abs}$', color=c, fontsize=14)
                        elif norms[count] == 'Rc3': 
                            col = df['Ba'].rolling(window=16).corr(df['Al'])#cov
                            col = df['Ba']/(df['Al']+0.1) *abs(col)
                            ax[count].set_xlabel(f'{element}'+'$_{RC=12}$', color=c, fontsize=14)
                        elif norms[count]: #Normalize
                            col = df[element]/df[norms[count]]
                            ax[count].set_xlabel(element+'/'+norms[count], color=c, fontsize=14)
                        else: #Or don't
                            col = df[element]
                            ax[count].set_xlabel(element+' counts', color=c, fontsize=14)
                            
                        ax[count].plot(col, df[xpos_col], color=c)
                        sns.despine(ax=ax[count], offset=10, bottom=False, trim=True) #Despining looks fancy
                        ax[count].tick_params(axis='x', which='both', bottom=True, labelbottom=True)
                        
                        if count != 0: 
                            ax[count].spines["left"].set_visible(False)
                            ax[count].get_yaxis().set_ticks([])
                        else:#Add the left spine only for the first plot
                            ax[count].tick_params(axis='y', which='both', left=True)
                            ax[count].set_ylabel('Depth (mm) \n (from sediment surface)', fontsize=14)

                        ax[count].set_ylim(min(df[xpos_col]), max(df[xpos_col])) #Make sure its expanded to match the image
                        ax[count].invert_yaxis()
                    else:
                        print(element+' not available in core '+sample_ID)
                        
                    
                
                
                ax[len(elements)].imshow(img_array, aspect='auto') #Plot the sediment core image
                ax[len(elements)].axis('off')
                
                if addmap == True: #Plot the map
                    inset_ax = inset_axes(ax[len(elements)], width="30%", height="20%", #The width height parameters don't really work like that. It's a pain in my butt (Even though I was the one who made the function :')
                                          loc='upper right',  # General location
                                          bbox_to_anchor=(-0.01, 0.1, 1, 1),  # Adjusted to be above the title
                                          bbox_transform=fig.transFigure,  # Position relative to the entire figure
                                          axes_class=GeoAxes,
                                          axes_kwargs=dict(projection=ccrs.PlateCarree()))
                    # Load the bathymetry data
                    depths_str, shp_dict = load_bathymetry(lon_min, lon_max, lat_min, lat_max, bathymetry)
                    # Plot base bathymetry map using function
                    mapp = plot_bathymetry(depths_str, shp_dict, lon_min, lon_max, lat_min, lat_max, inset_ax)
    
                    inset_ax.scatter([lon], [lat], color='red', s=30, zorder=5, transform=cartopy.crs.PlateCarree())

                #So as not to overwrite replicates
                if 'REP1' in folder:
                    sample_ID = sample_ID+'_01'
                elif 'REP2' in folder:
                    sample_ID = sample_ID+'_02'
                fig.suptitle(sample_ID+'\n'+str(depth_sw)+'m below sea level', fontsize=17)


                # plt.show() #Unhash this and hash the next line if you want to view the plots in your IDE without saving them
                plt.savefig(os.path.join(output_path,sample_ID+'-'+str(elements)), bbox_inches='tight') #Can add dpi=700 or some number like that if you want better quality
                plt.close()
                
                
                
                
#Should make so that it prints the percent complete after each iteration
