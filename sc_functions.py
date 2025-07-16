# -*- coding: utf-8 -*-
"""
Created on Wed Mar 12 12:57:48 2025

@author: kw24171
"""

def read_dfl_boolean_section(file_path, section="[Periodic table]"):
    periodic_table = {}  # Dictionary to store elements and their boolean values
    in_section = False  # Track whether in the target section
    
    with open(file_path, "r", encoding="utf-8") as file:
        for line in file:
            line = line.strip()

            # Detect section start
            if line == section:
                in_section = True
                continue
            
            # Exit section when encountering a new section
            if in_section and line.startswith("[") and line.endswith("]"):
                break

            # Parse element and boolean values
            if in_section and "=" in line:
                element, value = line.split("=")
                element = element.strip()
                value = value.strip().upper() == "TRUE"  # Convert to boolean
                periodic_table[element] = value

    return periodic_table




def load_bathymetry(ln_min, ln_max, lt_min, lt_max, directory="path_to_your_local_files"):
    """Load bathymetry shapefiles from a local directory."""
    import cartopy.io.shapereader as shpreader
    from glob import glob
    import numpy as np

    # Read shapefiles, sorted by depth
    shp_dict = {}
    files = glob(f'{directory}/*.shp')  # Adjust path if necessary
    assert len(files) > 0, "No shapefiles found in the specified directory."
    files.sort()
    depths = []
    
    for f in files:
        depth = '-' + f.split('_')[-1].split('.')[0]  # Extract depth from file name
        depths.append(depth)
        bbox = (ln_min, lt_min, ln_max, lt_max)  # Define bounding box
        nei = shpreader.Reader(f, bbox=bbox)  # Load shapefile with specified bounding box
        shp_dict[depth] = nei

    depths = np.array(depths)[::-1]  # Sort from surface to bottom
    return depths, shp_dict



#Plot the bathymetry and return the axes to continue to plot in the main file
def plot_bathymetry(depths_str, shp_dict, lon_min, lon_max, lat_min, lat_max, ax):
    import cartopy.feature as cfeature
    #import matplotlib.pyplot as plt
    import cartopy.crs as ccrs
    import matplotlib
    import matplotlib.patches as patches
    import matplotlib.path as mpath
    # Construct a discrete colormap with colors corresponding to each depth
    depths = depths_str.astype(int)
    N = len(depths)
    nudge = 0.01  # shift bin edge slightly to include data
    boundaries = [min(depths)] + sorted(depths+nudge)  # low to high
    norm = matplotlib.colors.BoundaryNorm(boundaries, N)
    blues_cm = matplotlib.colormaps['Blues_r'].resampled(N)
    colors_depths = blues_cm(norm(depths))

    # Set up plot
    #subplot_kw = {'projection': ccrs.LambertCylindrical()}
    #fig, ax = plt.subplots(subplot_kw=subplot_kw, figsize=(9, 7))
    ax.set_extent([lon_min, lon_max, lat_min, lat_max], crs=ccrs.PlateCarree())  # x0, x1, y0, y1

    # Iterate and plot feature for each depth level
    for i, depth_str in enumerate(depths_str):
        ax.add_geometries(shp_dict[depth_str].geometries(),
                          crs=ccrs.PlateCarree(),
                          color=colors_depths[i])

    # Add standard features
    ax.add_feature(cfeature.LAND, color='grey')
    ax.coastlines(lw=1, resolution='110m')
    gridlines = ax.gridlines(draw_labels=True, color='black', alpha=0.1, linestyle='-')
    ax.set_position([0.03, 0.05, 0.8, 0.9])

    gridlines.bottom_labels = False
    gridlines.right_labels = True
    gridlines.left_labels = False 
    gridlines.top_labels = True
    
    # Set font size for labels
    gridlines.xlabel_style = {'size': 12} 
    gridlines.ylabel_style = {'size': 12} 


    #Plot zebra borders
    # Define alternating colors
    colors = ['black', 'white']
    line_width = 7  # Adjust line width as needed
    
    # Generate zebra stripes for the bottom edge (left to right)
    for lon in range(int(lon_min), int(lon_max)):
        segment_vertices = [[lon, lat_min], [lon + 1, lat_min]]
        segment_codes = [mpath.Path.MOVETO, mpath.Path.LINETO]
        segment_path = mpath.Path(segment_vertices, segment_codes)
        color = colors[lon % len(colors)]
        segment_patch = patches.PathPatch(segment_path, facecolor='none', ec=color, lw=line_width, zorder=100)
        ax.add_patch(segment_patch)
    
    # Generate zebra stripes for the right edge (bottom to top)
    for lat in range(int(lat_min), int(lat_max)):
        segment_vertices = [[lon_max, lat], [lon_max, lat + 1]]
        segment_codes = [mpath.Path.MOVETO, mpath.Path.LINETO]
        segment_path = mpath.Path(segment_vertices, segment_codes)
        color = colors[lat % len(colors)]
        segment_patch = patches.PathPatch(segment_path, facecolor='none', ec=color, lw=line_width, zorder=100)
        ax.add_patch(segment_patch)
    
    # Generate zebra stripes for the top edge (right to left)
    for lon in range(int(lon_max), int(lon_min), -1):
        segment_vertices = [[lon, lat_max], [lon - 1, lat_max]]
        segment_codes = [mpath.Path.MOVETO, mpath.Path.LINETO]
        segment_path = mpath.Path(segment_vertices, segment_codes)
        color = colors[lon % len(colors)]
        segment_patch = patches.PathPatch(segment_path, facecolor='none', ec=color, lw=line_width, zorder=100)
        ax.add_patch(segment_patch)
    
    # Generate zebra stripes for the left edge (top to bottom)
    for lat in range(int(lat_max), int(lat_min), -1):
        segment_vertices = [[lon_min, lat], [lon_min, lat - 1]]
        segment_codes = [mpath.Path.MOVETO, mpath.Path.LINETO]
        segment_path = mpath.Path(segment_vertices, segment_codes)
        color = colors[lat % len(colors)]
        segment_patch = patches.PathPatch(segment_path, facecolor='none', ec=color, lw=line_width, zorder=100)
        ax.add_patch(segment_patch)
        
    
    
    
    # # Add custom colorbar
    # axi = fig.add_axes([0.85, 0.1, 0.025, 0.8])
    # ax.add_feature(cfeature.BORDERS, linestyle=':')
    # sm = plt.cm.ScalarMappable(cmap=blues_cm, norm=norm)
    # fig.colorbar(mappable=sm,
    #              cax=axi,
    #              spacing='proportional',
    #              extend='min',
    #              ticks=depths,
    #              label='Depth (m)')

    #need to add coordinate system

    # Convert vector bathymetries to raster (saves a lot of disk space)
    # while leaving labels as vectors
    ax.set_rasterized(True)
    
    return ax

