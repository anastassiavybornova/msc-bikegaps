#%% packages to import

import pickle

import time
from datetime import date

import os
from pathlib import Path

import numpy as np
import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt
import math

from collections import Counter

import shapely
import random
import itertools
import networkx as nx
import igraph as ig

from haversine import haversine, Unit, haversine_vector

import folium
from folium import Circle

from IPython.display import display

import pyproj

from ast import literal_eval # to convert str into lists (imported osmid data)

import winsound

### IMPORT PACKAGES

import pickle

import time
from datetime import date

import os
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import random
import itertools

import networkx as nx
import igraph as ig

from haversine import haversine, Unit, haversine_vector

import folium
from folium import Circle

from IPython.display import display

import pyproj
from shapely.geometry import Point
import shapely
import shapely.ops as ops
import geopandas as gpd

import winsound

#%% Functions to import

# custom functions for tfm workflow
# last update: 2021-07-20

# define function that creates attribute dictionary for nodes and edges
# (for input to nx.add_edges_from/add_nodes_from)
def make_attr_dict(*args, **kwargs): 
    
    argCount = len(kwargs)
    
    if argCount > 0:
        attributes = {}
        for kwarg in kwargs:
            attributes[kwarg] = kwargs.get(kwarg, None)
        return attributes
    else:
        return None # (if no attributes are given)
    
# adding a cost factor to car edges on nw
def add_cost(nw, costfactor):
    for edge in nw.edges:
        edgetype = nw.edges[edge]["category_edge"]
        if edgetype != "car":
            nw.edges[edge]["cost"] = nw.edges[edge]["length"] # cost of bikeable edge is its length
        elif edgetype == "car":
            # cost added because of type
            nw.edges[edge]["cost"] = nw.edges[edge]["length"] * costfactor     
            # i.e. cost of car edge is increased prop. to its length and the cost factor
    return nw

# function that finds length distribution among edge types
def type_distr(keys_edge, dict_lengths, dict_types):

    clen = 0
    blen = 0

    for edge in keys_edge:

        mylen = dict_lengths[edge]
        mytype = dict_types[edge]

        if mytype == "car":
            clen += mylen
        else:
            blen += mylen
    mylens = [clen, blen]
    return mylens

# function to get "bike percentage" of path from igraph object (returns value in %):
def get_bp_ig(nw, edgelist):
        tl = np.sum([nw.es[edge]["length"] for edge in edgelist]) # total length = sum of all edge lengths
        bl = np.sum([nw.es[edge]["length"] for edge in edgelist if nw.es[edge]["category_edge"]!="car"]) # bike length = sum of all non-car edges
        return np.round(100*bl/tl, 3)

# returns fraction of path length on test nw vs. path length on base nw;
# if nodes are not connected on test nw, returns positive infinity 
import networkx as nx
def path_on_other_nw_nx(testnw, origin, destination, pathlength):
    try:
        return(nx.dijkstra_path_length(testnw, origin, destination, "length")/pathlength)
    except:
        return(math.inf) 


# guessing by deviation whether nodes are actually connected by "parallel" bike edges;
# flexible tolerace (with e.g. 1.1 = 110%)
# parallel bikepath guess
import math
def pb_guess(perc, tolerance):
    if perc > tolerance:
        return(0)
    else:
        return(1)
def common_entries(*dcts):
    """Like zip() but for dicts.
    See: https://stackoverflow.com/questions/16458340/python-equivalent-of-zip-for-dictionaries
    """
    if not dcts:
        return
    for i in set(dcts[0]).intersection(*dcts[1:]):
        yield (i,) + tuple(d[i] for d in dcts)


# EDITED trying to clean up lonlat confusion
def project_nxpos_edit(G, map_center = False):
    """Take a spatial nx network G and projects its GPS coordinates to local azimuthal.
    Returns transformed positions, as used by nx.draw()
    """
    
    nodedict = nx.get_node_attributes(G, "coord")
    nids = nodedict.keys()
    lons = {nid:[c for c in nodedict[nid].coords][0][0] for nid in nids} # longitude (1st position in geopandas Points coord, "x" in OSM data; +12 for cph)
    lats = {nid:[c for c in nodedict[nid].coords][0][1] for nid in nids} # latitude (2nd position in geopandas Points coord, "y" in OSM data; +54 for cph)

    pos = {nid:(lon,lat) for (nid,lon,lat) in common_entries(lons,lats)}
    
    if map_center:
        loncenter = map_center[0]
        latcenter = map_center[1]
    else:
        loncenter = np.mean(list(lons.values())) # centering of longitude around mean of LON values (not LAT values)
        latcenter = np.mean(list(lats.values())) # centering of latitude around mean of LAT values (not LON values)
    
    local_azimuthal_projection = "+proj=aeqd +R=6371000 +units=m +lon_0={} +lat_0={}".format(loncenter, latcenter)
    
    wgs84_to_aeqd = pyproj.Transformer.from_proj(
        proj_from = pyproj.CRS('EPSG:4326'), # crs epsg 4326 is "WGS84" system
        proj_to = pyproj.Proj(local_azimuthal_projection), # as defined above
    always_xy = True) # !! set always_xy to True to define lonlat order for all objects (otherwise WGS84 assumes yx/latlon)
    
    pos_transformed = {nid:list(ops.transform(wgs84_to_aeqd.transform, Point(lonlat)).coords)[0] for nid, lonlat in pos.items()}
    
    return pos_transformed, (loncenter,latcenter)

def nxdraw(G_nx, networktype, plotparam, map_center = False, nnids = False, drawfunc = "nx.draw", nodesize = 0, weighted = False, maxwidthsquared = 0, simplified = False):
    
    if nnids is not False: # Restrict to nnids node ids
        #nnids_nx = [k for k,v in dict(G_nx.nodes(data=True)).items() if v['id'] in nnids]
        G_nx = G_nx.subgraph(nnids) # !! modified
        
    pos_transformed, map_center = project_nxpos_edit(G_nx, map_center)
    
    if weighted is True:
        # The max width should be the node diameter (=sqrt(nodesize))
        widths = list(nx.get_edge_attributes(G_nx, "weight").values())
        widthfactor = 1.1 * math.sqrt(maxwidthsquared) / max(widths)
        widths = [max(0.33, w * widthfactor) for w in widths]
        eval(drawfunc)(G_nx, pos_transformed, **plotparam[networktype], node_size = nodesize, width = widths)
    elif type(weighted) is float or type(weighted) is int and weighted > 0:
        eval(drawfunc)(G_nx, pos_transformed, **plotparam[networktype], node_size = nodesize, width = weighted)
    else:
        eval(drawfunc)(G_nx, pos_transformed, **plotparam[networktype], node_size = nodesize)
    return map_center

def initplot(plotparam):
    fig = plt.figure(figsize=(plotparam["bbox"][0]/plotparam["dpi"], plotparam["bbox"][1]/plotparam["dpi"]), dpi=plotparam["dpi"])
    plt.axes().set_aspect('equal')
    plt.axes().set_xmargin(0.01)
    plt.axes().set_ymargin(0.01)
    plt.axes().set_axis_off()
    return fig


#%% Parameters to set

# parameters for plotting (folium / nx)from packages import *

# create basemaps dict
basemaps = {
    'Google Maps': folium.TileLayer(
                tiles = 'https://mt1.google.com/vt/lyrs=m&x={x}&y={y}&z={z}',
                attr = 'Google',
                name = 'Google Maps',
                overlay = True,
                control = True
        ),
    'Google Satellite': folium.TileLayer(
                tiles = 'https://mt1.google.com/vt/lyrs=s&x={x}&y={y}&z={z}',
                attr = 'Google',
                name = 'Google Satellite',
                overlay = True,
                control = True
        ),
    'st_lines': folium.TileLayer(
        tiles = 'https://stamen-tiles.a.ssl.fastly.net/toner-lines/{z}/{x}/{y}.png',
        attr = 'Map tiles by <a href="http://stamen.com">Stamen Design</a>, <a href="http://creativecommons.org/licenses/by/3.0">CC BY 3.0</a> &mdash; Map data &copy; <a href="https://www.openstreetmap.org/copyright">OpenStreetMap</a> contributors',
        name = 'Stamen Toner Lines',
        overlay = True,
        control = True),
    'st_lite': folium.TileLayer(
        tiles = 'https://stamen-tiles.a.ssl.fastly.net/toner-lite/{z}/{x}/{y}.png',
        attr = 'Map tiles by <a href="http://stamen.com">Stamen Design</a>, <a href="http://creativecommons.org/licenses/by/3.0">CC BY 3.0</a> &mdash; Map data &copy; <a href="https://www.openstreetmap.org/copyright">OpenStreetMap</a> contributors',
        name = 'Stamen Toner Lite',
        overlay = True,
        control = True),
    'st_toner': folium.TileLayer(
        tiles = 'https://stamen-tiles.a.ssl.fastly.net/toner/{z}/{x}/{y}.png',
        attr = 'Map tiles by <a href="http://stamen.com">Stamen Design</a>, <a href="http://creativecommons.org/licenses/by/3.0">CC BY 3.0</a> &mdash; Map data &copy; <a href="https://www.openstreetmap.org/copyright">OpenStreetMap</a> contributors',
        name = 'Stamen Toner',
        overlay = True,
        control = True),
    'st_terrain': folium.TileLayer(
        tiles = 'https://stamen-tiles.a.ssl.fastly.net/terrain/{z}/{x}/{y}.png',
        attr = 'Map tiles by <a href="http://stamen.com">Stamen Design</a>, <a href="http://creativecommons.org/licenses/by/3.0">CC BY 3.0</a> &mdash; Map data &copy; <a href="https://www.openstreetmap.org/copyright">OpenStreetMap</a> contributors',
        name = 'Stamen Terrain',
        overlay = True,
        control = True),
    'osm': folium.TileLayer(
        tiles = "openstreetmap", 
        name = "OpenStreetMap",
        control = True, 
        overlay = True),
    'Blueprint-copy': folium.TileLayer(
        tiles = 'https://api.mapbox.com/styles/v1/krktalilu/ckrdjkf0r2jt217qyoai4ndws/tiles/256/{z}/{x}/{y}@2x?access_token=pk.eyJ1Ijoia3JrdGFsaWx1IiwiYSI6ImNrcmRqMXdycTB3NG8yb3BlcGpiM2JkczUifQ.gEfOn5ttzfH5BQTjqXMs3w',
        name = "whiteprint",
        attr = 'Mapbox',
        control = True,
        overlay = True
        )

}

plotparam = {
    "bbox": (1280,1280),
    "dpi": 96,
    
    "caredges": {"width": 0.5, "edge_color": '#999999'},
    "bikeedges": {"width": 0.5, "edge_color": '#2222ff'},
    
#    "biketrack_offstreet": {"width": 0.75, "edge_color": '#00aa22'},
#    "bikeable": {"width": 0.75, "edge_color": '#222222'},
#    "bikegrown": {"width": 3.75, "edge_color": '#0eb6d2', "node_color": '#0eb6d2'},
    
#    "highlight_biketrack": {"width": 3.75, "edge_color": '#2222ff', "node_color": '#2222ff'},
#    "highlight_bikeable": {"width": 3.75, "edge_color": '#222222', "node_color": '#222222'},
    
#    "poi_unreached": {"node_color": '#ff7338', "edgecolors": '#ffefe9'},
#    "poi_reached": {"node_color": '#0b8fa6', "edgecolors": '#f1fbff'},
    
#    "abstract": {"edge_color": '#000000', "alpha": 0.75},
    
#    "blcc": {"width": 1.25, "edge_color": '#33d1ff'}, # changed color to lighter blue
#    "clcc": {"width": 0.5, "edge_color": '#bb4126'}, # changed color dark red

#    "blcc_path": {"width": 1.25, "edge_color": '#33d1ff'}, # path connecting d1 nodes - pink
#    "clcc_path": {"width": 3, "edge_color": '#ffc733'}, # path connecting d1 nodes - orange
    
    "bikenodes": {"width": 0, "edge_color": '#0eb6d2', "node_color": '#b226bb'}, # to plot bike nodes (pink)
    "carnodes": {"width": 0, "edge_color": '#0eb6d2', "node_color": '#000000'}, # to plot car nodes (black)
    "bothnodes": {"width": 0, "edge_color": '#0eb6d2', "node_color": '#ffa53f'}, # to plot multi nodes (orange)
    "d2nodes": {"width": 0, "edge_color": '#0eb6d2', "node_color": '#ff3349'} # to plot d2 nodes (red)
    
}

gapcolors = {}
gapcolors["ml"] = "red"
gapcolors["br"] = "orange"
gapcolors["ra"] = "brown"
gapcolors["rt"] = "#ff00f7"
gapcolors["is"] = "yellow"
gapcolors["di"] = "black"
gapcolors["pe"] = "black"
gapcolors["x"] = "grey"
