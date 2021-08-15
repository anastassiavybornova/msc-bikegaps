#!/usr/bin/env python
# coding: utf-8

# In[1]:


# last update: 2021-06-25
# find gap sequences with cut-off length 1200


# In[2]:


# *make notebook trusted
# import packages and custom functions
get_ipython().run_line_magic('run', 'ppf.py')


# In[3]:


### import data

# nx
H = nx.read_gpickle("../data/pickle/H.gpickle") # read in H network (simplified from full data)
H_typedict_edge = nx.get_edge_attributes(H, "category_edge")
H_typedict_node = nx.get_node_attributes(H, "category_node")
H_lendict = nx.get_edge_attributes(H, "length")
H_coorddict_edge = nx.get_edge_attributes(H, "coord")
H_coorddict_node = nx.get_node_attributes(H, "coord")
B = nx.read_gpickle("../data/pickle/B.gpickle") # bnw

# ig
h = ig.Graph.Read_Pickle("../data/pickle/h.pickle")
b = ig.Graph.Read_Pickle("../data/pickle/b.pickle") # bnw

# conversion tables
nids_conv = pd.read_pickle("../data/pickle/nids_conv.pickle")
eids_conv = pd.read_pickle("../data/pickle/eids_conv.pickle")

# od is list of randomly sampled OD (origin-destination) pairs  
# od = pd.read_pickle("../data/pickle/odf.pickle")

# ebc df
ebc = pd.read_pickle("../data/pickle/ebc.pickle")


# In[4]:


get_ipython().run_cell_magic('time', '', '\n### get a list of "gaps" - deciding on cut-off length of 500m\n\ngap_cutoff = 1200 # look for car-edge stretchtes of a max length of 500m\n\n# APSP car sequences with non-car source & target NODES:\n\n# get all paths of specified max length into a dict\nmypaths = dict(nx.all_pairs_dijkstra_path(H, cutoff = gap_cutoff, weight = "length"))  \n\ncar_stretches = [] # initialize list of "car stretches"\nnodes_origin = []\nnodes_destination = []\n\nfor key_source in mypaths.keys(): # loop through all source nodes for paths in dict\n    \n    for key_target in mypaths[key_source].keys(): # loop through all target nodes for paths\n\n        if len(mypaths[key_source][key_target])>1: # if node list has at least 2 elements = 1 edge: \n            \n            myedges = [tuple(sorted(edge)) for edge in zip(mypaths[key_source][key_target], \n                                                           mypaths[key_source][key_target][1:])] # make edge ids from node list\n            mytypes = [H_typedict_edge[edge] for edge in myedges] # read out types of the edges from dict\n            if np.all([mytype == "car" for mytype in mytypes]) \\\n            and H_typedict_node[key_source] != "car" \\\n            and H_typedict_node[key_target] != "car":\n                car_stretches.append(myedges) # if yes: add the "gap" edges to car_stretches list\n                nodes_origin.append(key_source)\n                nodes_destination.append(key_target)\n                \n# store all found gaps in pandas df and process them there:\nmygaps = pd.DataFrame({"edges_nx" : car_stretches,\n                      "o_nx": nodes_origin,\n                      "d_nx": nodes_destination}) # make pd df to use apply functions\nmygaps["o_ig"] = mygaps.apply(lambda x: int(nids_conv[nids_conv["nx"]==x.o_nx]["ig"]), axis = 1)\nmygaps["d_ig"] = mygaps.apply(lambda x: int(nids_conv[nids_conv["nx"]==x.d_nx]["ig"]), axis = 1)\n\nmygaps["edges_ig"] = mygaps.apply(lambda x: [int(eids_conv[eids_conv["nx"]==edge_nx]["ig"]) for edge_nx in x.edges_nx], axis = 1)\nmygaps["edges_nx_set_str"] = mygaps.apply(lambda x: str(sorted(set(x.edges_nx))), axis = 1) # make string of sorted sets \nmygaps["edges_ig_set_str"] = mygaps.apply(lambda x: str(sorted(set(x.edges_ig))), axis = 1) # make string of sorted sets \n# (to find dupli/tripli...cated gaps - getting edges in the same order)\nmygaps = mygaps[-mygaps.duplicated(subset = ["edges_nx_set_str", "edges_ig_set_str"], keep = "first")].reset_index(drop = True)\nmygaps = mygaps.drop(columns = ["edges_nx_set_str", "edges_ig_set_str"]) # remove str col used for duplicates finding\nmygaps["length"] = mygaps.apply(lambda x: np.sum([H_lendict[edge] for edge in x.edges_nx]), axis = 1)\nmygaps = mygaps.sort_values(by = "length").reset_index(drop = True) # sort by length and...\nmygaps["gapid"] = mygaps.index # add gap ids\n\n# add length info for edges in gap\nmygaps["lengths"] = mygaps.apply(lambda x: [h.es[ig_id]["length"] for ig_id in x.edges_ig], axis = 1)\n# add ebc info for edges in gap (4 different factors) # ebc index = ig id!\nmygaps["ebcs_length"] = mygaps.apply(lambda x: [ebc.loc[ig_id, "ebc_length"] for ig_id in x.edges_ig], axis = 1)\nmygaps["ebcs_cutoff"] = mygaps.apply(lambda x: [ebc.loc[ig_id, "ebc_cutoff"] for ig_id in x.edges_ig], axis = 1)\nmygaps["ebcs_cost"] = mygaps.apply(lambda x: [ebc.loc[ig_id, "ebc_cost"] for ig_id in x.edges_ig], axis = 1)\nmygaps["ebcs_r"] = mygaps.apply(lambda x: [ebc.loc[ig_id, "ebc_r"] for ig_id in x.edges_ig], axis = 1)\n\nmygaps["ms_length"] = mygaps.apply(lambda x: np.sum(np.array(x.lengths) * np.array(x.ebcs_length)), axis = 1)\nmygaps["ms_length_n"] = mygaps["ms_length"] / mygaps["length"] # meters on network per meter gap (investment...)\nmygaps["ms_r"] = mygaps.apply(lambda x: np.sum(np.array(x.lengths) * np.array(x.ebcs_r)), axis = 1)\nmygaps["ms_r_n"] = mygaps["ms_r"] / mygaps["length"] # meters on network per meter gap (investment...)\n\n# print info on gaps found\nprint("gap_cutoff", gap_cutoff, "m")\nprint(len(car_stretches), "total paths")\nprint(len(mygaps), "unique gaps")\n\n# investigating "parallel" bike//gap edges:\nmygaps["connected_otherwise"] = mygaps.apply(lambda x: path_on_other_nw_nx(B, x.o_nx, x.d_nx, x.length), axis = 1)\nmygaps = mygaps[mygaps["connected_otherwise"]>1.5].reset_index(drop = True) # remove gaps that are connected otherwise with a factor <=1.5\n# save mygaps df to pickle\nmygaps.to_pickle("../data/pickle/mygaps.pickle")\n\nwinsound.Beep(300,700)')


# In[6]:


# make df of top1000 gaps for ebc within radius normed to length
mygaps_ttrn = mygaps.sort_values(by = "ms_r_n", ascending = False).reset_index(drop = True)[0:1000]

# save as pickle (for later use) and as csv (for ma - manual analysis)
mygaps_ttrn.to_pickle("../data/pickle/mygaps_ttrn.pickle")
mygaps_ttrn[["gapid", "connected_otherwise"]].to_csv("../analysis/top1000/mygaps_ttrn_ma.csv", index = False)


# # from here on: folium plots

# In[16]:


# plotting each gap separately (by ID) - 100gaps per plot
# Plotting with folium

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
        )
}

### add citizen input
# conversion form UTM32 to lonlat
mycrs = 4326

# os.listdir("../data/kk/citizen_input")
ci = gpd.read_file("../data/kk/citizen_input/citizen input - cykelstiprioriteringsplan.shp")
ci["geometry"] = ci["geometry"].to_crs(mycrs)
ci["locations"] = ci.apply(lambda x: [(point.y, point.x) for point in x.geometry], axis = 1)
ci["locations"] = ci.apply(lambda x: x.locations[0], axis = 1)
# add citizen input
ci1_fg = folium.FeatureGroup(name = "CI Cykelsti mangler", show = False)
ci2_fg = folium.FeatureGroup(name = "CI Cykelsti er for smal", show = False)
ci3_fg = folium.FeatureGroup(name = "CI Signalregulerede kryds med stor trængsel", show = False)


# make ci geometries
for i in range(len(ci["locations"])):
    if ci["input_kate"][i] == "Cykelsti mangler":
        folium.Circle(location = ci["locations"][i], radius = 1, color = "#f333ff").add_to(ci1_fg)
    if ci["input_kate"][i] == "Cykelsti er for smal":
        folium.Circle(location = ci["locations"][i], radius = 1, color = "#33f3ff").add_to(ci2_fg)
    if ci["input_kate"][i] == "Signalregulerede kryds med stor trængsel":
            folium.Circle(location = ci["locations"][i], radius = 1, color = "#9633ff").add_to(ci3_fg)
            
# make H network geometries
mynw_fg = folium.FeatureGroup(name = "H", show = False)
# car-only edges
cnw = folium.FeatureGroup(name = "H car edges", show = True)
# bike edges
bnw = folium.FeatureGroup(name = "H bike edges", show = True)

hloc = []
cloc = []
bloc = []

for edge in H.edges:
    myloc = [(c[1], c[0]) for c in H_coorddict_edge[tuple(sorted(edge))].coords]
    hloc.append(myloc)
    if H_typedict_edge[edge] == "car":
        cloc.append(myloc)
    if H_typedict_edge[edge] != "car":
        bloc.append(myloc)
            
mynw_line = folium.PolyLine(locations = hloc, weight = 2, color = "black").add_to(mynw_fg)
cnw_line = folium.PolyLine(locations = cloc, weight = 2, color = "grey").add_to(cnw)
bnw_line = folium.PolyLine(locations = bloc, weight = 2, color = "blue").add_to(bnw)
    
count = 0

for b in range(10):
    
    # b for bins - making bins of 100 gaps each
    mygaps_plot = mygaps_ttrn.loc[range(100*b,100*(b+1))].copy().reset_index(drop = True)

    # create map object
    m = folium.Map(location=[55.6761, 12.5683], zoom_start = 12) # cph

    # add satellite to m
    basemaps['Google Satellite'].add_to(m)
    
    # add citizen input
    ci1_fg.add_to(m)
    ci2_fg.add_to(m)
    ci3_fg.add_to(m)

    # add our networks
    mynw_fg.add_to(m)
    cnw.add_to(m)
    bnw.add_to(m)
    
    
    for i in mygaps_plot.index:
        my_fg = folium.FeatureGroup(name = "Gap " + str(mygaps_plot.iloc[i]["gapid"]), show = False)
        my_seqs = [H_coorddict_edge[edge].coords for edge in mygaps_plot.iloc[i]["edges_nx"]]
        my_coords = [[tuple([c[1], c[0]]) for c in seq] for seq in my_seqs] # unpack coord seqs
        my_line = folium.PolyLine(locations = my_coords, 
                    weight = 3, 
                    color = "red").add_to(my_fg)
        my_fg.add_to(m)
        
    marker_fg = folium.FeatureGroup(name = "Gap markers", show = False)
    for i in mygaps_plot.index:
        myloc = [(c[1], c[0]) for c in H_coorddict_node[mygaps_plot["o_nx"][i]].coords][0]
        myid = mygaps_plot["gapid"][i]
        folium.map.Marker(location=myloc, popup=myid).add_to(marker_fg)

    folium.LayerControl().add_to(m)

    m.save("../analysis/top1000/2021-06-25-gapids-" + str(b) + "-hide.html")


# In[17]:


# plotting each gap separately (by ID) - 100gaps per plot
# Plotting with folium

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
        )
}

### add citizen input
# conversion form UTM32 to lonlat
mycrs = 4326

# os.listdir("../data/kk/citizen_input")
ci = gpd.read_file("../data/kk/citizen_input/citizen input - cykelstiprioriteringsplan.shp")
ci["geometry"] = ci["geometry"].to_crs(mycrs)
ci["locations"] = ci.apply(lambda x: [(point.y, point.x) for point in x.geometry], axis = 1)
ci["locations"] = ci.apply(lambda x: x.locations[0], axis = 1)
# add citizen input
ci1_fg = folium.FeatureGroup(name = "CI Cykelsti mangler", show = False)
ci2_fg = folium.FeatureGroup(name = "CI Cykelsti er for smal", show = False)
ci3_fg = folium.FeatureGroup(name = "CI Signalregulerede kryds med stor trængsel", show = False)


# make ci geometries
for i in range(len(ci["locations"])):
    if ci["input_kate"][i] == "Cykelsti mangler":
        folium.Circle(location = ci["locations"][i], radius = 1, color = "#f333ff").add_to(ci1_fg)
    if ci["input_kate"][i] == "Cykelsti er for smal":
        folium.Circle(location = ci["locations"][i], radius = 1, color = "#33f3ff").add_to(ci2_fg)
    if ci["input_kate"][i] == "Signalregulerede kryds med stor trængsel":
            folium.Circle(location = ci["locations"][i], radius = 1, color = "#9633ff").add_to(ci3_fg)
            
# make H network geometries
mynw_fg = folium.FeatureGroup(name = "H", show = False)
# car-only edges
cnw = folium.FeatureGroup(name = "H car edges", show = True)
# bike edges
bnw = folium.FeatureGroup(name = "H bike edges", show = True)

hloc = []
cloc = []
bloc = []

for edge in H.edges:
    myloc = [(c[1], c[0]) for c in H_coorddict_edge[tuple(sorted(edge))].coords]
    hloc.append(myloc)
    if H_typedict_edge[edge] == "car":
        cloc.append(myloc)
    if H_typedict_edge[edge] != "car":
        bloc.append(myloc)
            
mynw_line = folium.PolyLine(locations = hloc, weight = 2, color = "black").add_to(mynw_fg)
cnw_line = folium.PolyLine(locations = cloc, weight = 2, color = "grey").add_to(cnw)
bnw_line = folium.PolyLine(locations = bloc, weight = 2, color = "blue").add_to(bnw)
    
count = 0

for b in range(10):
    
    # b for bins - making bins of 100 gaps each
    mygaps_plot = mygaps_ttrn.loc[range(100*b,100*(b+1))].copy().reset_index(drop = True)

    # create map object
    m = folium.Map(location=[55.6761, 12.5683], zoom_start = 12) # cph

    # add satellite to m
    basemaps['Google Satellite'].add_to(m)
    
    # add citizen input
    ci1_fg.add_to(m)
    ci2_fg.add_to(m)
    ci3_fg.add_to(m)

    # add our networks
    mynw_fg.add_to(m)
    cnw.add_to(m)
    bnw.add_to(m)
    
    

    for i in mygaps_plot.index:
        my_fg = folium.FeatureGroup(name = "Gap " + str(mygaps_plot.iloc[i]["gapid"]), show = True)
        my_seqs = [H_coorddict_edge[edge].coords for edge in mygaps_plot.iloc[i]["edges_nx"]]
        my_coords = [[tuple([c[1], c[0]]) for c in seq] for seq in my_seqs] # unpack coord seqs
        my_line = folium.PolyLine(locations = my_coords, 
                    weight = 3, 
                    color = "red").add_to(my_fg)
        my_fg.add_to(m)
        # add markers with gap ids
        
    marker_fg = folium.FeatureGroup(name = "Gap markers", show = False)
    for i in mygaps_plot.index:
        myloc = [(c[1], c[0]) for c in H_coorddict_node[mygaps_plot["o_nx"][i]].coords][0]
        myid = mygaps_plot["gapid"][i]
        folium.map.Marker(location=myloc, popup=myid).add_to(marker_fg)

    marker_fg.add_to(m)

    folium.LayerControl().add_to(m)

    m.save("../analysis/top1000/2021-06-25-gapids-" + str(b) + "-show.html")

