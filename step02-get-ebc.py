#!/usr/bin/env python
# coding: utf-8

# # Workflow step 1: make ebc df from H/h
# 
# last update: 2021-06-22
# 
# * input: H.gpickle (nx) and h.pickle (ig) from step 00 "makeH"
# * output: ebc data frame (ebc.pickle) with ebcs calculated on length, cost

# In[1]:


# import packages and custom functions
get_ipython().run_line_magic('run', 'ppf.py')


# In[2]:


# read in H network (simplified from full data)
H = nx.read_gpickle("../data/pickle/H.gpickle")
h = ig.Graph.Read_Pickle("../data/pickle/h.pickle")

nids_conv = pd.read_pickle("../data/pickle/nids_conv.pickle")
eids_conv = pd.read_pickle("../data/pickle/eids_conv.pickle")

# extract edge and node attributes as dictionaries
typedict = nx.get_edge_attributes(H, "category_edge")
lendict = nx.get_edge_attributes(H, "length")
ncoorddict = nx.get_node_attributes(H, "coord")

H_coorddict_edge = nx.get_edge_attributes(H, "coord")


# In[3]:


get_ipython().run_cell_magic('time', '', 'ebc = pd.DataFrame({"edge_ig": [e.index for e in h.es]})\nebc["edge_nx"] = ebc.apply(lambda x: eids_conv.at[x.edge_ig, "nx"], axis = 1)\nebc["length"] = ebc.apply(lambda x: h.es[x.edge_ig]["length"], axis = 1)\nebc["ebc_length"] = h.edge_betweenness(directed = False, cutoff = None, weights = "length")\nebc["ebc_cutoff"] = h.edge_betweenness(directed = False, cutoff = 2500, weights = "length")\nebc["ebc_cost"] = h.edge_betweenness(directed = False, cutoff = None, weights = "cost")\nebc')


# In[4]:


get_ipython().run_cell_magic('time', '', '\n# https://igraph.org/python/doc/api/igraph._igraph.GraphBase.html#get_all_shortest_paths\n# ebc on shortest paths in certain radius:\n\nr = 2500 # define radius\ncrs_lonlat = 4326 # coordinate systems for conversion \ncrs_utm32 = 32632 # coordinate systems for conversion\n\n# get node coordinates (as geopandas df... for euclidean distance calc)\nmynodeids = [node for node in H.nodes]\nmynodes = [ncoorddict[node] for node in H.nodes]\nndf = gpd.GeoDataFrame(geometry = mynodes, crs = crs_lonlat)\nndf["node_nx"] = mynodeids\nndf["node_ig"] = ndf.apply(lambda x: int(nids_conv[nids_conv["nx"]==x.node_nx]["ig"]), axis = 1)\nndf["geometry_meters"] = ndf["geometry"].to_crs(crs_utm32)\nndf["points"] = ndf.apply(lambda x: x.geometry_meters.coords[0], axis = 1)\nndf = ndf.drop(columns = ["geometry", "geometry_meters"])\nnodes_array = np.array(list(ndf["points"])) # nodes array')


# In[5]:


get_ipython().run_cell_magic('time', '', '\nndict = {}\n# nids_conv["ig"]\nfor i in range(len(ndf)):\n    \n    node_ig = ndf.loc[i, "node_ig"]\n    node_nx = ndf.loc[i, "node_nx"]\n    \n    ndict[node_ig] = {} # dictionary by key = igraph node\n\n    # distance calculations\n    node_start = np.array(ndf.loc[i, "points"])\n    distances = np.sqrt(np.sum((nodes_array - node_start)**2, axis = 1))\n    inrange_list_ig = list(ndf[distances < r]["node_ig"])\n    inrange_list_nx = list(ndf[distances < r]["node_nx"])\n    # save results: inrange_array and its length\n    ndict[node_ig]["nodes_id_nx"] = inrange_list_nx\n    ndict[node_ig]["nodes_id_ig"] = inrange_list_ig\n    ndict[node_ig]["nodes_nr"] = len(inrange_list_ig)')


# In[6]:


get_ipython().run_cell_magic('time', '', '\nbiglist = [] # initialize list\n# for each node,\nfor i in range(len(ndf)):\n    for node_ig in ndict[i]["nodes_id_ig"]: # look at its neighbours within r (as calculated in ndict)\n        mytuple = tuple(sorted([i, node_ig]))\n        biglist.append(mytuple) # save node-neighbour tuple pairs\nbiglist = [item for item in biglist if item[0]!=item[1]] # drop tuples that contain same node twice\nbiglist = list(dict.fromkeys(biglist)) # drop duplicates (od vs. do)\nlen(biglist)\nwinsound.Beep(350, 700)')


# In[9]:


get_ipython().run_cell_magic('time', '', 'ebcr = pd.DataFrame({"od": biglist})\ndel(biglist) # free memory')


# In[10]:


get_ipython().run_cell_magic('time', '', '# should take approx. 2h\nebcr["path"] = ebcr.apply(lambda x: h.get_shortest_paths(v = x.od[0], \n                                                         to = x.od[1], \n                                                         weights = "length", \n                                                         output = "epath"), \n                          axis = 1)')


# In[16]:


get_ipython().run_line_magic('time', 'ebcr["pathlist"] = ebcr.apply(lambda x: [item for sublist in x.path for item in sublist], axis = 1)')


# In[17]:


winsound.Beep(400, 7000)


# In[18]:


get_ipython().run_line_magic('time', 'hugelist = [item for sublist in ebcr["pathlist"] for item in sublist]')


# In[19]:


from collections import Counter


# In[20]:


get_ipython().run_line_magic('time', 'ebcr_dict = Counter(hugelist)')


# In[21]:


get_ipython().run_line_magic('time', 'ebc["ebc_r"] = ebc.apply(lambda x: ebcr_dict[x.edge_ig], axis = 1)')


# In[29]:


ebc


# In[24]:


ebc.to_pickle("../data/pickle/ebc.pickle")

