#!/usr/bin/env python
# coding: utf-8

# # Workflow step 0: make H (simplified network)
# last update: 2021-06-17
# * input: osm of bike and car nodes and edges as csv (using geopandas!)
# * outputs: simplified nw H (nx) and h (ig), subgraph B (nx) and b (ig) - without car edges ("bikeable"), all pickled

# In[2]:


# import packages and custom functions
get_ipython().run_line_magic('run', 'ppf.py')


# In[4]:


### DEFINE OVERALL PARAMETERS
path_data = "..\data" # define path to data folder


# In[5]:


### IMPORT CSV FILES

### importing data step 1: get folder & subfolders 
myfolders = os.listdir(path_data) # subfolders with csv files inside

folder_be = os.path.join(path_data, "copenhagen_biketrack_edges")
folder_bn = os.path.join(path_data, "copenhagen_biketrack_nodes")
folder_ce = os.path.join(path_data, "copenhagen_carall_edges")
folder_cn = os.path.join(path_data, "copenhagen_carall_nodes")

# defining filenames (=paths to files)
path_be = os.path.join(folder_be, os.listdir(folder_be)[0])
path_bn = os.path.join(folder_bn, os.listdir(folder_bn)[0])
path_ce = os.path.join(folder_ce, os.listdir(folder_ce)[0])
path_cn = os.path.join(folder_cn, os.listdir(folder_cn)[0])

# read in data

# bike edges
be = pd.read_csv(path_be).drop(columns = ["key", "lanes", "name", "highway", "maxspeed", "bridge", "tunnel", "junction", "width", "access", "ref", "service", "area"])
be["geometry"] = be.apply(lambda x: shapely.wkt.loads(x.geometry), axis = 1)
be = gpd.GeoDataFrame(be, geometry = "geometry") 

# car edges
ce = pd.read_csv(path_ce).drop(columns = ["key", "lanes", "name", "highway", "maxspeed", "bridge", "tunnel", "junction", "width", "access", "ref", "service"])
ce["geometry"] = ce.apply(lambda x: shapely.wkt.loads(x.geometry), axis = 1)
ce = gpd.GeoDataFrame(ce, geometry = "geometry") 

# bike nodes
bn = pd.read_csv(path_bn).drop(columns = ["highway", "ref"])
bn["geometry"] = bn.apply(lambda x: shapely.wkt.loads(x.geometry), axis = 1)
bn = gpd.GeoDataFrame(bn, geometry = "geometry")

# car nodes
cn = pd.read_csv(path_cn).drop(columns = ["highway", "ref"])
cn["geometry"] = cn.apply(lambda x: shapely.wkt.loads(x.geometry), axis = 1)
cn = gpd.GeoDataFrame(cn, geometry = "geometry")

del(path_be, path_bn, 
    path_ce, path_cn, 
    folder_be, folder_bn, 
    folder_ce, folder_cn,
    myfolders)


# In[6]:


### NODES: 
# make df with all nodes (an) to pass it to nx
an = pd.merge(bn, cn, how = "outer", indicator = True) # merging
an["type"] = an["_merge"].cat.rename_categories(["bike", "car", "multi"]) # adding info on type
an = an.drop(columns = "_merge")
an = an.sort_values(by = "osmid").reset_index(drop = True) # sort by osmid
an["attr_dict"] = an.apply(lambda x: make_attr_dict(category_node = x.type, coord = x.geometry), axis = 1) # add attr_dict


# In[7]:


### EDGES: 
# make df with all edges (ae) to pass it to nx

# add edge ids (strings with "id1, id2" sorted (id1 < id2))
be["edge_id"] = be.apply(lambda x: str(sorted([x["u"], x["v"]])), axis = 1)
ce["edge_id"] = ce.apply(lambda x: str(sorted([x["u"], x["v"]])), axis = 1)
# (edge ids are set as strings; converting back: with "from ast import literal_eval" fct)
# finding duplicates by ["osmid", "oneway", "edge_id", "length"]

# simplifying network into undirected - beu and ceu contain the "undirected" edges
# (removing all parallel edges)

beu = be.drop_duplicates(subset = ["osmid", "oneway", "edge_id", "length"],
                  keep = "first",
                  inplace = False,
                  ignore_index = True).copy()
ceu = ce.drop_duplicates(subset = ["osmid", "oneway", "edge_id", "length"],
                  keep = "first",
                  inplace = False,
                  ignore_index = True).copy()

# add type info prior to merging
beu["type"] = "bike"
ceu["type"] = "car"

# concatenate
ae = pd.concat([beu, ceu]).reset_index(drop = True)

ae.loc[ae.duplicated(subset = ["u", "v", "osmid", "oneway", "length", "edge_id"], keep = False), "type"] = "multi"

ae = ae.drop_duplicates(subset = ["u", "v", "osmid", "oneway", "length", "edge_id", "type"], 
                          keep = "first",
                          ignore_index = True, 
                          inplace = False)

ae_tokeep = ae[ae.duplicated("edge_id", keep = False) & (ae["type"]=="bike")].index
ae_todrop = ae[ae.duplicated("edge_id", keep = False) & (ae["type"] == "car")].index

ae.loc[ae_tokeep, "type"] = "multi"
ae = ae.drop(ae_todrop)

# add attribute dictionary (for nx)
ae["attr_dict"] = ae.apply(lambda x: make_attr_dict(length = x.length, 
                                                    category_edge = x.type,
                                                    edge_id = x.edge_id,
                                                    coord = x.geometry,
                                                    intnodes = []), # intnodes attribute: for storing simplification info on interstitial nodes 
                             axis = 1)

# sort by "left" node (id1 < id2 - to control order of tuple keys in nx)
ae["order"] = ae.apply(lambda x: np.min([x["u"], x["v"]]), axis = 1)
ae = ae.sort_values(by = "order").reset_index(drop = True)
ae["orig"] = ae.apply(lambda x: np.min([x["u"], x["v"]]), axis = 1)
ae["dest"] = ae.apply(lambda x: np.max([x["u"], x["v"]]), axis = 1)
ae = ae.drop(columns = ["order", "u", "v"]) # instead of "u" and "v",
# we will use "origin" and "destination" where osmid(origin) < osmid (destination)!

del(ae_todrop, ae_tokeep, beu, ceu)

print("ae df: " + str(len(ae)))
print("unique edge ids: " + str(len(np.unique(ae["edge_id"]))))


# In[8]:


### MERGING BY NX

# bicycle network
bnw = nx.Graph()
bnw.add_nodes_from(an[an["type"]!="car"].loc[:,["osmid", "attr_dict"]].itertuples(index = False))
bnw.add_edges_from(ae[ae["type"]!="car"].loc[:,["orig", "dest", "attr_dict"]].itertuples(index = False))

# car network
cnw = nx.Graph()
cnw.add_nodes_from(an[an["type"]!="bike"].loc[:,["osmid", "attr_dict"]].itertuples(index = False))
cnw.add_edges_from(ae[ae["type"]!="bike"].loc[:,["orig", "dest", "attr_dict"]].itertuples(index = False))

# make multinetwork containing ALL edges
mnw = nx.Graph()
mnw.add_nodes_from(an.loc[:,["osmid", "attr_dict"]].itertuples(index = False))
mnw.add_edges_from(ae.loc[:,["orig", "dest", "attr_dict"]].itertuples(index = False))


# In[9]:


get_ipython().run_cell_magic('time', '', '\n# KEEPING ONLY THE LARGEST CONNECTED COMPONENT ON MNW\n\ncd_nodeset = []\n\nfor comp in nx.connected_components(mnw):\n    cd_nodeset = cd_nodeset + [comp]\n\nn = len(cd_nodeset)\n\nprint("number of disconnected components on mnw: " + str(n))\n\ncd_size = [None]*n\ncd_network = [None]*n\ncd_coord_dict = [None]*n\ncd_coord_list = [None]*n\ncd_types = [None]*n\n\nfor i in range(n):\n    cd_size[i] = len(cd_nodeset[i])\n    cd_network[i] = nx.subgraph(mnw, cd_nodeset[i])\n    cd_coord_dict[i] = nx.get_edge_attributes(cd_network[i], "coord")\n    cd_coord_list[i] = [cd_coord_dict[i][key] for key in cd_coord_dict[i].keys()]\n    cd_types[i] = nx.get_edge_attributes(cd_network[i], "category_edge")\n    \ncomps = pd.DataFrame({\n    \'nodeset\': cd_nodeset, \n    \'size\': cd_size,\n    \'network\': cd_network,\n    \'coord\': cd_coord_list,\n    \'type\': cd_types\n             })\n\ndel(cd_nodeset, cd_size, cd_network, cd_coord_list, cd_types, cd_coord_dict)\n\nlcc = np.max(comps["size"])\nprint("size of lcc: " + str(lcc))\n\ncomps = comps.sort_values(by = "size", ascending = False).reset_index(drop = True)')


# In[10]:


# DEFINE MNWL as largest connected component
# (drop all others)
mnwl_nodes = comps["nodeset"][0]
mnwl_edges = ae.loc[ae.apply(lambda x: x.orig in mnwl_nodes, axis = 1),:].copy().reset_index(drop = True)
mnwl = nx.subgraph(mnw, mnwl_nodes)


# In[9]:


# Plotting disc comp on folium:
m = folium.Map(location=[55.6761, 12.5683], zoom_start = 12) 

# get node coordinate dict
mnw_node_coorddict = nx.get_node_attributes(mnw, "coord")

# plot all disconnected components with n_nodes > 1
comp_lcc = folium.FeatureGroup(name = "largest connected comp", show = False)
comp_disc = folium.FeatureGroup(name = "disconnected comps", show = False)
comp_node = folium.FeatureGroup(name = "disconnected nodes", show = False)

# largest connected component
lloc_ls = [ls.coords for ls in comps["coord"][0]] # unpack linestrings
lloc = [[tuple([c[1], c[0]]) for c in seq] for seq in lloc_ls] # unpack coord seqs

# disconnected (smaller) component coordinates
dloc = []
# disconnected (node) coordinates
nloc = []

for i in range(1, 72):
    my_ls = [ls.coords for ls in comps["coord"][i]] # unpack linestrings
    my_cs = [[tuple([c[1], c[0]]) for c in seq] for seq in my_ls] # unpack coord seqs
    dloc = dloc + my_cs
for i in range(72, 77):
    my_node = comps["nodeset"][i].copy()
    my_c = [c for c in mnw_node_coorddict[my_node.pop()].coords][0]
    my_c = tuple([my_c[1], my_c[0]]) # reverse order
    nloc.append(my_c)

lcc_line = folium.PolyLine(locations = lloc, weight = 1, color = "grey").add_to(comp_lcc)
disc_line = folium.PolyLine(locations = dloc, weight = 2, color = "red").add_to(comp_disc)

# single nodes with d=0 as Circles: 
for loc in nloc:
    Circle(location = loc, radius = 3, color = "purple").add_to(comp_node)
    
comp_lcc.add_to(m)
comp_disc.add_to(m)
comp_node.add_to(m)

folium.LayerControl().add_to(m)

#display(m)

m.save("../plots/00f_disccomp.html")


# In[61]:


# Plotting disc comp on folium:
m = folium.Map(location=[55.6761, 12.5683], zoom_start = 12) 

# get node coordinate dict
mnw_node_coorddict = nx.get_node_attributes(mnw, "coord")
H_node_coorddict = nx.get_node_attributes(H, "coord")

# plot all disconnected components with n_nodes > 1
e = folium.FeatureGroup(name = "edge", show = True)
n_end = folium.FeatureGroup(name = "nodes before", show = True)
n_int = folium.FeatureGroup(name = "nodes after", show = True)

e_simp = [key for key in H_intnode_dict.keys() if len(H_intnode_dict[key])==16][0]
e_loc = [(c[1], c[0]) for c in H.edges[e_simp]["coord"].coords]
e_line = folium.PolyLine(locations = e1_loc, weight = 2, color = "blue").add_to(e)

# nodes 
n_end_loc = []
n_int_loc = []
for n in e_simp:
    n_end_loc.append([(c[1], c[0]) for c in H_node_coorddict[n].coords][0])
for n in H.edges[e_simp]["intnodes"]:
    n_int_loc.append([(c[1], c[0]) for c in mnw_node_coorddict[n].coords][0])

for loc in n_end_loc:
    Circle(location = loc, radius = 3, color = "black").add_to(n_end)
for loc in n_int_loc:
    Circle(location = loc, radius = 1, color = "red").add_to(n_int)
    
e.add_to(m)
n_end.add_to(m)
n_int.add_to(m)

folium.LayerControl().add_to(m)

#display(m)

# m.save("../plots/00f_intnodes.html")


# In[11]:


del(comps, mnwl_nodes, mnw)


# In[13]:


mnwl_typedict = nx.get_edge_attributes(mnwl, "category_edge")


# In[15]:


get_ipython().run_cell_magic('time', '', '\n# SIMPLIFICATION\n\n# make a copy of mnwl - H will be simplified and manipulated throughout while loop\nH = mnwl.copy()\n\n# set parameters for the while loop\nsimplify_further = True\nrun = 0\n\n# loop runs while there are interstitial nodes on the nw\nwhile simplify_further:\n    \n    run += 1\n    print("Run " + str(run) + ", " + time.ctime())\n    \n    # get all nodes from nw\n    points_all_list = sorted(list(H.nodes))\n\n    # get all node degrees\n    degrees_all_list = [None]*len(points_all_list)\n    for i in range(len(points_all_list)):\n        degrees_all_list[i] = H.degree(points_all_list[i])\n\n    # make df with node + degree info + remove (T/F) + types (of incident edges)\n    pointsall = pd.DataFrame({\n        "osmid": points_all_list, \n        "d": degrees_all_list, \n        "remove": None, \n        "types": None})\n    \n    # get edge attributes as dict\n    catdict = nx.get_edge_attributes(H, "category_edge")\n    \n    # get edge type information (car/bike/multi) from attribute dictionary\n    pointsall["types"] = pointsall.apply(lambda x: \n                                         [ catdict[tuple(sorted(edge))] for edge in H.edges(x.osmid) ], \n                                         axis = 1)\n\n    # split df in "endpoints" and d2 nodes\n    pointsend = pointsall[pointsall["d"]!=2].copy().reset_index(drop = True)\n    pointsd2 = pointsall[pointsall["d"]==2].copy().reset_index(drop = True)\n\n    # non-d2 nodes: all of them are remove=False (to keep)\n    pointsend.at[:,"remove"] = False\n    # d2 nodes: the ones that have same 2 edge types incident are remove=True\n    pointsd2["remove"] = pointsd2.apply(lambda x: x.types[0]==x.types[1], axis = 1)\n\n    # final result: 2 dfs - nodes_final and nodes_interstitial\n\n    # nodes_final = nodes to keep (either they have d!=2 or they have d==2 but 2 different edge types)\n    nodes_final = pd.concat([pointsend, pointsd2[pointsd2["remove"]==False].copy()]).reset_index(drop = True)\n\n    # nodes_interstitial = nodes to remove (d2 nodes with same 2 edge types incident)\n    nodes_interstitial = pointsd2[pointsd2["remove"]==True].copy().reset_index(drop = True)\n    nodes_interstitial["types"] = nodes_interstitial.apply(lambda x: x.types[0], axis = 1) # remove second-edge info (is same as first)\n\n    print(str(len(nodes_final)) + "+" + str(len(nodes_interstitial)) + "=" + str(len(H)) + " (final + interstitial = total?)")\n    print(len(nodes_final) + len(nodes_interstitial) == len(H)) # check - should be the same (give "True")\n\n    del(pointsall, catdict, degrees_all_list, points_all_list, pointsend, pointsd2)\n\n    # save info about endpoint/interstitial to node attributes on mnwl\n    for i in range(len(nodes_interstitial)):\n        H.nodes[nodes_interstitial.loc[i, "osmid"]]["category_point"] = "int"\n    for i in range(len(nodes_final)):\n        H.nodes[nodes_final.loc[i, "osmid"]]["category_point"] = "end"\n\n    # make df with interstitial edges\n    eint = nodes_interstitial.copy() \n    eint["orig"] = eint.apply(lambda x: sorted([n for n in H.neighbors(x.osmid)])[0], axis = 1)\n    eint["dest"] = eint.apply(lambda x: sorted([n for n in H.neighbors(x.osmid)])[1], axis = 1)\n\n    # add info on edge lengths\n    lendict = nx.get_edge_attributes(H, "length")\n    eint["length_new"] = eint.apply(lambda x: \n                                    np.sum(\n                                        [lendict[tuple(sorted(edge))] for edge in H.edges(x.osmid)]\n                                    ), \n                                    axis = 1)\n\n    stack = list(np.unique(eint["osmid"]))\n\n    print("edge len before: " + str(len(H.edges)))\n    print("stack len before: " + str(len(stack)))\n    \n    Hprior = H.copy() # make a copy of the nw in each simplification step\n    # to use for checking for neighbours for removing from stack\n    \n    # interstitial nodes dictionary - to keep track of nodes that are removed by "while stack"\n    intnodesdict = nx.get_edge_attributes(H, "intnodes")\n    # edge coordinate dictionary - to merge linestrings of aggregated edges\n    edgecoorddict = nx.get_edge_attributes(H, "coord")\n    \n    while stack:\n\n        mynode = stack.pop()\n        \n        for n in nx.neighbors(Hprior, mynode): # remove neighbors from ORIGINAL nw\n            if n in stack:\n                stack.remove(n)\n                #print("removed "+ str(n))\n                \n        # u and v are the neighbors of "mynode"\n        u = eint.loc[eint["osmid"]==mynode]["orig"].values[0]\n        v = eint.loc[eint["osmid"]==mynode]["dest"].values[0]\n        \n        # counter (to break out of loop if it is not increased)\n        nodes_removed = 0\n        \n        if (u,v) not in H.edges: # only if neighbors are not neighbors themselves - \n            # to avoid roundabouts from disappearing\n            \n            # get info on interstitional nodes (for deriving edge coordinates later on)\n            myintnodes = [intnodesdict[tuple(sorted(edge))] for edge in H.edges(mynode)]\n            myintnodes.append([mynode])\n            myintnodes = [x for x in list(itertools.chain.from_iterable(myintnodes)) if x]\n            \n            H.add_edge(u_of_edge = u,\n                        v_of_edge = v,\n                        length = eint.loc[eint["osmid"]==mynode]["length_new"].values[0],\n                        category_edge = eint.loc[eint["osmid"]==mynode]["types"].values[0],\n                        intnodes = myintnodes,\n                        edge_id = str(sorted([u, v])),\n                        coord = shapely.ops.linemerge( [ edgecoorddict[tuple(sorted([u,mynode]))],\n                                                         edgecoorddict[tuple(sorted([v,mynode]))] ]\n                                                     ) \n                      )\n\n            H.remove_node(mynode)\n            nodes_removed += 1\n            \n    print("edge len after: " + str(len(H.edges)))\n    print("stack len after: " + str(len(stack)))\n    winsound.Beep(350, 700)\n    \n    if nodes_removed == 0:\n        \n        simplify_further = False # to break out of loop\n        \n        # add cost factor to car edges\n        mycost = 1.25\n        H = add_cost(H, mycost)\n        \n        # save simplified network to H gpickle\n        nx.write_gpickle(H, "../data/pickle/H.gpickle") \n        \n        print("Done")')


# In[52]:


# make "bikeable" network from H (excluding car edges)
bikeable_nodes = [node for node in H.nodes if H.nodes[node]["category_node"]!="car"]
H_noncar_induced = H.subgraph(bikeable_nodes).copy() 
# induced subgraph - still contains the car edges that lie between multi nodes; - exclude them:
banw = H_noncar_induced.copy()
banw.remove_edges_from([edge for edge in banw.edges if banw.edges[edge]["category_edge"]=="car"])
nx.write_gpickle(banw, "../data/pickle/B.gpickle") 
       


# In[53]:


# conversion to igraph!
h = ig.Graph.from_networkx(H)
h.write_pickle("../data/pickle/h.pickle")
b = ig.Graph.from_networkx(banw)
b.write_pickle("../data/pickle/b.pickle")
# to read in again: Graph.Read_Pickle()

# eids: "conversion table" for edge ids from igraph to nx 
eids_nx = [tuple(sorted(literal_eval(h.es(i)["edge_id"][0]))) for i in range(len(h.es))]
eids_ig = [i for i in range(len(h.es))]
eids_conv = pd.DataFrame({"nx": eids_nx, "ig": eids_ig})

# nids: "conversion table" for node ids from igraph to nx
nids_nx = [h.vs(i)["_nx_name"][0] for i in range(len(h.vs))]
nids_ig = [i for i in range(len(h.vs))]
nids_conv = pd.DataFrame({"nx": nids_nx, "ig": nids_ig})

eids_conv.to_pickle("../data/pickle/eids_conv.pickle")
nids_conv.to_pickle("../data/pickle/nids_conv.pickle")


# In[142]:


# plot on folium - 3 cols
m = folium.Map(location=[55.6761, 12.5683], zoom_start = 11.5) 

H_edge_coorddict = nx.get_edge_attributes(H, "coord")
H_typedict = nx.get_edge_attributes(H, "category_edge")

# entire H network
mynw_fg = folium.FeatureGroup(name = "H", show = False)
# car-only edges
cnw = folium.FeatureGroup(name = "car edges", show = False)
# bike edges
bnw = folium.FeatureGroup(name = "bike edges", show = False)
# multi edges
multi_nw = folium.FeatureGroup(name = "multi edges", show = False)

hloc = []
cloc = []
bloc = []
mloc = []

for edge in H.edges:
    myloc = [(c[1], c[0]) for c in H_edge_coorddict[tuple(sorted(edge))].coords]
    hloc.append(myloc)
    if H_typedict[edge] == "car":
        cloc.append(myloc)
    if H_typedict[edge] == "bike":
        bloc.append(myloc)
    if H_typedict[edge] == "multi":
        mloc.append(myloc)

mynw_line = folium.PolyLine(locations = hloc, weight = 2, color = "black").add_to(mynw_fg)
cnw_line = folium.PolyLine(locations = cloc, weight = 2, color = "orange").add_to(cnw)
bnw_line = folium.PolyLine(locations = bloc, weight = 2, color = "blue").add_to(bnw)
multi_nw_line = folium.PolyLine(locations = mloc, weight = 2, color = "green").add_to(multi_nw)

mynw_fg.add_to(m)
cnw.add_to(m)
bnw.add_to(m)
multi_nw.add_to(m)


folium.LayerControl().add_to(m)

#display(m)

# save map
m.save("../plots/00f_H_simp_3col.html")
# del(m)
winsound.Beep(350, 700)


# In[73]:


# plot on folium - 2 cols
m = folium.Map(location=[55.6761, 12.5683], zoom_start = 11.5) 

H_edge_coorddict = nx.get_edge_attributes(H, "coord")
H_typedict = nx.get_edge_attributes(H, "category_edge")

# entire H network
mynw_fg = folium.FeatureGroup(name = "H", show = False)
# car-only edges
cnw = folium.FeatureGroup(name = "car edges", show = False)
# bike edges
bnw = folium.FeatureGroup(name = "bike edges", show = False)
# multi edges
multi_nw = folium.FeatureGroup(name = "multi edges", show = False)

hloc = []
cloc = []
bloc = []
mloc = []

for edge in H.edges:
    myloc = [(c[1], c[0]) for c in H_edge_coorddict[tuple(sorted(edge))].coords]
    hloc.append(myloc)
    if H_typedict[edge] == "car":
        cloc.append(myloc)
    if H_typedict[edge] == "bike":
        bloc.append(myloc)
    if H_typedict[edge] == "multi":
        mloc.append(myloc)

mynw_line = folium.PolyLine(locations = hloc, weight = 2, color = "black").add_to(mynw_fg)
cnw_line = folium.PolyLine(locations = cloc, weight = 2, color = "orange").add_to(cnw)
bnw_line = folium.PolyLine(locations = bloc, weight = 2, color = "blue").add_to(bnw)
multi_nw_line = folium.PolyLine(locations = mloc, weight = 2, color = "blue").add_to(multi_nw)

mynw_fg.add_to(m)
cnw.add_to(m)
bnw.add_to(m)
multi_nw.add_to(m)


folium.LayerControl().add_to(m)

#display(m)

# save map
m.save("../plots/00f_H_simp_2col.html")
# del(m)
winsound.Beep(350, 700)


# In[146]:


d2nodes = [node for node in H.nodes if H.degree[node]==2]


# In[148]:


# plot on folium - d2 nodes
m = folium.Map(location=[55.6761, 12.5683], zoom_start = 11.5) 

H_edge_coorddict = nx.get_edge_attributes(H, "coord")
H_typedict = nx.get_edge_attributes(H, "category_edge")

# entire H network
mynw_fg = folium.FeatureGroup(name = "H", show = False)
# car-only edges
cnw = folium.FeatureGroup(name = "car edges", show = False)
# bike edges
bnw = folium.FeatureGroup(name = "bike edges", show = False)
# multi edges
multi_nw = folium.FeatureGroup(name = "multi edges", show = False)

hloc = []
cloc = []
bloc = []
mloc = []

for edge in H.edges:
    myloc = [(c[1], c[0]) for c in H_edge_coorddict[tuple(sorted(edge))].coords]
    hloc.append(myloc)
    if H_typedict[edge] == "car":
        cloc.append(myloc)
    if H_typedict[edge] == "bike":
        bloc.append(myloc)
    if H_typedict[edge] == "multi":
        mloc.append(myloc)

mynw_line = folium.PolyLine(locations = hloc, weight = 3, color = "black").add_to(mynw_fg)
cnw_line = folium.PolyLine(locations = cloc, weight = 3, color = "orange").add_to(cnw)
bnw_line = folium.PolyLine(locations = bloc, weight = 3, color = "blue").add_to(bnw)
multi_nw_line = folium.PolyLine(locations = mloc, weight = 3, color = "green").add_to(multi_nw)

mynw_fg.add_to(m)
cnw.add_to(m)
bnw.add_to(m)
multi_nw.add_to(m)


nodes_fg = folium.FeatureGroup(name = "nodes", show = False)
nodes_loc = []
for n in d2nodes:
    nodes_loc.append([(c[1], c[0]) for c in H_node_coorddict[n].coords][0])
for loc in nodes_loc:
    Circle(location = loc, radius = 3, color = "black").add_to(nodes_fg)
    
nodes_fg.add_to(m)

folium.LayerControl().add_to(m)

# display(m)

# save map
m.save("../plots/00f_H_d2_nodes.html")
# del(m)
winsound.Beep(350, 700)


# In[49]:


# plot on folium - 2 cols + banw (to make sure the correct subgraph has been induced)
m = folium.Map(location=[55.6761, 12.5683], zoom_start = 11.5) 

H_edge_coorddict = nx.get_edge_attributes(H, "coord")
H_typedict = nx.get_edge_attributes(H, "category_edge")
H_noncar_induced_edge_coorddict = nx.get_edge_attributes(H_noncar_induced, "coord")
banw_edge_coorddict = nx.get_edge_attributes(banw, "coord")

# entire H network
mynw_fg = folium.FeatureGroup(name = "H (all edges)", show = False)
# car-only edges
cnw = folium.FeatureGroup(name = "car only", show = False)
# bike edges
bnw = folium.FeatureGroup(name = "bike only", show = False)
# multi edges
multi_nw = folium.FeatureGroup(name = "mixed (car+bike)", show = False)
# induced subgraph
induced_nw = folium.FeatureGroup(name = "induced subgraph (non-car nodes)", show = True)
# bikeable edges
bikeable_nw = folium.FeatureGroup(name = "bikeable (bikeonly+mixed)", show = True)

hloc = []
cloc = []
bloc = []
mloc = []
# bikeable
baloc = []
# induced
inloc = []

for edge in H.edges:
    myloc = [(c[1], c[0]) for c in H_edge_coorddict[tuple(sorted(edge))].coords]
    hloc.append(myloc)
    if H_typedict[edge] == "car":
        cloc.append(myloc)
    if H_typedict[edge] == "bike":
        bloc.append(myloc)
    if H_typedict[edge] == "multi":
        mloc.append(myloc)

for edge in banw.edges:
    myloc = [(c[1], c[0]) for c in banw_edge_coorddict[edge].coords]
    baloc.append(myloc)
for edge in H_noncar_induced.edges:
    myloc = [(c[1], c[0]) for c in H_noncar_induced_edge_coorddict[edge].coords]
    inloc.append(myloc)

    
mynw_line = folium.PolyLine(locations = hloc, weight = 2, color = "black").add_to(mynw_fg)
cnw_line = folium.PolyLine(locations = cloc, weight = 2, color = "grey").add_to(cnw)
bnw_line = folium.PolyLine(locations = bloc, weight = 2, color = "green").add_to(bnw)
multi_nw_line = folium.PolyLine(locations = mloc, weight = 2, color = "yellow").add_to(multi_nw)
bikeable_nw_line = folium.PolyLine(locations = baloc, weight = 2, color = "blue").add_to(bikeable_nw)
induced_nw_line = folium.PolyLine(locations = inloc, weight = 2, color = "red").add_to(induced_nw)

mynw_fg.add_to(m)
cnw.add_to(m)
bnw.add_to(m)
multi_nw.add_to(m)
induced_nw.add_to(m)
bikeable_nw.add_to(m)

folium.LayerControl().add_to(m)

#display(m)

# save map
m.save("../plots/00_H_banw.html")
# del(m)
winsound.Beep(350, 700)

