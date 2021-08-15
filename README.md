# msc-bikegaps readme

import parameters and custom functions: run ppf.py

step 1 - make simplified network: 
input = raw OSM data (car and bicycle nodes and links, 1 csv file each)
output = simplified network *H* (nx) resp. *h* (ig), subgraph *B* (nx) resp. *b* (ig)

step 2 - get edge betweenness centrality values: 
input = *H* and *h* from step 1
output = *ebc* df of edge betweenness centralities

step 3 - get list of gaps:
input = *H* and *h* from step 1, *ebc* from step 2
output = *mygaps* list of gaps ranked by mbar_c_gap
