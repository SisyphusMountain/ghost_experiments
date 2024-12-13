#!/usr/bin/env python3
# by Theo Tricou

from ete3 import Tree as tr
import sys, os
import glob
import re
import pandas as pd




if len(sys.argv) < 3:
    ale_version = "aletree_"
else:
    ale_version = sys.argv[2] + "_"

ale = tr("aletree", format = 1)
ext = tr("SampledSpeciesTree.nwk", format = 1)

map = dict()

map[ale.name] = ext.name

# postorder traversal of the ale tree to find the appropriate mapping
for node in ale.traverse("postorder"):
    # if the node is a leaf, we will know the right mapping between them. Otherwise, its left child will already have been visited.
    # in this case, to find the correpondance, we find the left child in the ale tree, find the corresponding node in the ext tree, go up in the ext 
    # tree, and the name of this node in the ext tree will be the mapping of the node in the ale tree.
    if node.is_leaf():
        map[node.name] = node.name
    else:
        # find the left child of this node
        left_child = node.children[0]
        # find the corresponding node in the ext tree
        corresponding_node = ext.search_nodes(name = map[left_child.name])[0]
        # go up in the ext tree
        corresponding_node = corresponding_node.up
        # the name of this node in the ext tree is the mapping of the node in the ale tree
        map[node.name] = corresponding_node.name




all_uTs = []
uTs = glob.glob(ale_version + "*sampledtree.nwk.ale.uTs")
for i in uTs:
    if os.stat(i).st_size != 0:
        try:
            all_uTs.append(pd.read_csv(i, delim_whitespace=True, header = None))
        except:
            print("Currently in folder: ", os.getcwd())
            raise ValueError("Error in file: ", i)

full_uTs = pd.concat(all_uTs, ignore_index=True)
full_uTs[[0, 1]] = full_uTs[[0, 1]].astype(str)


def is_to_strict_descendant(name_from ,name_to, tree):
    n_from = tree.search_nodes(name = name_from)[0]
    n_to = tree.search_nodes(name = name_to)[0]
    if n_to in n_from.get_descendants():
        return(True)
    else:
        return(False)

def is_to_time_descendant(name_from ,name_to, tree):
    n_from = tree.search_nodes(name = name_from)[0]
    n_to = tree.search_nodes(name = name_to)[0].up
    if n_from.get_distance(tree) >= n_to.get_distance(tree):
        return(False)
    elif n_from.get_distance(tree) <= n_to.get_distance(tree):
        return(True)
    else:
        return("ERROR")


res_df = pd.DataFrame()
res_df["ALE"] = [i.name for i in ale.traverse()]
try:
    res_df["Node"] = [map[i.name] for i in ale.traverse()]
except:
    raise ValueError(f"Error in mapping. {map}")
res_df["br_length"] = res_df.apply(lambda x: ale.search_nodes(name = x["ALE"])[0].dist, axis = 1)
#res_df["dist_to_root"] = res_df.apply(lambda x: round(ale.search_nodes(name = x["ALE"])[0].get_distance(ale),6), axis = 1)
res_df["N_transfers_donor"] = res_df.apply(lambda x: round(full_uTs.loc[full_uTs[0] == x["ALE"]][2].sum(), 6), axis = 1)
#res_df["N_transfers_recip"] = res_df.apply(lambda x: round(full_uTs.loc[full_uTs[1] == x["ALE"]][2].sum(), 6), axis = 1)
full_uTs["to_direct_descendant"] = full_uTs.apply(lambda x: is_to_strict_descendant(x[0], x[1], ale), axis = 1)
full_uTs["to_descendant"] = full_uTs.apply(lambda x: is_to_time_descendant(x[0], x[1], ale), axis = 1)
#res_df["to_direct_descendant"] = res_df.apply(lambda x: round(full_uTs.loc[(full_uTs[0] == x["ALE"]) & (full_uTs.to_direct_descendant)][2].sum(), 6), axis = 1)
res_df["to_descendant"] = res_df.apply(lambda x: round(full_uTs.loc[(full_uTs[0] == x["ALE"]) & (full_uTs.to_descendant)][2].sum(), 6), axis = 1)


def get_ghost_prediction(com_node, back_bone_nodes, back_bone_tree):
    # com_node = com.search_nodes(name = node.name)[0]
    while com_node.name not in back_bone_nodes:
        com_node = com_node.up
    ext_node = back_bone_tree.search_nodes(name = com_node.name)[0]
    while len(ext_node.get_children()) not in [0,2]:
        ext_node = ext_node.get_children()[0]
    return(ext_node.name)

def get_ghost_banch_length(node, nodes_contemporary, com):
    if node.name in nodes_contemporary:
        return(node.dist)
    else:
        dist_to_root = node.get_distance(com)
        distup_to_root = node.up.get_distance(com)
        if dist_to_root > extroot_to_comroot and distup_to_root > extroot_to_comroot:
            ext_bl = round(node.dist, 6)
        elif dist_to_root > extroot_to_comroot and distup_to_root < extroot_to_comroot:
            ext_bl = round((dist_to_root - distup_to_root) - (extroot_to_comroot - distup_to_root), 6)
        else:
            ext_bl = 0
        return(ext_bl)



com = tr("../T/CompleteTree.nwk",format = 1)
samp = tr("SampledSpeciesTree.nwk", format = 1)
back_bone_nodes = []
for i in samp:
    node = com.search_nodes(name = i.name)[0]
    while node:
        back_bone_nodes.append(node.name)
        node = node.up


back_bone_nodes = list(set(back_bone_nodes))
back_bone_tree = com.copy()
back_bone_tree.prune(list(set(back_bone_nodes)))


anc = com.search_nodes(name=samp.name)[0]
extroot_to_comroot = anc.get_distance(com)
nodes_contemporary = [node.name for node in anc.iter_descendants()]
stat_by_node = dict()
for i in com.traverse():
    prediction = get_ghost_prediction(i, back_bone_nodes, back_bone_tree)
    if not i.is_root():
        i.add_features(ghost_prediction = prediction, ghost_dist = get_ghost_banch_length(i, nodes_contemporary, com))
    else:
        i.add_features(ghost_prediction = prediction, ghost_dist = 0)
    if prediction in stat_by_node:
        if i.name not in back_bone_tree:
            stat_by_node[prediction][0] += 1
            stat_by_node[prediction][1] += i.ghost_dist
    else:
        stat_by_node[prediction] = [0, 0]

#res_df["N_ghost_node"] = res_df.apply(lambda x: stat_by_node[x["Node"]][0], axis = 1)
res_df["L_ghost_branch"] = res_df.apply(lambda x: round(stat_by_node[x["Node"]][1],6), axis = 1)





res_df.to_csv(sys.argv[1], sep=' ', index = False)

# GNU Ghost