#!/usr/bin/env python


import os
import gzip
import re
import pickle
import hashlib

import toytree
import numpy as np


def digest(obj, algo="md5"):
    """
    Meant to emulate R's digest::digest(obj, algo).

    algo: 'md5', 'sha1', 'sha256', etc.
    """
    data = pickle.dumps(obj, protocol=pickle.HIGHEST_PROTOCOL)
    h = hashlib.new(algo)
    h.update(data)
    return h.hexdigest()

def load_split_heights(folders_list, sampnums = [5000]):
    """
    Load split heights for all species pairs
    from ARG inference
    """
    out_dict = {'cd': [], 'cr': [], 'dr': []}
    for folder in folders_list:
        files = [folder+'/out.{}.smc.gz'.format(num) for num in sampnums]
        m = re.search('\.\d+_\d+_\d+\.', os.path.basename(folder))
        region_str = m.group()[1:-1]
        region_start = int(region_str.split('_')[1])
        region_stop = int(region_str.split('_')[2])
        region = [region_start,region_stop]
        positions = range(region[0],region[1],2500)
        len(positions)
        
        post_cd_heights = []
        post_cr_heights = []
        post_dr_heights = []
        
        for pos in positions:
            post_cd_heights_pos = []
            post_cr_heights_pos = []
            post_dr_heights_pos = []
            for file in files:
                f=gzip.open(file,'rt')
                names=f.readline().split()[1:]
                names_dict = dict(zip(range(len(names)),names))
                
                cusp_keys = np.array(list(range(len(names))))[[i.startswith('pcusp') for i in names_dict.values()]]
                drum_keys = np.array(list(range(len(names))))[[i.startswith('pdrum') for i in names_dict.values()]]
                roem_keys = np.array(list(range(len(names))))[[i.startswith('proe') for i in names_dict.values()]]
                
                f=gzip.open(file,'rt')
                for line in f.readlines():
                    if line.startswith('TREE'):
                        linesplit = line.split()
                        if (pos>=int(linesplit[1])) & (pos<=int(linesplit[2])):
                            break
                newick = linesplit[3]
                
                t = toytree.tree(newick)
                
                cd_heights = []
                cr_heights = []
                dr_heights = []

                for ck in cusp_keys:
                    for dk in drum_keys:
                        cd_heights.append(t.mod.prune(str(ck),str(dk)).treenode.height)
                for ck in cusp_keys:
                    for rk in roem_keys:
                        cr_heights.append(t.mod.prune(str(ck),str(rk)).treenode.height)
                for dk in drum_keys:
                    for rk in roem_keys:
                        dr_heights.append(t.mod.prune(str(dk),str(rk)).treenode.height)
                
                post_cd_heights_pos.append(cd_heights)
                post_cr_heights_pos.append(cr_heights)
                post_dr_heights_pos.append(dr_heights)
                
            post_cd_heights.append(post_cd_heights_pos)
            post_cr_heights.append(post_cr_heights_pos)
            post_dr_heights.append(post_dr_heights_pos)
            
        out_dict['cd'].append(np.array(post_cd_heights))
        out_dict['cr'].append(np.array(post_cr_heights))
        out_dict['dr'].append(np.array(post_dr_heights))

    return out_dict
