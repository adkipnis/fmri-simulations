#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Get all image labels from cryptic synset IDs (mostly by querying imagenet.stanford.edu)
Generate tar shell commands to extract imagenet subdirectories
"""

import pandas as pd
from selenium import webdriver

# Directories and Intitialization
map_clsloc_dir = "/home/alex/ds001246/map_clsloc.txt"
beta_ids_dir = "/home/alex/ds001246/derivatives/SPM/sub-01/perceptionTest01/run-01/" \
    + "sub-01-perceptionTest01-run-1_stim_ids.txt"
n_stim = 50
class_id = []
find_label = []
tar_commands = []

# load the goddamn image ids
map_clsloc = pd.read_csv(map_clsloc_dir, sep=" ", header=None)
map_clsloc.columns = ["class_id", "category_id", "label"]
map_clsloc = map_clsloc.drop(columns=["category_id"])
map_clsloc_dict = map_clsloc.set_index('class_id').T.to_dict('list')
beta_ids = pd.read_csv(beta_ids_dir, sep=" ")
label_dict = {}

# extract corresponding class ids
for i in range(n_stim):
    class_id.append('n0'+beta_ids['RegressorNames'][:n_stim][i].split('.')[0])
    if class_id[i] in map_clsloc_dict.keys():
        label_dict.update({class_id[i]: map_clsloc_dict[class_id[i]][0]})
    else:
        find_label.append(class_id[i])
    tar_commands.append("tar -xvf fall11_whole.tar " + class_id[i] +".tar ")

# save tar commands to .txt
with open("tar_commands.txt", "w") as outfile:
    outfile.write("\n".join(tar_commands))

# query remaining synsets with selenium
for synset in find_label:  
    driver = webdriver.Firefox()
    driver.implicitly_wait(30)
    url = 'http://imagenet.stanford.edu/synset?wnid=' + synset
    driver.get(url)
    element = driver.find_element_by_id('detailsPanel') 
    panel_text = element.text
    label = panel_text.partition('\n')[0]
    label_dict.update({synset:label})   
    driver.quit()
    
# for synset in label_dict.keys():  
#     driver = webdriver.Firefox()
#     driver.implicitly_wait(30)
#     url = 'http://imagenet.stanford.edu/synset?wnid=' + synset
#     driver.get(url)
#     element = driver.find_element_by_id('detailsPanel') 
#     panel_text = element.text
#     label = panel_text.partition('\n')[0]
#     label_dict.update({synset:label})   
#     driver.quit()

# save label_dictionary 
import numpy as np
np.save("/home/alex/ds001246/stimulus_label_dictionary", label_dict)

# # Overwrite existing labels
# label_dict_1= np.load("/home/alex/ds001246/stimulus_label_dictionary_1.npy",allow_pickle='TRUE').item()
# label_dict_2= np.load("/home/alex/ds001246/stimulus_label_dictionary_2.npy",allow_pickle='TRUE').item()

# for synset in label_dict_2.keys():
#     label_dict_1[synset] = label_dict_2[synset]


