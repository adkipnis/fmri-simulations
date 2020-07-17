#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Imagenet image label utilities:
 1) Use a prebuilt dictionary to map synset IDs to image labels 
 2) For any synset ID: get image label by querying imagenet.stanford.edu
 3) Generate tar shell commands to extract imagenet synset subdirectories

@author: alex
"""

class SynsetIDContainer():   
    os = __import__('os')
    np = __import__('numpy')
    pd = __import__('pandas')
    
    def __init__(self, n_stim = None):       
        self.n_stim = n_stim
        
    def import_synset_list(self, synset_list_dir):
        self.synset_list_dir = synset_list_dir
        self.synset_list = self.pd.read_csv(self.synset_list_dir, sep=" ")
        self.synset_column = self.synset_list.columns[0]
        self.synset_list = self.synset_list[self.synset_column].tolist()
    
    def import_events_tsv(self, events_dir, synset_colname = "stim_id", decimals = 6):
        self.events_dir = events_dir
        self.events = self.pd.read_csv(self.events_dir, sep='\t')
        self.stimulus_ids = self.np.unique(self.np.around(self.events[synset_colname].tolist(), decimals))
        self.synset_list = self.stimulus_ids[~self.np.isnan(self.stimulus_ids)].astype(str)
            
    def rename_synsets(self, remove_zeros = True):
        if self.n_stim is None:
            self.n_stim = len(self.synset_list)
        self.synset_ids = []
        self.image_numbers = []
        for i in range(self.n_stim):
            self.synset_ids.append('n0'+self.synset_list[:self.n_stim][i].split('.')[0])
            self.image_numbers.append(self.synset_list[:self.n_stim][i].split('.')[1])
        if remove_zeros:
            self.image_numbers = self.np.array(self.image_numbers).astype(int).astype(str).tolist()
            
### official map from imagenet
    def build_from_official_map(self):
        import urllib
        self.online_map_url = "http://image-net.org/archive/words.txt"
        self.online_map_html = urllib.request.urlopen(self.online_map_url)
        all_synset_ids = []
        all_label_ids = []
        for line in self.online_map_html:
            all_synset_ids.append(line.decode("utf-8").split('\t')[0])
            all_label_ids.append(line.decode("utf-8").split('\t')[1].split('\n')[0])
        self.official_map = dict(zip(all_synset_ids, all_label_ids))     
                
    def import_official_map(self, official_map_dir):      
        self.official_map_dir = official_map_dir
        self.official_map_df = self.pd.read_csv(self.official_map_dir, sep=" ", header=None)
        self.official_map_df.columns = ["synset_id", "label"]
        self.official_map = self.official_map_df.set_index('synset_id').T.to_dict('list')

    def use_official_map(self):
        self.sub_dict = {}
        self.missing_labels = []
        for i in range(self.n_stim):
            if self.synset_ids[i] in self.official_map.keys():
                self.sub_dict.update({self.synset_ids[i]: self.official_map[self.synset_ids[i]]})
            else:
                self.missing_labels.append(self.synset_ids[i])
        if len(self.missing_labels) > 0:
            print("Warning: official map does not contain all synset IDs of interest. See get_missing_labels() method for more.")

        
### map_clsloc       
    def import_map_clsloc(self, map_clsloc_dir):
        self.map_clsloc_dir = map_clsloc_dir
        self.map_clsloc_df = self.pd.read_csv(self.map_clsloc_dir, sep=" ", header=None)
        self.map_clsloc_df.columns = ["synset_id", "category_id", "label"]
        self.map_clsloc_df = self.map_clsloc_df.drop(columns=["category_id"])
        self.map_clsloc = self.map_clsloc_df.set_index('synset_id').T.to_dict('list')

    def use_map_clsloc(self):
        self.sub_dict = {}
        self.missing_labels = []
        for i in range(self.n_stim):
            if self.synset_ids[i] in self.map_clsloc.keys():
                self.sub_dict.update({self.synset_ids[i]: self.map_clsloc[self.synset_ids[i]][0]})
            else:
                self.missing_labels.append(self.synset_ids[i])
        if len(self.missing_labels) > 0:
            print("Warning: map_clsloc does not contain all synset IDs of interest. See get_missing_labels() method for more.")

### Web scraping (not recommended)
    def query_labels(self, base_url=None, html_element = 'detailsPanel'):
        self.base_url = base_url
        if self.base_url is None:
            self.base_url = 'http://imagenet.stanford.edu/synset?wnid='
            
        from selenium import webdriver
        self.sub_dict = {}
        for synset in self.synset_ids:  
            driver = webdriver.Firefox()
            driver.implicitly_wait(30)
            url = self.base_url + synset
            driver.get(url)
            element = driver.find_element_by_id(html_element) 
            panel_text = element.text
            label = panel_text.partition('\n')[0]
            self.sub_dict.update({synset:label})   
            driver.quit()

    def query_missing_labels(self, base_url=None, html_element = 'detailsPanel'):
        self.base_url = base_url
        if self.base_url is None:
            self.base_url = 'http://imagenet.stanford.edu/synset?wnid='
            
        from selenium import webdriver
        for synset in self.missing_labels:  
            driver = webdriver.Firefox()
            driver.implicitly_wait(30)
            url = self.base_url + synset
            driver.get(url)
            element = driver.find_element_by_id(html_element) 
            panel_text = element.text
            label = panel_text.partition('\n')[0]
            self.sub_dict.update({synset:label})   
            driver.quit()

    def query_image_urls(self):
        import urllib
        self.image_url_base = "http://www.image-net.org/api/text/imagenet.synset.geturls?wnid="
        self.url_dict = {}        
        for synset in self.synset_ids:  
            collection_url = self.image_url_base + synset
            collection_html = urllib.request.urlopen(collection_url)
            urllist = []
            for line in collection_html:               
                urllist.append(line.decode("utf-8").split('\n')[0])
            self.url_dict.update({synset : urllist})
                
    def query_image_urls_sparse(self):
        import urllib
        self.image_url_base = "http://www.image-net.org/api/text/imagenet.synset.geturls?wnid="
        self.url_dict = {}        
        for i in range(len(self.synset_ids)):            
            collection_url = self.image_url_base + self.synset_ids[i]
            collection_html = urllib.request.urlopen(collection_url)
            line_num = int(self.image_numbers[i])
            n = 0            
            
            for line in collection_html:               
                n+=1
                if n == line_num:
                    self.url_dict.update({self.synset_ids[i]: line.decode("utf-8").split('\n')[0]})   
            
                      
### Crop labels
    def crop_labels(self):
        self.full_subdict = self.sub_dict
        cropped_dict = {}
        for synset in self.sub_dict.keys():
            cropped_dict.update({synset : self.sub_dict[synset].split(',')[0]})
        self.sub_dict = cropped_dict       

### Shell commands
    def compile_untar_commands(self, tar_name = None):
        self.tar_name = tar_name
        if self.tar_name is None:
            self.tar_name = "fall11_whole"
            
        self.tar_commands = []
        for i in range(self.n_stim):
            self.tar_commands.append("tar -xvf "+ self.tar_name + ".tar " + self.synset_ids[i] +".tar ")
   
    def relative_jpeg_filenames(self, drop_first = 0):
        self.rel_jpeg = []
        for i in range(len(self.synset_ids)):
            self.rel_jpeg.append(self.os.path.join(self.synset_ids[i], self.synset_ids[i]+ "_" + self.image_numbers[i][drop_first:] + ".JPEG"))
     
    def absolute_jpeg_filenames(self, abs_directory, drop_first = 0):
        self.abs_jpeg = []
        for i in range(len(self.synset_ids)):
            self.abs_jpeg.append(self.os.path.join(abs_directory, self.synset_ids[i], self.synset_ids[i]+ "_" + self.image_numbers[i][drop_first:] + ".JPEG"))

    def compile_cp_commands(self, cp_command = 'cp', goal_dir = None):
            self.cp_command = cp_command
            self.goald_dir = goal_dir
            if self.goald_dir is None:
                self.goald_dir = self.os.getcwd()
            self.cp_commands = []
            for abs_path in self.abs_jpeg:
                self.cp_commands.append(self.cp_command + abs_path + ' ' + self.goald_dir)
                            
### Saver methods
    def save_new_dict(self, filename = None):   
        self.filename = filename
        if self.filename is None:
            self.filename = self.os.path.join(self.os.getcwd(), 'custom_synset_dictionary')                  
        self.np.save(self.filename, self.get_new_dict())
        print("Saved new dictionary to:", self.filename + ".npy")
    
    def save_untar_commands(self, filename = None): 
        self.filename = filename
        if self.filename is None:
            self.filename = self.os.path.join(self.os.getcwd(), 'untar_commands')     
        with open(self.filename + ".txt", "w") as outfile:
            outfile.write("\n".join(self.tar_commands))
        print("Saved untar commands to:", self.filename + ".txt")
     
    def save_jpeg_filenames(self, filename = None, path_type = 'rel'): 
        self.filename = filename
        if self.filename is None:
            self.filename = self.os.path.join(self.os.getcwd(), 'jpeg_filenames')            
        if path_type == 'rel':
            with open(self.filename + ".txt", "w") as outfile:
                outfile.write("\n".join(self.get_rel_jpeg_paths()))
        elif path_type == 'abs':
            with open(self.filename + ".txt", "w") as outfile:
                outfile.write("\n".join(self.get_abs_jpeg_paths()))
        print("Saved jpeg paths to:", self.filename + ".txt")
    
    def save_cp_commands(self, filename = None): 
        self.filename = filename
        if self.filename is None:
            self.filename = self.os.path.join(self.os.getcwd(), 'cp_commands')            
        with open(self.filename + ".txt", "w") as outfile:
            outfile.write("\n".join(self.get_cp_commands()))
        print("Saved cp commands to:", self.filename + ".txt")    
        
### Getter and setter methods    
    def set_new_dict(self, sub_dict):
        assert type(sub_dict) == dict, "imported object must be of type dict"    
        self.sub_dict = sub_dict    
    def set_official_map(self, official_map):
        assert type(official_map) == dict, "imported object must be of type dict"    
        self.official_map = official_map 
    def get_synset_ids(self):
        return self.synset_ids.copy()
    def get_image_numbers(self):
        return self.image_numbers.copy()
    def get_map_clsloc(self):
        return self.map_clsloc_dict.copy()
    def get_official_map(self):
        return self.official_map.copy()
    def get_url_dict(self):
        return self.url_dict.copy()
    def get_untar_commands(self):
        return self.tar_commands.copy()   
    def get_new_dict(self):
        return self.sub_dict.copy()
    def get_uncropped_dict(self):
        return self.full_subdict.copy()    
    def get_missing_labels(self):
        return self.missing_labels.copy()
    def get_rel_jpeg_paths(self):
        return self.rel_jpeg.copy()
    def get_abs_jpeg_paths(self):
        return self.abs_jpeg.copy()
    def get_cp_commands(self):
        return self.cp_commands.copy()
    
###############################################################################

import os

# Directories and Intitialization
ds_dir = "/home/alex/Datasets/ds001246/"
map_clsloc_dir = os.path.join(ds_dir, "map_clsloc.txt")
beta_ids_dir = os.path.join(ds_dir, "derivatives/SPM/sub-01/perceptionTest01/run-01/sub-01-perceptionTest01-run-1_stim_ids.txt")
events_dir = os.path.join(ds_dir, "sub-01/ses-perceptionTest01/func/sub-01_ses-perceptionTest01_task-perception_run-01_events.tsv")


# Intilialize synset ID container
# synsets = SynsetIDContainer(n_stim)
# synsets.import_synset_list(beta_ids_dir)
synsets = SynsetIDContainer()
synsets.import_events_tsv(events_dir)
synsets.rename_synsets()
synsets.get_synset_ids()
synsets.get_image_numbers()


# # Try with map_clsloc from https://gist.github.com/aaronpolhamus/964a4411c0906315deb9f4a3723aac57
# synsets.import_map_clsloc(map_clsloc_dir)
# synsets.use_map_clsloc()
# synsets.get_new_dict()
# synsets.get_missing_labels()

# # The previous dictionary is incomplete - get remaining labels by webscraping (not recommended, slow)
# synsets.query_missing_labels()

# Try using the official map by imagenet
if 'official_map' in locals():
    synsets.set_official_map(official_map)
else:
    synsets.build_from_official_map()
    official_map = synsets.get_official_map()
synsets.use_official_map()

# Get full url dictionary for all synset IDs
# synsets.query_image_urls()
# url_dict = synsets.get_url_dict()

# Create untar commands and paths to images
synsets.compile_untar_commands()
# synsets.relative_jpeg_filenames()
synsets.absolute_jpeg_filenames("/media/heiko/Disk/imagenet/ILSVRC2012_img_val")
jpg_paths = synsets.get_abs_jpeg_paths()
synsets.compile_cp_commands(cp_command = "scp", goal_dir = None)
synsets.get_cp_commands()
# synsets.save_cp_commands()
    
# Crop new labels
synsets.get_new_dict()
synsets.crop_labels()

# # You can still access the original labels and reset your new_dict
# uncropped_dict = synsets.get_uncropped_dict()
# synsets.set_new_dict(uncropped_dict)

# Save results
synsets.save_untar_commands(ds_dir + "untar_commands")
synsets.save_jpeg_filenames(ds_dir + "jpeg_filenames")
synsets.save_new_dict(ds_dir + "custom_synset_dictionary")
