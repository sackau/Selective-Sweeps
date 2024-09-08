import sys
import pandas as pd
import numpy as np
from collections import Counter
class sfs():
    def __init__(self,data_file):
        self.data_file=data_file
        self.loaded_data1=[]
        self.count_ones=[]
        self.frequencies=None
        self.hist_scaled=None
        self.N1=None
        self.max_count_ones=None
    def load_data(self):
        mutation_data=[]
        all_lines=[]
        with open(self.data_file) as file:
            for line in file:
                all_lines.append(line.strip())
            for i, line in enumerate(all_lines): # enumerate returns index and value
                if line.startswith('positions:'):
                    start_index=i+1
                    break
            if all_lines[-1]=="1.0: freq":
                end_index=-2
            else:
                end_index=-1
            for line in all_lines[start_index:end_index]:
                mutation_data.append(pd.Series(list(line)).astype(int))
        mutation_data = pd.DataFrame(mutation_data)
        return mutation_data
    def load_data_1(self,loaded_data):
        self.loaded_data_1=loaded_data[loaded_data[0]==1]
        
        return self.loaded_data_1
    def count_ones_func(self, loaded_data_1):
        self.count_ones=loaded_data_1.sum(axis=0).tolist()
        return self.count_ones
    def create_histogram(self, count_ones):
        return dict(Counter(count_ones)) # dict function gives inside values as 1857:2 
    def normalize_histogram(self,histogram,N1):
        max_key=max(histogram.keys()) # maximum x value in histogram
        # initialise the normalized histogram list with zeros
        normalized_histogram=[0]*(2*N1+1)
        for key, count in histogram.items(): # run the for loop for key and their freq
            # normalize between 0 and 1
            normalized_key=key/max_key 
            bin_index=int(normalized_key*2*N1) # determine which bin the key will fall into
            bin_index=min(bin_index, 2*N1-1)
            normalized_histogram[bin_index]+=count*2*N1
        return normalized_histogram
    def add_histograms(self,hist1,hist2):
        result=[]
        for a, b in zip(hist1,hist2): # zip command combine the two lists element wise
            result.append(a+b) # add them element wise and append them
        return result

if __name__=="__main__":
	i = int(sys.argv[1]) 
	N1=1000
	total_files=0
	total_histogram=[0]*(2*N1+1) # initialize
	for j in range(1,1000):
		file_path=f"/work/users/s/a/sackau/UNC_work/N_1000_s_0.01_h_0.5/Post_fixation/h_0.5/txt_files/{i}/Ns_100_output_t_100_subs_off_{(i-1)*1000+j}.txt"
		try:
			with open(file_path,'r'):
				obj=sfs(file_path)
				loaded_data=obj.load_data() # get the data
				loaded_data_1=obj.load_data_1(loaded_data) # filter subpop from the data
				count_ones=obj.count_ones_func(loaded_data_1) # count # of ones at each locus
				histogram=obj.create_histogram(count_ones) # make the histogram (SFS) of the above
				normalized_histogram=obj.normalize_histogram(histogram,N1) # do normalization of SFS
				for k in range(2*N1+1):
					total_histogram[k]+=normalized_histogram[k]
				total_files+=1
				averaged_histogram=[x/total_files for x in total_histogram]
				x_values = [i / (2*N1) for i in range(2*N1 + 1)]
				output_folder=f"/work/users/s/a/sackau/UNC_work/N_1000_s_0.01_h_0.5/Post_fixation/h_0.5/Ns_100_average_data"
				with open(f"{output_folder}/Ns_10_averaged_data_{i}.txt", "w") as outfile:
					for x, y in zip(x_values, averaged_histogram):
						outfile.write(f"{x}\t{y}\n")
		except FileNotFoundError:
				print(f" file not found : {file_path}. Skipping" )