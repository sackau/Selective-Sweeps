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
        position_index=None
        threshold=0.092
        # my .txt file will have data between position line and last line is Fixed
        with open(self.data_file) as file:
            for line in file:
                all_lines.append(line.strip())
            for i, line in enumerate(all_lines): # enumerate returns index and value
                if line.startswith('positions:'):
                    positions = [float(pos) for pos in line.split()[1:]]
                    for idx, value in enumerate(positions):
                          if value > threshold:
                                position_index=idx
                                break
                    start_index=i+1
                    break
            end_index=-2
            for line in all_lines[start_index:end_index]:
                mutation_data.append(pd.Series(list(line)).astype(int))
        mutation_data = pd.DataFrame(mutation_data)
        mutation_data=mutation_data.iloc[:,position_index:]  # here filtered the position after which i want data
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
	for j in range(1,250):
		file_path=f"/work/users/s/a/sackau/UNC_work/time_dependent/Ns_100/Ns_100_CI/Parul_simulations/faster_simulations/no_CI/dist_files/{i}/output_postfixation_{j}.ms"
		try:
			with open(file_path,'r'):
				obj=sfs(file_path)
				loaded_data=obj.load_data() # get the data
				#loaded_data_1=obj.load_data_1(loaded_data) # filter subpop from the data
				count_ones=obj.count_ones_func(loaded_data) # count # of ones at each locus
				histogram=obj.create_histogram(count_ones) # make the histogram (SFS) of the above
				normalized_histogram=obj.normalize_histogram(histogram,N1) # do normalization of SFS
				for k in range(2*N1+1):
					total_histogram[k]+=normalized_histogram[k]
				total_files+=1
				averaged_histogram=[x/total_files for x in total_histogram]
				x_values = [i / (2*N1) for i in range(2*N1 + 1)]
				output_folder=f"/work/users/s/a/sackau/UNC_work/time_dependent/Ns_100/Ns_100_CI/Parul_simulations/faster_simulations/no_CI/average_data_t_100"
				with open(f"{output_folder}/averaged_data_t_100_{i}.txt", "w") as outfile:
					for x, y in zip(x_values, averaged_histogram):
						outfile.write(f"{x}\t{y}\n")
		except FileNotFoundError:
				print(f" file not found : {file_path}. Skipping" )
