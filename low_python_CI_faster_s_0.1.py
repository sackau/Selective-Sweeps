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
        with open(self.data_file) as file:
            for line in file:
                all_lines.append(line.strip())
            for i, line in enumerate(all_lines): # enumerate returns index and value
                if line.startswith('Mutations:'):
                    start_index=i+1
                    break
            for i, line in enumerate(all_lines): # enumerate returns index and value
                if line.startswith('Individuals:'):
                    end_index=i
                    break
            for line in all_lines[start_index:end_index]:
                # Split line into components
                parts = line.split()
                # Convert numeric values where applicable, keep others as strings
                processed_parts = [int(part) if part.isdigit() else float(part) if part.replace('.', '', 1).isdigit() else part for part in parts]
                mutation_data.append(processed_parts) # Append processed parts to data list
        mutation_data = pd.DataFrame(mutation_data, columns=['ID', 'ID2', 'm', 'pos', 'Zero', 'h', 'Population', 'time_arise', 'freq'])
        return mutation_data
    def create_histogram(self, count_ones):
        return dict(Counter(count_ones)) # dict function gives inside values as 1857:2 
    def normalize_histogram(self,histogram,N1):
        #max_key=max(histogram.keys()) # maximum x value in histogram
        max_key=2*N1 
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
        file_path=f"/work/users/s/a/sackau/UNC_work/time_dependent/Ns_100/Ns_100_CI/Parul_simulations/faster_simulations/low_CI_paper/s_0.1/dist_files/{i}/output_t_100_subs_off_{(i-1)*1000+j}.txt"
        try:
            with open(file_path,'r'):
                obj=sfs(file_path)
                loaded_data=obj.load_data() # get the data
                loaded_data_2=loaded_data[loaded_data['m']=='m1'] # filtered the m1s (neutrals) from the output
                count_ones=loaded_data_2['freq'].tolist() # count # of ones at each locus
                histogram=obj.create_histogram(count_ones) # make the histogram (SFS) of the above
                normalized_histogram=obj.normalize_histogram(histogram,N1) # do normalization of SFS
                for k in range(2*N1+1):
                    total_histogram[k]+=normalized_histogram[k]
                total_files+=1
                averaged_histogram=[x/total_files for x in total_histogram]
                x_values = [i / (2*N1) for i in range(2*N1 + 1)]
                output_folder=f"/work/users/s/a/sackau/UNC_work/time_dependent/Ns_100/Ns_100_CI/Parul_simulations/faster_simulations/low_CI_paper/s_0.1/average_data"
                with open(f"{output_folder}/averaged_data_{i}.txt", "w") as outfile:
                    for x, y in zip(x_values, averaged_histogram):
                        outfile.write(f"{x}\t{y}\n")
        except FileNotFoundError:
            print(f" file not found : {file_path}. Skipping" )