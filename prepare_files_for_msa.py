#importing libraries

import pandas as pd

#reading csv file input

FILENAME_CSV = ["/Users/User/Documents/ihf_hu_project/ihf_hu_project/data/HUA_1030198_bacterial_protein_id_nums_only.csv",
                "/Users/User/Documents/ihf_hu_project/ihf_hu_project/data/HUB_1030195_bacterial_protein_id_nums_only.csv",
                "/Users/User/Documents/ihf_hu_project/ihf_hu_project/data/IHFA_1030196_bacterial_protein_id_nums_only.csv",
                "/Users/User/Documents/ihf_hu_project/ihf_hu_project/data/IHFB_1030197_bacterial_protein_id_nums_only.csv"]

PREFIX = ["HU_alpha","HU_beta","IHF_alpha","IHF_beta"]

FASTA = ["/Users/User/Documents/ihf_hu_project/ihf_hu_project/data/HUA_1030198_all_390_members_including_bacteria.fasta",
         "/Users/User/Documents/ihf_hu_project/ihf_hu_project/data/HUB_1030195_all_993_members_including_bacteria.fasta",
         "/Users/User/Documents/ihf_hu_project/ihf_hu_project/data/IHFA_1030196_all_766_members_including_bacteria.fasta",
         "/Users/User/Documents/ihf_hu_project/ihf_hu_project/data/IHFB_1030197_all_607_members_including_bacteria.fasta"]

WRITE=["/Users/User/Documents/ihf_hu_project/ihf_hu_project/results/HUA_1030198_bacteria_nr_species.fasta",
       "/Users/User/Documents/ihf_hu_project/ihf_hu_project/results/HUB_1030195_bacteria_nr_species.fasta",
       "/Users/User/Documents/ihf_hu_project/ihf_hu_project/results/IHFA_1030196_bacteria_nr_species.fasta",
       "/Users/User/Documents/ihf_hu_project/ihf_hu_project/results/IHFB_1030197_bacteria_nr_species.fasta"]

MERGED_FASTA="/Users/User/Documents/ihf_hu_project/ihf_hu_project/results/hua_hub_ihfa_ihfb_merged_bacteria_nr_species.fasta"

__prog_name__ = 'prepare_files_for_msa.py'
__prog_desc__ = 'Takes a csv, fasta and output file name as input and filters non-redundant species level protein sequences in a fasta format.'

__author__ = 'Mohak Sharda'
__copyright__ = 'Copyright 2022'
__credits__ = ['Mohak Sharda']
__license__ = 'GPL3'
__version__ = '0.0.1'
__maintainer__ = 'Mohak Sharda'
__email__ = 'mohaks@ncbs.res.in'
__status__ = 'Development'

class PrepareForMsa(object):

    def __init__(self):
        "initialisation"
        pass

    def main(self):
        for filename_csv,prefix,fasta,write in zip(FILENAME_CSV,PREFIX,FASTA,WRITE):
            protid_list = self.make_protid_list(filename_csv)
            fasta_dict = self.read_fasta(fasta)
            filtered_bacterial_fasta_dictionary = self.filter_bacteria_prot_sequences(protid_list,fasta_dict)
            filtered_bacterial_fasta_dataframe = self.remove_species_level_redundancy(filtered_bacterial_fasta_dictionary,prefix)
            self.write_filtered_fasta_file(filtered_bacterial_fasta_dataframe,write)
        self.write_merged_filtered_fasta_file(WRITE,MERGED_FASTA)
    
    def make_protid_list(self,filename_csv):
        csv_bacteria = pd.read_csv(filename_csv)
        protid_list = csv_bacteria['Protein ID'].tolist()
        
        return protid_list

    #reading fasta file input
    def read_fasta(self,fasta):
        fasta_dict={}
        with open(fasta) as fh:
            seq=""
            for line in fh:
                if line.startswith(">"):
                    organism_name = line.rstrip().split("|")[3].split("[")[1].split("]")[0]
                    protein_id = line.rstrip().split("|")[0].split(">")[1].rstrip()
                elif not line.startswith(">") and not line.startswith("\n"):
                    seq = seq + line.rstrip()
                elif line.startswith("\n"):
                    fasta_dict[protein_id]=[seq,organism_name]
                    seq=""
        return fasta_dict
    
    def filter_bacteria_prot_sequences(self,taxon_list,fasta_dict):   
        #filter fasta sequences belonging to only bacteria
        filtered_bacterial_fasta_dictionary={}

        for protid_list in taxon_list:
            if protid_list in fasta_dict:
                bacteria_name = fasta_dict[protid_list][1]
                filtered_bacterial_fasta_dictionary[bacteria_name]=fasta_dict[protid_list][0]
        return filtered_bacterial_fasta_dictionary
    
    def remove_species_level_redundancy(self,filtered_bacterial_fasta_dictionary,prefix):
        #pd.DataFrame.from_dict(filtered_bacterial_fasta_dictionary)
        filtered_bacterial_fasta_dataframe = pd.DataFrame.from_dict(filtered_bacterial_fasta_dictionary, orient = 'Index', columns=['Sequences'])
        filtered_bacterial_fasta_dataframe.index.name = 'Bacteria_name'
        filtered_bacterial_fasta_dataframe.reset_index(inplace=True)
        filtered_bacterial_fasta_dataframe[['Genus','species','strain_info']] = filtered_bacterial_fasta_dataframe["Bacteria_name"].str.split(" ", 2, expand=True)
        filtered_bacterial_fasta_dataframe['Genus_species'] = filtered_bacterial_fasta_dataframe['Genus']+"_"+filtered_bacterial_fasta_dataframe['species']

        #remove duplicate rows
        filtered_bacterial_fasta_dataframe = filtered_bacterial_fasta_dataframe.drop_duplicates(subset='Genus_species', keep="first")
        filtered_bacterial_fasta_dataframe['Genus_species'] = prefix + "_" + filtered_bacterial_fasta_dataframe['Genus_species']

        return filtered_bacterial_fasta_dataframe

    def write_filtered_fasta_file(self,filtered_bacterial_fasta_dataframe,write):
        with open(write,"w") as fh:
            for ind in filtered_bacterial_fasta_dataframe.index:
                fh.write(">%s\n%s\n"%(filtered_bacterial_fasta_dataframe['Genus_species'][ind], filtered_bacterial_fasta_dataframe['Sequences'][ind]))
                
    def write_merged_filtered_fasta_file(self,WRITE,MERGED_FASTA):
        with open(MERGED_FASTA,"a+") as fh:
            for filename in WRITE:
                fr = open(filename,"r")
                fh.write(fr.read())
                fr.seek(0)
                fr.close()

if __name__=="__main__":
    
    print(__prog_name__ + ' v' + __version__ + ': ' + __prog_desc__)
    print('  by ' + __author__ + ' (' + __email__ + ')' + '\n')
    
    preprocess = PrepareForMsa()
    preprocess.main()
