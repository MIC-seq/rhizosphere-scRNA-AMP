import os
import phylophlan_checked as checked
engine=checked.engine
original=engine.gene_markers_identification
def mapped(configs,key,inputs,output_folder,database_name,database,nproc=1,verbose=False):
    return original(configs,key,inputs,output_folder,database_name,database,nproc=8,verbose=verbose)
engine.gene_markers_identification=mapped
if __name__=='__main__':engine.phylophlan_main()
