This part explains the quantification of the rawreads

1) First to MAP the raw reads, a data index had to be generated: 
	
In "index_build" you find the necessary steps to do that. I generated a conda environment named salmon-env. Because we use later for the quantification the tool called salmon, we also have to generate a specific index. 

The index have been saved at the shared drive found under "jlab/ValentinR/RNA_mapping/salmon_index_build/salmon_TAIR10_decoy_index_NEW2026".

2) The index can be used to now quanitify the rawreads with the "bulk_RNA_salmon_quantification_workflow_NEW_2026" script

The rawread data is stored on the big Computer in the informatic room and the code has been running there, because it also needs more computational power to generate it for all 100 samples. 


