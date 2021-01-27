# Phylogenetic pipeline
Project for Comaparative Genomics classes

## Requirements
* MMseqs2
* ClustalW2 (should be placed in tools)
* Seqret
* PhyML
* Fasturec (should be placed in tools)

## Usage
To download and cluster protein sequences use argument cluster and provide path to text file containg informations about 
RefSeq assembly accession number and meaningful organism name (up to 9 chars) in one line.  
```python3 main.py -cluster accessions.txt```    

To calculate multialignment and ML trees, depending on approach use following arguments  
```python3 main.py -approach one_to_one```    
```python3 main.py -approach allow_paralogs```    
```python3 main.py -approach bootstrap```  
  
To calculate consensus tree provide path to dir containing PhyML ML trees and name for file    
```python3 main.py -consensus_tree ./data/clusters/cluster_one_to_one/ml one_to_one```  
  
To calcuate supertree, provide path to dir containing PhyML ML trees        
```python3 main.py -supertree ./data/clusters/cluster_one_to_one/ml```    
To calcuate supertree for bootstrap, provide path to dir containing PhyML ML trees and bootsrap cutoff (percent)      
```python3 main.py -supertree_bootstrap ./data/clusters/cluster_bootstrap/ml 80```  
      
To select best supertree, provide path to Fasturec output file,  name to the file and index of best tree    
```python3 -select_best_supertree 20210126.010008.fu.txt one_to_one 16```      

To visualize a tree in newick format provide path to a file    
```python3 main.py -visualize_newick best_supertree_one_to_one_score4455.txt```      


## Steps of the analysis
1. Download reference protein sequences from RefSeq NCBI  
2. MMseqs2 clustering  
3. Calculating multialignment
4. Calculating ML trees
5. Calculating consensus and supertrees