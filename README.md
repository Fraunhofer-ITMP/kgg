<p align="center">
  <img style="width: 200; height: 200px;" src="data/misc/KGG.png">
</p>

# Knowledge Graph Generator
**A fully automated workflow to create disease-specific Knowledge Graphs**

The Knowledge Graph Generator (KGG) workflow allows users to create KGs representing chemotype-phenotype of diseases of interest. The KGG is developed such that it is able to generate KGs with a minimum input (i.e., standard disease id) which users are prompted to identify at the beginning of the workflow. Additionally, the users can customize the size and content of KG with options to choose number of proteins and clinical trial phase of chemicals to be represented in the KG. The final KG is composed of disease-associated entities such as proteins, protein-related pathways, biological processes and functions, chemicals, mechanism of actions, assays and adverse effects. This is achieved by embedding underlying schema of curated databases (such as OpenTargets, Uniprot, ChEMBL and so on) which resemble a clockwork-esque mechanism ([Full paper](https://doi.org/10.1093/bioinformatics/btaf383)).  

# Workflow

The workflow is divided into 3 main phases as shown below:

![KGGworkflow](https://github.com/Fraunhofer-ITMP/kgg/blob/main/data/manuscript%20figures%20and%20files/Figure%201.png)

# KGG schema

The workflow can capture upto 10 types of entities and 11 types of relationships with entity specific annotations.

![KGGschema](https://github.com/Fraunhofer-ITMP/kgg/blob/main/data/misc/kggSchema.png)

# KGG demo

A demo of generating a KG is shown below.

![KGG_demo](https://github.com/Fraunhofer-ITMP/kgg/blob/main/data/misc/kgg_gif.gif)

# Usage guidelines

### Method 1: Dashboard version
#### *No pre-installations required, recommended for non-programmers* ####
This deployment (beta-version) of KGG is available at [SciLifeLab Serve](https://fraunhofer-itmp-ds-toolkit.serve.scilifelab.se/KGG) and ***does not require installation*** of python and relevant packages. Please select the "KG Generator" tab and follow the step-wise process to generate disease-specific KGs.

### Method 2: Local setup 
#### *Pre-installations required, intended for users with advanced programming skills* ####

#### Software requirements

    Operating system(s): Windows/Linux/Mac
    
    Programming language: Python 3.9.1 or higher
    
    Other requirements: Pre-installed Visual Studio Code (version 1.100.2, tested and stable) 
    
    License: MIT license

#### Cloning the repository and setting up environment

    git clone https://github.com/Fraunhofer-ITMP/kgg.git
    cd kgg
    conda create --name=kgg python=3.9
    conda activate kgg
    pip install -r requirements.txt

#### Quickstart
*Note: Please ensure that the **kgg** environment is activated.*

#### 1. Import required packages and dependencies. Do the following in VS Code: 

    from utils_v2 import *
    from kg_gen_5 import *
    
#### 2. Execute the `createKG` function which encapsulates multiple operations necessary for constructing a KG. It is a user-input driven multi-step workflow. Saving files and plots is possible at the end. 

    kg = createKG()

#### 3. Get summary of KG once files are saved.

    kg.summarize

#### 4. Visualize a sub-graph of random 250 edges. 

*Note: Please avoid visualizing entire KG in IPython Notebook. Only specific tools such as neo4j and cytoscape can handle large KGs.*

    
    to_jupyter(pybel.struct.mutation.induction.get_random_subgraph(kg))

#### 5. For analysis and evaluation of KGs, please refer to [pybel](https://pybel.readthedocs.io/en/latest/index.html) and [PyKEEN](https://github.com/pykeen/pykeen) documentations. 

#### Manuscript Results

The results included in the KGG manuscript are generated from the final KG files with `.pkl` format. Their usage in each of results are provided as indiviual IPython Notebook files in `src` folder.
1. [Comparison of AD-COVID19 KGs](https://github.com/Fraunhofer-ITMP/kgg/blob/main/src/Comparison%20of%20AD-COVID19%20KGs.ipynb)  
2. [Comparison of Depression KGs](https://github.com/Fraunhofer-ITMP/kgg/blob/main/src/Comparison%20of%20Depression%20KGs.ipynb)
3. [Comparison of Parkinson KGs](https://github.com/Fraunhofer-ITMP/kgg/blob/main/src/Comparison%20of%20Parkinson%20KGs.ipynb)

####  List of important and useful functions


1. Retrieve mechanism of action for drugs/chemicals

*Input: A list of ChEMBL identifiers :::: Output: A dictionary of mechanism of actions and target proteins*

    RetMech(chembl_ids)
    
2. Retrieve active assays (biological/functional, pChEMBL > 6) and target proteins for drugs/chemicals

*Input: A list of ChEMBL identifiers :::: Output: A dictionary*
    
    RetAct(chembl_ids) 

3. Map proteins represented as ChEMBL identifiers with UniProt identifiers and approved names

*Input: A list ChEMBL identifiers :::: Output: A dictionary of Uniprot ids and HGNC names*
    
    chembl2uniprot(chembl_ids)

4. Retrieve biological process, molecular functions and pathways for proteins

*Input: A list UniProt identifiers :::: Output: A dictionary*
    
    ExtractFromUniProt(uniprot_ids)
    
5. Get SMILES for drugs/chemicals

*Input: A list ChEMBL identifiers :::: Output: A dataframe of canonical SMILES*
    
    GetSmiles(chembl_ids)

6. Perform druglikeness assessment (Lipinski ro5, Ghose, Veber, REOS and QED properties) of drugs/chemicals
   
*Input: A dataframe from GetSmiles :::: Output: A dataframe with various physicochemical properties and flags for druglikeness*
    
    calculate_filters(dataframe,chembl_id_colname)

7. Convert CAS ids to CIDs (i.e. PubChem compound identifiers)

*Input: A list CAS ids :::: Output: A list of CIDs*
    
    cas2cid(cas_ids)

8. Convert CIDs to ChEMBL identifiers

*Input: A list CIDs  :::: Output: A list of ChEMBL ids*
    
    cid2chembl(cid_ids)

9. Create sub-graph

*Input: A list of desired entities i.e., protein, drug, etc.  :::: Output: A sub-graph with input entities and their 1st neighbors*
    
    filter_graph(mainGraph,list_of_entities)

10. Get drugs (FDA approved + clinical trials) and associated diseases for proteins

*Input: A list HGNC symbols  :::: Output: A dataframe of drugs diseases*

    getDrugsforProteins(protein_list)
