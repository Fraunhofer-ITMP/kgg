# -*- coding: utf-8 -*-

"""Utils files with all functions relevant to generation of KG."""
import os
from IPython.display import Image
import logging
import pickle
from collections import defaultdict
import networkx as nx
import pandas as pd
import pubchempy as pcp
import requests
from chembl_webresource_client.http_errors import HttpBadRequest, HttpApplicationError
from chembl_webresource_client.new_client import new_client
from pybel import BELGraph
from pybel.dsl import Protein, Abundance, Pathology, BiologicalProcess, Gene
from tqdm.auto import tqdm
import time
import seaborn as sns
import pybel
import json
from utils import *

import matplotlib.pyplot as plt
import matplotlib.image as mpimg

from IPython.display import display, HTML

from IPython.display import Markdown, display

logger = logging.getLogger("__name__")

def searchDisease(name): 
    
    disease_name = str(name)

    query_string = """
        query searchAnything ($disname:String!){
          search(queryString:$disname,entityNames:"disease",page:{size:20,index:0}){
            total
            hits {
              id
              entity
              name
              description
            }
          }
        }
    """

    variables = {"disname": disease_name}

    # Set base URL of GraphQL API endpoint
    base_url = "https://api.platform.opentargets.org/api/v4/graphql"

    # Perform POST request and check status code of response
    r = requests.post(base_url, json={"query": query_string, "variables": variables})
    #print(r.status_code)

    # Transform API response from JSON into Python dictionary and print in console
    api_response = json.loads(r.text)
    #print(api_response)

    df = pd.DataFrame(api_response['data']['search']['hits'])
    
    if not df.empty:
        df = df.drop(columns=['entity', 'description'])
        df['index'] = df.index
        df.style.hide(axis='index')
        df = df[['index','id','name']]
        df.style.hide(axis='index')  
        #print('\n')
        #print(df.head(20))
        display(HTML(df.to_html(index=False)))
        
    # else:
        # print('\n')
        # print('Ooops!! Did you have a typo in the name. Please try again!')
        # Generate_KG()
        
    return(df)
     
def printmd(string, color=None):
    colorstr = "<span style='color:{}'>{}</span>".format(color, string)
    display(Markdown(colorstr))

def GetImage():
    
    from PIL import Image
    import requests
    import matplotlib.pyplot as plt

    img = "https://github.com/Fraunhofer-ITMP/kgg/blob/main/data/KGG.png?raw=true"
    response = requests.get(img, stream=True)
    img = Image.open(response.raw)
    
    fig = plt.figure(figsize=(3,3))
    plt.imshow(img)
    plt.axis('off')
    plt.show()

def GetQuery():
    
    GetImage()

    printmd("**Welcome to the KG Generator tool. In the following steps, we will need some inputs from your side.**",color = "blue")

    query = input('Please enter the disease you are interested in and we will try to find the best matches for you.' +'\n' + '\n' + 'Input: ')
    return(query)

def Generate_KG():
 
    #this is required to get above printed before the input is asked
    time.sleep(0.05)
    
    #dis = input('Please enter the disease you are interested in and we will try to find the best matches for you.' +'\n' + '\n' + 'Input: ')
    
    dis = GetQuery()
    
    #temp = searchDisease(query_disease)
    temp = searchDisease(dis)
    #print(temp)
    #print('\n')
    if not temp.empty:
        printmd('**Here you go! Hopefully your disease of interest is in the list. If so, let\'s get started.**')
        #print('\n')
        #print(temp)
        return(dis,temp)
    else: 
        print('Ooops!! Did you have a typo in the name. Please try again!')
        return(Generate_KG())
        
    #return(temp)
        
def GetDiseaseAssociatedProteins(disease_id):   
    
    #efo_id = input('Please enter the id of your disease of interest. Input: ')
    #print('\n')
    #print('Just a second please!')
    #print('\n')
    
    #import numpy as np
    import matplotlib.pyplot as plt
    efo_id = str(disease_id)

    query_string = """
        query associatedTargets{
          disease(efoId: $efo_id){
            id
            name
            associatedTargets(page:{size:15000,index:0}){
              count
              rows {
                target {
                  id
                  approvedSymbol
                  proteinIds {
                    id
                    source
                  }
                }
                score
              }
            }
          }
        }

    """
    
    #replace $efo_id with value from efo_id
    query_string = query_string.replace("$efo_id",f'"{efo_id}"')
    
    #variables = {"$efo_id":efo_id}

    # Set base URL of GraphQL API endpoint
    base_url = "https://api.platform.opentargets.org/api/v4/graphql"

    # Perform POST request and check status code of response
    #r = requests.post(base_url, json={"query": query_string, "variables": variables})
    r = requests.post(base_url, json={"query": query_string})
    #print(r.status_code)

    # Transform API response from JSON into Python dictionary and print in console
    api_response = json.loads(r.text)
    
    temp_list = []
    for item in api_response['data']['disease']['associatedTargets']['rows']:
        #print(item['target'])
        #break
        for obj in item['target']['proteinIds']:
            if obj['source'] == 'uniprot_swissprot':
                #print(obj)
                uprot = obj['id']
                source = obj['source']
                score = item['score']
                ensg = item['target']['id']
                name = item['target']['approvedSymbol']
                temp = {'Protein':name,'ENSG':ensg,'UniProt':uprot,'Source':source,'Score':score}
                temp_list.append(temp)
    
    df = pd.DataFrame(temp_list)
    
    
    print('We have identified ' + str(len(df)) + ' proteins (Swiss-Prot) associated with the disease. Following is a histogram that shows '
         + 'distribution of proteins based on scores provided by OpenTargets. The scores are influenced by various factors '
         + 'such as genetic associations, expression, mutations, known pathways, targeting drugs and so on.'+'\n')
    
    print('Displaying top 20 genes')
    df_display = df.head(20)
    display(HTML(df_display.to_html(index=False)))
    #print(df,'\n')
    #print('\n')
    
    fig, ax = plt.subplots()
    ax.hist(df['Score'])
    ax.set_title('Distribution of proteins based on OpenTargets score')
    ax.set_xlabel('Score')
    ax.set_ylabel('No. of proteins')

    fig.tight_layout()
    plt.show()
    
    print('\n')
    
    time.sleep(0.05)
    
    score = input('We recommend taking a threshold above 0.3 to exclude loosely associated proteins. ' + '\n' +'Please enter your desired threshold: ')
    
    df = df.loc[df['Score'] >= float(score),:]
    
    print('\n')
    print('Alright, we are good to go now. Your KG is now being generated! Sit back and relax!!')
    
    
    print('\n','Total no. of proteins: ',len(df))
   
    #display(HTML(df.to_html(index=False)))
    print('\n',df)
    #print('\n')
    return(df)
    
def getDrugCount(disease_id):
    
    efo_id = disease_id

    query_string = """
        query associatedTargets($my_efo_id: String!){
          disease(efoId: $my_efo_id){
            id
            name
            knownDrugs{
                uniqueTargets
                uniqueDrugs
                count
            }
          }
        }

    """

    # Set variables object of arguments to be passed to endpoint
    variables = {"my_efo_id": efo_id}

    # Set base URL of GraphQL API endpoint
    base_url = "https://api.platform.opentargets.org/api/v4/graphql"

    # Perform POST request and check status code of response
    r = requests.post(base_url, json={"query": query_string, "variables": variables})

    # Transform API response from JSON into Python dictionary and print in console
    api_response = json.loads(r.text)
    
    #get the count value from api_repsonse dict 
    api_response = api_response['data']['disease']['knownDrugs']['count']
    return(api_response)
    
def GetDiseaseAssociatedDrugs(disease_id,CT_phase):

    efo_id = disease_id
    size = getDrugCount(efo_id)

    query_string = """
        query associatedTargets($my_efo_id: String!, $my_size: Int){
          disease(efoId: $my_efo_id){
            id
            name
            knownDrugs(size:$my_size){
                uniqueTargets
                uniqueDrugs
                count
                rows{
                    approvedSymbol
                    approvedName
                    prefName
                    drugType
                    drugId
                    phase
                    ctIds
                }

            }
          }
        }

    """

    #replace $efo_id with value from efo_id
    #query_string = query_string.replace("$efo_id",f'"{efo_id}"')
    #query_string = query_string.replace("$efo_id",f'"{efo_id}"')

    # Set variables object of arguments to be passed to endpoint
    variables = {"my_efo_id": efo_id, "my_size": size}

    # Set base URL of GraphQL API endpoint
    base_url = "https://api.platform.opentargets.org/api/v4/graphql"

    # Perform POST request and check status code of response
    r = requests.post(base_url, json={"query": query_string, "variables": variables})
    #r = requests.post(base_url, json={"query": query_string})
    #print(r.status_code)

    # Transform API response from JSON into Python dictionary and print in console
    api_response = json.loads(r.text)

    df = pd.DataFrame(api_response['data']['disease']['knownDrugs']['rows'])
    df = df.loc[df['phase'] >= int(CT_phase),:]
    
    if not df.empty:
        df['id'] = efo_id
        df['disease'] = api_response['data']['disease']['name']
        #print('Your dataframe is ready')
        return(df)
    
    else:
        print('No drugs found in clinical trials')
    
    #chembl_list = list(set(df['drugId']))
    #return(chembl_list)
    
def KG_namespace_plot(final_kg,kg_name):
    
    import matplotlib.pyplot as plt
    nspace_count = pybel.struct.summary.count_namespaces(final_kg)
    nspace_count = dict(nspace_count)

    nspace_data = {'Namespace':list(nspace_count.keys()),'Number':list(nspace_count.values())}
    nspace = pd.DataFrame(nspace_data)
    plt.figure()

    a = sns.barplot(x="Number", y="Namespace", data=nspace_data)
    a.set(xlabel='Number',ylabel='Namespace',title= 'KG Namespace in numbers')

    plt.tight_layout()
    plt.savefig(kg_name+'_namespace.png',dpi=600)
    #plt.show()    

def getAdverseEffectCount(chembl_id):
    
    get_id = chembl_id

    query_string = """
        query AdverseEventsQuery(
          $chemblId: String!
          $index: Int = 0
          $size: Int = 10
        ) {
          drug(chemblId: $chemblId) {
            id
            maxLlr: adverseEvents(page: { index: 0, size: 1 }) {
              rows {
                logLR
              }
            }
            adverseEvents(page: { index: $index, size: $size }) {
              criticalValue
              count
              rows {
                name
                count
                logLR
                meddraCode
              }
            }
          }
        }

    """

    # Set variables object of arguments to be passed to endpoint
    variables = {"chemblId": get_id}

    # Set base URL of GraphQL API endpoint
    base_url = "https://api.platform.opentargets.org/api/v4/graphql"

    # Perform POST request and check status code of response
    r = requests.post(base_url, json={"query": query_string, "variables": variables})

    # Transform API response from JSON into Python dictionary and print in console
    api_response = json.loads(r.text)
    
    #get the count value from api_repsonse dict 
    api_response = api_response['data']['drug']['adverseEvents']['count']
    return(api_response)
    
def GetAdverseEvents(chem_list):
    
    api_response = pd.DataFrame()
    
    for chem in tqdm(chem_list, desc = 'Retrieving Adverse Effects for each drug'):
        
        chembl_id = chem
        
        try:
        
            #get total no. of adverse effects for a given drug
            count = getAdverseEffectCount(chembl_id)

            query_string = """
                query AdverseEventsQuery(
                  $chemblId: String!
                  $index: Int = 0
                  $size: Int!
                ) {
                  drug(chemblId: $chemblId) {
                    id
                    maxLlr: adverseEvents(page: { index: 0, size: 1 }) {
                      rows {
                        logLR
                      }
                    }
                    adverseEvents(page: { index: $index, size: $size }) {
                      criticalValue
                      count
                      rows {
                        name
                        count
                        logLR
                        meddraCode
                      }
                    }
                  }
                }

        """
            # Set variables object of arguments to be passed to endpoint
            variables = {"chemblId": chembl_id, "size": count}

            # Set base URL of GraphQL API endpoint
            base_url = "https://api.platform.opentargets.org/api/v4/graphql"

            # Perform POST request and check status code of response
            r = requests.post(base_url, json={"query": query_string, "variables": variables})
            #r = requests.post(base_url, json={"query": query_string})
            #print(r.status_code)

            # Transform API response from JSON into Python dictionary and print in console
            api_response_temp = json.loads(r.text)

            api_response_temp = api_response_temp['data']['drug']['adverseEvents']['rows']
            api_response_temp = pd.DataFrame(api_response_temp)
            api_response_temp ['chembl_id'] = chembl_id

            api_response = pd.concat([api_response,api_response_temp])
            
        except:
            continue
            
    api_response.reset_index(drop=True, inplace=True)
    return(api_response)
        
def chembl2adverseEffect_rel(
    chembl_adveff_df,
    graph: BELGraph
) -> BELGraph:
    """

    :param chembl_adveff_df:
    :param graph:
    :return:
    """

    for i in range(len(chembl_adveff_df)):

        graph.add_association(
            Abundance(namespace='ChEMBL', name= str(chembl_adveff_df['chembl_id'][i])),
            Pathology(namespace='SideEffect', name= str(chembl_adveff_df['name'][i])),  # TODO: Fix namespace
            citation="OpenTargets Platform",
            evidence='DrugReactions'
        )

    return graph

def GetViralProteins(query_disease):
    
    # file downloaded from https://www.genome.jp/ftp/db/virushostdb Dated: 12/09/2023
    virus = pd.read_csv('https://raw.githubusercontent.com/Fraunhofer-ITMP/kgg/main/data/virushostdb.csv')
    
    cols = ['virus tax id','virus name','DISEASE','host tax id']
    virus = virus[cols]
    
    #print(virus)

    #filter virus with host humans 
    virus = virus.loc[virus['host tax id'] == 9606.0,:]
    virus = virus.reset_index(drop=True)
    
    #replace 9606.0 to 9606
    virus["host tax id"] = pd.to_numeric(virus["host tax id"], downcast='integer')
    
    #get the initial keyword for disease search
    #disease = GetQuery()
    
    #subset df with disease keyword
    virus_subset_1 = virus[virus['DISEASE'].str.contains(query_disease,na=False,case=False)]
    
    if not virus_subset_1.empty:

        #print(disease)
        print('\n')
        print('The workflow has identified your query as a viral disease. Its proteins (SWISS-Prot) will be now represented in the KG.','\n')

        time.sleep(0.1)

        virus_name = input('Do you want to look further for a specific virus? Please type its name or skip it by typing \'no\': ')

        #subset df with virus name
        if virus_name.lower() != 'no':

            virus_subset_2 = virus[virus['virus name'].str.contains(virus_name,na=False,case=False)]
            #print(virus_subset_2)
            #break
        else:
            virus_subset_2 = pd.DataFrame()

        #merge subsets of df_1 and df_2
        virus_subset_merge = pd.concat([virus_subset_1,virus_subset_2])

        virus_subset_merge = virus_subset_merge.drop_duplicates(keep='first')

        virus_subset_merge = virus_subset_merge.reset_index(drop=True)

        virus_subset_merge['index'] = virus_subset_merge.index

        virus_subset_merge.style.hide(axis='index')  

        virus_subset_merge = virus_subset_merge[['index','virus tax id','virus name','DISEASE','host tax id']]
        display(HTML(virus_subset_merge.to_html(index=False)))
        

        time.sleep(0.1)
        temp_id = input('Enter the index value(s). If multiple, use space, for example -> 0 1 3: ')

        print('\n')

        temp_id = temp_id.split(' ')
        temp_id = [int(x) for x in temp_id]
        #print(virus_subset_merge.loc[0]['virus tax id'])

        uprot_list = []

        for item in temp_id:

            tax_id = virus_subset_merge.loc[item]['virus tax id']
            #print(tax_id)

            #fetch tax id related proteins from Uniprot
            #the link can be created from downloads option in uniprot
            query_string = 'https://rest.uniprot.org/uniprotkb/stream?fields=accession%2Creviewed%2Cid%2Cgene_names%2Corganism_name%2Clength%2Cgene_primary%2Cprotein_name&format=tsv&query=%28%28taxonomy_id%3A'+str(tax_id)+'%29+AND+%28reviewed%3Atrue%29%29'

            #query_string = 'https://rest.uniprot.org/uniprotkb/stream?fields=accession%2Creviewed%2Cid%2Cgene_names%2Corganism_name%2Clength%2Cgene_primary%2Cprotein_name&format=tsv&query=%28%28taxonomy_id%3A11676%29%29+AND+%28reviewed%3Atrue%29'

            query_uniprot = requests.get(query_string)
            query_uniprot = query_uniprot.text.split('\n')

            query_uniprot_df = pd.DataFrame([x.strip().split('\t') for x in query_uniprot])
            cols = query_uniprot_df.iloc[0]
            #print(cols)
            query_uniprot_df = query_uniprot_df[1:len(query_uniprot_df)-1]
            query_uniprot_df.columns = cols
            temp = list(query_uniprot_df['Entry'])
            #print(len(temp))
            uprot_list.append(temp)

        uprot_list = [item for sublist in uprot_list for item in sublist]
        
        print('A total of',str(len(uprot_list)), 'viral proteins have been identified.','\n')
        
        return(uprot_list)

def saveFiles(kgName, disease2protein, drugAdvEffect, final_kg, drug_df):

    print('Now let\'s save all the files that were created in the process.','\n')
    
    print('Please enter the location (e.g. \'C:\\Users\\rkarki\\Documents\\kg\\\' ) where KG files should be stored. A folder will be created automatically.','\n')
    
    time.sleep(0.2)
    
    path = input('Input: ')
    
    path = path+kgName
    
    os.makedirs(path,exist_ok=True)
    
    os.chdir(path)
    
    disease2protein.to_csv('diseaseAssociatedProteins.csv',sep=',')
    
    drugAdvEffect.to_csv('adverseEffects.csv',sep=',')
    
    drug_df.to_csv('diseaseAssociatedDrugs.csv',sep=',')
    
    #to cytoscape compatible graphml 
    pybel.to_graphml(final_kg,kgName+'.graphml')

    #to regular BEL format
    pybel.dump(final_kg,kgName+'.bel')

    #to neo4j
    pybel.to_csv(final_kg,kgName+'.csv')
    
    filename = kgName + '.pkl'
    outfile = open(filename,'wb')
    pickle.dump(final_kg,outfile)
    outfile.close()
    
    #plot for namespace dist
    KG_namespace_plot(final_kg,kgName)

def GetDiseaseSNPs(disease_id): 
    
    import requests
    
    doid = disease_id.split("_")
    #print(doid)
    
    #For this example we are going to use the python default http library
    

    #Build a dict with the following format, change the value of the two keys your DisGeNET account credentials, if you don't have an account you can create one here https://www.disgenet.org/signup/ 
    auth_params = {"email":"reagonkarki@gmail.com","password":"Bhunti.87"}

    api_host = "https://www.disgenet.org/api"

    api_key = 'e25cb13382cb9b016247822c49f325f75991e607'
    s = requests.Session()

    if api_key:
        #Add the api key to the requests headers of the requests Session object in order to use the restricted endpoints.
        s.headers.update({"Authorization": "Bearer %s" % api_key}) 
        #Lets get all the diseases associated to a gene eg. APP (EntrezID 351) and restricted by a source.

        #https://www.disgenet.org/api/vda/disease/D000544
    
        gda_response = s.get(api_host+'/vda/disease/'+str(doid[0]).lower()+ "/" +str(doid[1]) +'?format=json')
        #gda_response = s.get(api_host+'/gda/gene/351', params={'source':'UNIPROT'})

        if gda_response:
        
            gda_response = gda_response.json()
            gda_response = pd.DataFrame(gda_response)
            return(gda_response)

    if s:
        s.close()

def snp2gene_rel(snp_df,graph): 

    
    if 'errors' in snp_df.columns:
        print('No SNPs found')
    
    else:
        
        print('A total of ' + str(len(snp_df)) + ' SNPs have been identified from DisGeNET. Now adding relevant data')
        print('\n')
        #convert col to datatype == str to remove rows that have 'None' in gene_symbol
        snp_df[['gene_symbol']] =  snp_df[['gene_symbol']].astype(str)
        snp_df = snp_df.loc[snp_df['gene_symbol'] != "None"]
        snp_df = snp_df.reset_index(drop=True)
        
        #print(len(snp_df))
        
        for i in tqdm(range(len(snp_df)),desc='Adding disease associated SNPs'):
            genes = snp_df['gene_symbol'][i].split(';')
            
            for j in range(len(genes)):
            
                graph.add_association(
                    Gene(namespace="dbSNP",name=snp_df['variantid'][i]),
                    Protein(namespace = "HGNC", name = genes[j]),
                    citation = "DisGeNet",
                    evidence = "SNPs for queried disease"
                )
        
        return(graph)
    
def createKG():

    #import matplotlib.pyplot as plt
    #import matplotlib.image as mpimg
    # image = mpimg.imread("data/KGG.png")
    # plt.imshow(image)
    # plt.axis('off')
    # plt.show()
    
    #query = GetQuery()
    
    #doid = Generate_KG(query)
    query,doid = Generate_KG()
    #break
    
    #print(doid)
    
    time.sleep(0.1)
    
    efo_id = int(input('Please enter the index value of your disease of interest. Input: '))
    #print(efo_id)
    
    vir_prot = GetViralProteins(query)
    
    print('\n')

    print('Please enter the clinical trial phase of chemicals which should be identified by the workflow. Use a number between 1 (early phase) and 4 (FDA approved). For example, if you use 3, the KG will fetch chemicals that are in phase 3. Also, remember that lower the input value, higher will be the number of identified chemicals and therefore the running time of workflow also increases.')
    
    time.sleep(0.05)
    
    print('\n')
    ct_phase = input('Your desired clinical trial phase: ')
    print('\n')
   
    kg_name = input('Please provide a name for you KG. Input: ')
    
    print('\n')
    
    #print(doid['id'][efo_id])
    
    #df_dis2prot = GetDiseaseAssociatedProteins(efo_id)
    
    df_dis2prot = GetDiseaseAssociatedProteins(doid['id'][efo_id])
    
    #chembl_list = GetDiseaseAssociatedDrugs(efo_id,ct_phase)
    
    chembl_df = GetDiseaseAssociatedDrugs(doid['id'][efo_id],ct_phase)
       
    #create empty KG
    kg = pybel.BELGraph(name=kg_name, version="0.0.1")
    
    uprot_ext = ExtractFromUniProt(df_dis2prot['UniProt'])
    
    if vir_prot:
        vir_uprot_ext = ExtractFromUniProt(vir_prot)
    
    print('A total of ' + str(len(list(set(chembl_df['drugId'])))) + ' drugs have been identified. Now fetching relevant data')
    
    chembl2mech = RetMech(list(set(chembl_df['drugId'])))
    chembl2act = RetAct(list(set(chembl_df['drugId'])))
    
    prtn_as_chembl = Ret_chembl_protein(chembl2act) + Ret_chembl_protein(chembl2mech)
    prtn_as_chembl = set(prtn_as_chembl)
    prtn_as_chembl = list(prtn_as_chembl)
    chembl2uprot = chembl2uniprot(prtn_as_chembl)
    
    chembl2act = chembl2gene2path(chembl2uprot, chembl2act)
    chembl2mech = chembl2gene2path(chembl2uprot, chembl2mech)
    
    
    kg = uniprot_rel(uprot_ext, 'HGNC', kg)
    
    if vir_prot:
        kg = uniprot_rel(vir_uprot_ext,'VP',kg)
    
    kg = chem2moa_rel(chembl2mech, 'HGNC', kg)
    kg = chem2act_rel(chembl2act, 'HGNC', kg)
    kg = gene2path_rel(chembl2uprot, 'HGNC', kg)
    
    adv_effect = GetAdverseEvents(list(set(chembl_df['drugId'])))
    kg = chembl2adverseEffect_rel(adv_effect,kg)
    
    snp_dgnet = GetDiseaseSNPs(doid['id'][efo_id])  
    #print('A total of ' + str(len(snp_dgnet)) + ' SNPs have been identified from DisGeNET. Now adding relevant data')
    
    if snp_dgnet:
        kg = snp2gene_rel(snp_dgnet,kg)
    
    print('Your KG is now generated!','\n')
    
    saveFiles(kg_name, df_dis2prot, adv_effect, kg, chembl_df)
    
    return(kg)
    

# def getNodeList(nodeName,graph):
    # # import pybel
    # node_list = []
    # for node in graph.nodes():
        # if isinstance(node,pybel.dsl.Gene):
            # if node.namespace == nodeName:
                # node_list.append(node.name)
    # return(node_list)

# def dbSNP_annotation(graph):
    # snps = getNodeList('dbSNP',graph)
    # print(snps)
    # for item in snps:
        # nx.set_node_attributes(graph,{Gene(namespace='dbSNP',
        # name=item):'https://www.ncbi.nlm.nih.gov/snp/'+item},'source')
        
    # return(graph)

    
           
