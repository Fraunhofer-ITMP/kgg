# -*- coding: utf-8 -*-

"""Utils files with all functions relevant to generation of KG."""
import io
import os
import sys
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
from utils_v2 import *
from rdkit.Chem.SaltRemover import SaltRemover

import matplotlib.pyplot as plt
import matplotlib.image as mpimg

from IPython.display import display, HTML

from IPython.display import Markdown, display
from PIL import Image
import requests
import matplotlib.pyplot as plt
import requests
import json
import pandas as pd
import matplotlib.pyplot as plt
import os
from tqdm.auto import tqdm

import pybel
from pybel.io.jupyter import to_jupyter
from IPython.display import Markdown, display

from rdkit import Chem
from rdkit.Chem import Descriptors

from pandasgwas import get_variants
from pandasgwas.get_variants import get_variants_by_efo_id

import animation
from PIL import Image
import plotly
import plotly.graph_objects as go
import plotly.express as px

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

    img = "https://github.com/Fraunhofer-ITMP/kgg/blob/main/data/misc/KGG.png?raw=true"
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
        
        
wheel = ('-', '/', '|', '\\')
@animation.wait(wheel,text = '    Fetching Disease-specific Proteins ')        
# def GetDiseaseAssociatedProteins(disease_id):   
    
    # #efo_id = input('Please enter the id of your disease of interest. Input: ')
    # #print('\n')
    # #print('Just a second please!')
    # #print('\n')
    
    # #import numpy as np
    # import matplotlib.pyplot as plt
    # efo_id = str(disease_id)

    # query_string = """
        # query associatedTargets{
          # disease(efoId: $efo_id){
            # id
            # name
            # associatedTargets(page:{size:15000,index:0}){
              # count
              # rows {
                # target {
                  # id
                  # approvedSymbol
                  # proteinIds {
                    # id
                    # source
                  # }
                # }
                # score
              # }
            # }
          # }
        # }

    # """
    
    # #replace $efo_id with value from efo_id
    # query_string = query_string.replace("$efo_id",f'"{efo_id}"')
    
    # #variables = {"$efo_id":efo_id}

    # # Set base URL of GraphQL API endpoint
    # base_url = "https://api.platform.opentargets.org/api/v4/graphql"

    # # Perform POST request and check status code of response
    # #r = requests.post(base_url, json={"query": query_string, "variables": variables})
    # r = requests.post(base_url, json={"query": query_string})
    # #print(r.status_code)

    # # Transform API response from JSON into Python dictionary and print in console
    # api_response = json.loads(r.text)
    
    # temp_list = []
    # for item in api_response['data']['disease']['associatedTargets']['rows']:
        # #print(item['target'])
        # #break
        # for obj in item['target']['proteinIds']:
            # if obj['source'] == 'uniprot_swissprot':
                # #print(obj)
                # uprot = obj['id']
                # source = obj['source']
                # score = item['score']
                # ensg = item['target']['id']
                # name = item['target']['approvedSymbol']
                # temp = {'Protein':name,'ENSG':ensg,'UniProt':uprot,'Source':source,'Score':score}
                # temp_list.append(temp)
    
    # df = pd.DataFrame(temp_list)
    # df['disease_id'] = efo_id
    
    # return(df)
    
def GetDiseaseAssociatedProteins(disease_id,index_counter=0,merged_list= []):   
    
    efo_id = str(disease_id)

    query_string = """
        query associatedTargets($efoId: String!,$index:Int!){
          disease(efoId: $efoId){
            id
            name
            associatedTargets(page:{size:3000,index:$index}){
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
    #query_string = query_string.replace("$efoId",f'"{efo_id}"')
    
    variables = {"efoId":efo_id,"index":index_counter}

    # Set base URL of GraphQL API endpoint
    base_url = "https://api.platform.opentargets.org/api/v4/graphql"

    # Perform POST request and check status code of response
    r = requests.post(base_url, json={"query": query_string, "variables": variables})
    #r = requests.post(base_url, json={"query": query_string})
    #print(r.status_code)

    # Transform API response from JSON into Python dictionary and print in console
    api_response = json.loads(r.text)
    
    result = api_response['data']['disease']['associatedTargets']['rows']
    #print(len(result))
    
    merged_list.extend(result)
    
    if result:
        counter = index_counter+1
        #print('Counter',counter)
        GetDiseaseAssociatedProteins(disease_id,counter,merged_list)
    
    #return(merged_list)
    
    temp_list = []
    for item in merged_list:
        #api_response['data']['disease']['associatedTargets']['rows']
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
    df['disease_id'] = efo_id
    
    return(df)


def GetDiseaseAssociatedProteinsPlot(df):

    print('We have identified ' + str(len(df)) + ' proteins (Swiss-Prot) associated with the disease. Please note that the proteins identified may not be unique if you combined two or more diseases. Following is a histogram that shows ' + 'distribution of proteins based on scores provided by OpenTargets. The scores are influenced by various factors ' + 'such as genetic associations, expression, mutations, known pathways, targeting drugs and so on.'+'\n')
    
    print('A total of ' +str(len(list(set(df['UniProt'])))) + ' unique proteins have been identified.')
    
    print('\n')
    
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

#key_gen_3 version    
# wheel = ('-', '/', '|', '\\')
# @animation.wait(wheel,text = '    Fetching Disease-specific Drugs ')
# def GetDiseaseAssociatedDrugs(disease_id,CT_phase):

    # efo_id = disease_id
    # size = getDrugCount(efo_id)

    # query_string = """
        # query associatedTargets($my_efo_id: String!, $my_size: Int){
          # disease(efoId: $my_efo_id){
            # id
            # name
            # knownDrugs(size:$my_size){
                # uniqueTargets
                # uniqueDrugs
                # count
                # rows{
                    # approvedSymbol
                    # approvedName
                    # prefName
                    # drugType
                    # drugId
                    # phase
                    # ctIds
                # }

            # }
          # }
        # }

    # """

    # #replace $efo_id with value from efo_id
    # #query_string = query_string.replace("$efo_id",f'"{efo_id}"')
    # #query_string = query_string.replace("$efo_id",f'"{efo_id}"')

    # # Set variables object of arguments to be passed to endpoint
    # variables = {"my_efo_id": efo_id, "my_size": size}

    # # Set base URL of GraphQL API endpoint
    # base_url = "https://api.platform.opentargets.org/api/v4/graphql"

    # # Perform POST request and check status code of response
    # r = requests.post(base_url, json={"query": query_string, "variables": variables})
    # #r = requests.post(base_url, json={"query": query_string})
    # #print(r.status_code)

    # # Transform API response from JSON into Python dictionary and print in console
    # api_response = json.loads(r.text)

    # df = pd.DataFrame(api_response['data']['disease']['knownDrugs']['rows'])
    # df = df.loc[df['phase'] >= int(CT_phase),:]
    
    # if not df.empty:
        # df['id'] = efo_id
        # df['disease'] = api_response['data']['disease']['name']
        # #print('Your dataframe is ready')
        # return(df)
    
    # else:
        # print('No drugs found in clinical trials')
        
   
wheel = ('-', '/', '|', '\\')
@animation.wait(wheel,text = '    Fetching Disease-specific Drugs ')
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
    
    #df = df.loc[df['phase'] >= int(CT_phase),:]
    
    if not df.empty:
        df = df.loc[df['phase'] >= int(CT_phase),:]
        df['id'] = efo_id
        df['disease'] = api_response['data']['disease']['name']
        #print('Your dataframe is ready')
        return(df)
    
    else:
        print('No drugs found in clinical trials')
        
        
# def KG_namespace_plot(final_kg,kg_name):
#
#     import matplotlib.pyplot as plt
#     nspace_count = pybel.struct.summary.count_namespaces(final_kg)
#     nspace_count = dict(nspace_count)
#
#     nspace_data = {'Namespace':list(nspace_count.keys()),'Number':list(nspace_count.values())}
#     nspace = pd.DataFrame(nspace_data)
#     plt.figure()
#
#     a = sns.barplot(x="Number", y="Namespace", data=nspace_data)
#     a.set(xlabel='Number',ylabel='Namespace',title= 'KG Namespace in numbers')
#
#     plt.tight_layout()
#     plt.savefig(kg_name+'_namespace.png',dpi=600)

def KG_namespace_plot(final_kg, kg_name):
    df = dict(final_kg.count.namespaces())
    df = pd.DataFrame([df]).T
    df = df.reset_index()
    df.columns = ['Namespace', 'Numbers']
    fig = px.bar(df, x='Namespace', y='Numbers', color='Namespace')
    fig.update_layout(title_text='KG namespace in numbers', title_x=0.5)

    path_img = kg_name + '_namespace.png'
    fig.write_image(path_img, scale=5)

    path_html = kg_name + '_namespace.html'
    fig.write_html(path_html)

    fig.show()

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

        graph.add_qualified_edge(
            Abundance(namespace='ChEMBL', name= str(chembl_adveff_df['chembl_id'][i])),
            Pathology(namespace='SideEffect', name= str(chembl_adveff_df['name'][i])),  
            relation='hasSideEffect',
            citation="OpenTargets Platform",
            evidence='DrugReactions'
        )

    return graph

def GetViralProteins(query_disease):
    
    # file downloaded from https://www.genome.jp/ftp/db/virushostdb Dated: 12/09/2023
    #virus = pd.read_csv('https://raw.githubusercontent.com/Fraunhofer-ITMP/kgg/main/data/misc/virushostdb.csv')
    virus = pd.read_csv('../data/misc/virushostdb.csv')
    
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

def saveFiles(kgName, disease2protein, drugAdvEffect, final_kg, drug_df,snp_df,uprot_dict):

    print('Now let\'s save all the files that were created in the process.','\n')

    print('Please enter the location where KG files should be stored. A folder will be created automatically.','\n')
    print('For local folders, use this example: C:\\Users\\rkarki\\Documents\\kg\\','\n')
    print('For servers like mybinder or EGI notebook, use this: /home/jovyan/data/kgs/test/','\n')
    
    time.sleep(0.2)
    
    path = input('Input: ')
    
    path = path+kgName
    
    #original_path = sys.path[0]
    original_path = os.getcwd()
    
    os.makedirs(path,exist_ok=True)
    #print(original_path)
    
    os.chdir(path)
    
    disease2protein.to_csv('diseaseAssociatedProteins.csv',sep=',')
    
    drugAdvEffect.to_csv('adverseEffects.csv',sep=',')
    
    drug_df.to_csv('diseaseAssociatedDrugs.csv',sep=',')
    
    snp_df.to_csv('diseaseAssociatedSNPs.csv',sep=',')
    
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
    
    filename = kgName+ '_prot_dict' + '.pkl'
    outfile = open(filename,'wb')
    pickle.dump(uprot_dict,outfile)
    outfile.close()
    
    #plot for namespace dist
    KG_namespace_plot(final_kg,kgName)
    print('Process Completed. Please check your folder.')

    os.chdir(original_path)

wheel = ('-', '/', '|', '\\')
@animation.wait(wheel,text = '    Fetching SNPs. Can take a while!')
def GetDiseaseSNPs(disease_id): 
    
    try:
        snps = get_variants_by_efo_id(disease_id)
        snps_df = snps.genomic_contexts
        snps_df['disease_id'] = disease_id
        #snps_functional_class_df = snps.variants
        #snps_functional_class_df['disease_id'] = disease_id

        return(snps_df)
    
    except:
        print('No SNPs found. This could be either because 1. The identifier of your disease of interest is not compatible with Experimental Factor Ontology (EFO) or 2. No SNPs have been reported for the disease')
    

def snp2gene_rel(snp_df,graph): 
    
    #print(snp_df.head(2))
    
    kg_prots = getProtfromKG(graph)

    unique_prots_df = pd.DataFrame(kg_prots,columns=['Proteins'])
    
    snp_df = unique_prots_df.merge(snp_df, how='inner', left_on='Proteins', right_on='gene.geneName')
    
    
    #take SNPs which are within the gene sequence
    snp_df = snp_df.loc[snp_df['distance'] == 0]
    snp_df = snp_df.reset_index(drop=True)

        
    print('A total of ' + str(len(snp_df)) + ' SNPs have been identified from GWAS Central. Now adding relevant data')
    print('\n')

    
    #graph = pybel.BELGraph(name='test', version="0.0.1")
    
    
    for i in tqdm(range(len(snp_df)),desc='Adding disease associated SNPs'):
    
        graph.add_qualified_edge(
        Gene(namespace="dbSNP",name=snp_df['rsId'][i]),
        Protein(namespace = "HGNC", name = snp_df['Proteins'][i]),
        relation='hasGeneticVariant',
        citation = "GWAS Central",
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
    
    efo_id = input('Please enter the index value(s) of your disease(s) of interest. For multiple, separate by a space. Input: ')
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
    
    print('Now fetching real-time data from databases. Be patient!')
    print('\n')
    
    temp_id = efo_id.split(' ')
    #print(temp_id)
    temp_id = [int(x) for x in temp_id]
    #print(temp_id)
    
    drugs_df = pd.DataFrame()
    dis2prot_df = pd.DataFrame()
    dis2snp_df = pd.DataFrame()
    
    for id in temp_id:
        
        
        chembl_list = GetDiseaseAssociatedDrugs(doid['id'][id],ct_phase)
        drugs_df = pd.concat([drugs_df,chembl_list])
        
        prot_list = GetDiseaseAssociatedProteins(doid['id'][id])
        dis2prot_df = pd.concat([dis2prot_df,prot_list])
        
        snp_dgnet = GetDiseaseSNPs(doid['id'][id]) 
        dis2snp_df = pd.concat([dis2snp_df,snp_dgnet])
        
    
    drugs_df = drugs_df.reset_index(drop=True)
    
    dis2prot_df = dis2prot_df.reset_index(drop=True)
    dis2snp_df = dis2snp_df.reset_index(drop=True)
    
    uprot_df = GetDiseaseAssociatedProteinsPlot(dis2prot_df)
    adv_effect = pd.DataFrame()
          
    #create empty KG
    kg = pybel.BELGraph(name=kg_name, version="0.0.1")
    
    uprot_list = list(set(uprot_df['UniProt']))
    
    uprot_ext = ExtractFromUniProt(uprot_list)
    
    if vir_prot:
        vir_uprot_ext = ExtractFromUniProt(vir_prot)

    if not drugs_df.empty:
    
        print('A total of ' + str(len(list(set(drugs_df['drugId'])))) + ' drugs have been identified. Now fetching relevant data')
        
        chembl2mech = RetMech(list(set(drugs_df['drugId'])))
        chembl2act = RetAct(list(set(drugs_df['drugId'])))
        
        prtn_as_chembl = Ret_chembl_protein(chembl2act) + Ret_chembl_protein(chembl2mech)
        prtn_as_chembl = set(prtn_as_chembl)
        prtn_as_chembl = list(prtn_as_chembl)
        chembl2uprot = chembl2uniprot(prtn_as_chembl)
        
        chembl2act = chembl2gene2path(chembl2uprot, chembl2act)
        chembl2mech = chembl2gene2path(chembl2uprot, chembl2mech)
    
        kg = chem2moa_rel(chembl2mech,'HGNC',uprot_list, kg)
        kg = chem2act_rel(chembl2act,'HGNC',uprot_list, kg)
        kg = gene2path_rel(chembl2uprot, 'HGNC',uprot_list, kg)
    
        adv_effect = GetAdverseEvents(list(set(drugs_df['drugId'])))
        kg = chembl2adverseEffect_rel(adv_effect,kg)
        
        kg = chembl_annotation(kg)
        kg = chembl_name_annotation(kg,drugs_df)
        
    
    kg = uniprot_rel(uprot_ext, 'HGNC', kg)
    
    kg = gene_ontology_annotation(kg,uprot_ext)
    kg = protein_annotation_druggability(kg)
    
    if vir_prot:
        kg = uniprot_rel(vir_uprot_ext,'VP',kg)
    
    #snp_dgnet = GetDiseaseSNPs(doid['id'][efo_id])  
    
    # if snp_dgnet != None:
        # kg = snp2gene_rel(snp_dgnet,kg)
    
    if not dis2snp_df.empty:
        kg = snp2gene_rel(dis2snp_df,kg)
    
    print('Your KG is now generated!','\n')
    
    saveFiles(kg_name, uprot_df, adv_effect, kg, drugs_df, dis2snp_df,uprot_ext)
    
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

def GetSmiles(drugs_df,colname):

    #chembl ids can be antibodies which have sequence instead of smiles
    
    drugs=drugs_df
    #drugs = pd.read_csv('data/kgs/metabolic diseases/t1dm/diseaseAssociatedDrugs.csv')
    drugs_list = set(list(drugs[colname]))
    #print(drugs_list)
    molecule = new_client.molecule
    temp_list = []
    
    for item in tqdm(drugs_list,desc='Getting SMILES for CHEMBL ids and generating descriptors'):
        
        try:
            #print(item)
            mol = molecule.filter(chembl_id=str(item)).only(['molecule_chembl_id', 'molecule_structures'])
            mol_smiles = mol[0]['molecule_structures']['canonical_smiles']
            temp = [item,mol_smiles]
            temp_list.append(temp)
        
        except:
            print('This ChEMBL id could not be parsed: ' + str(item))
            #print(item)
            continue
    
    print('Info: Some drug has no SMILES representation because its type is either antibody, protein or unknown')
    
    temp_df = pd.DataFrame(temp_list,columns=['drugId','smiles'])     
       
                              
    return(temp_df)
                          
    
def lipinski_rule_of_5(df):
    
    df_ro5 = []
    
    pass_ro5_binary = []
    for item in df['smiles']:
    
        mol = Chem.MolFromSmiles(item) 

        # Ro5 descriptors

        MW = Descriptors.MolWt(mol)

        HBA = Descriptors.NOCount(mol)

        HBD = Descriptors.NHOHCount(mol)

        LogP = Descriptors.MolLogP(mol)

        conditions = [MW <= 500, HBA <= 10, HBD <= 5, LogP <= 5]
        
        violation = conditions.count(False)
                
        pass_ro5 = conditions.count(True) >= 3
        
        if pass_ro5 == True:
            pass_ro5_binary.append(0)
        
        else:
            pass_ro5_binary.append(1)
        
        temp = [violation]
        
        #flag_list.append(pass_ro5)
        
        descriptors = [MW,HBA,HBD,LogP]
        
        descriptors = descriptors + temp 
        
        df_ro5.append(descriptors)
    
    names = 'MW<500 HBA<10 HBD<5 LogP<5 Violation'.split(' ')
    df_ro5 = pd.DataFrame(df_ro5,columns=names)
    df_ro5['Druglikeness'] = pass_ro5_binary
    
    
    df = pd.concat([df, df_ro5], axis=1) # using ignore_index=True removes colnames
    return(df)
                                  
           
def GetDrugWarnings(chem_list):
    
    api_response = pd.DataFrame()
    
    for chem in tqdm(chem_list, desc = 'Retrieving drug warnings for each drug'):
        
        chembl_id = chem
        
        try:
        
            #get total no. of adverse effects for a given drug
            #count = getAdverseEffectCount(chembl_id)
            
            query_string = """
                query DrugWarningsQuery($chemblId: String!) {
                  drug(chemblId: $chemblId) {
                    id
                    drugWarnings {
                      warningType
                      toxicityClass
                    }
                  }
                }

            """

            # Set variables object of arguments to be passed to endpoint
            variables = {"chemblId": chembl_id}

            # Set base URL of GraphQL API endpoint
            base_url = "https://api.platform.opentargets.org/api/v4/graphql"

            # Perform POST request and check status code of response
            r = requests.post(base_url, json={"query": query_string, "variables": variables})
            #r = requests.post(base_url, json={"query": query_string})
            #print(r.status_code)

            # Transform API response from JSON into Python dictionary and print in console
            api_response_temp = json.loads(r.text)
            
            #return(api_response_temp)

            api_response_temp = api_response_temp['data']['drug']['drugWarnings']
            api_response_temp = pd.DataFrame(api_response_temp)
            api_response_temp ['chembl_id'] = chembl_id

            api_response = pd.concat([api_response,api_response_temp])
            
        except:
            continue
          
    #remove rows with 'None' values in col 'toxicityClass'       
    api_response = api_response[api_response.astype(str).ne('None').all(1)]    
    api_response.reset_index(drop=True, inplace=True)
    return(api_response)
    
def chembl2drugWarnings_rel(
    chembl_drugWarnings_df,
    graph: BELGraph
) -> BELGraph:
    """

    :param chembl_adveff_df:
    :param graph:
    :return:
    """

    for i in range(len(chembl_drugWarnings_df)):

        graph.add_qualified_edge(
            Abundance(namespace='ChEMBL', name= str(chembl_drugWarnings_df['chembl_id'][i])),
            Pathology(namespace='DrugWarning', name= str(chembl_drugWarnings_df['toxicityClass'][i])),  
            relation='hasDrugWarning',
            citation="OpenTargets Platform",
            evidence='DrugReactions'
        )

    return graph
    
def getProtfromKG(mainGraph):

    prot_list = []
    for u, v, data in tqdm(mainGraph.edges(data=True),desc='Filtering Proteins/Genes'):
        
        if 'HGNC' in u.namespace:
            if u.name not in prot_list:
                prot_list.append(u.name)

        if 'HGNC' in v.namespace:
            if v.name not in prot_list:
                prot_list.append(v.name)
                
    return(prot_list)
    
def getChemfromKG(mainGraph):

    chem_list = []
    for u, v, data in tqdm(mainGraph.edges(data=True),desc='Filtering Chemicals'):

        if 'CHEMBL' in u.name:
            print(u)
            print(u.name)
            print(u.namespace)
            if u.name not in chem_list:
                chem_list.append(u.name)

        if 'CHEMBL' in v.name:
            if v.name not in chem_list:
                chem_list.append(v.name)
                
    return(chem_list)
    
def getEntityFromKG(mainGraph,namespaceType):
    
    entity_list = []
    
    for u, v, data in tqdm(mainGraph.edges(data=True),desc='Filtering ' + namespaceType):
        
        if namespaceType in u.namespace:
            if u.name not in entity_list:
                entity_list.append(u.name)

        if namespaceType in v.namespace:
            if v.name not in entity_list:
                entity_list.append(v.name)
                
    return(entity_list)

def ro5_filter(df):
    
    temp_df = []
    
    violation_counts_list = []
    
    pass_list = []
    
    for item in df['smiles']:
    
        molecule = Chem.MolFromSmiles(item) 

#         Lipinski Ro5
        
#         Moleculer Weight <= 500
#         LogP <= 5
#         H-Bond Donor Count <= 5
#         H-Bond Acceptor Count <= 10

        molecular_weight = Descriptors.MolWt(molecule)

        h_bond_acceptors = Descriptors.NOCount(molecule)

        h_bond_donor = Descriptors.NHOHCount(molecule)

        logp = Descriptors.MolLogP(molecule)
        
        conditions = [molecular_weight <= 500, h_bond_acceptors <= 10, h_bond_donor <= 5, logp <= 5]
        
        violation_counts = conditions.count(False)
        
        violation_counts_list.append(violation_counts)
        
        descriptors = [molecular_weight,h_bond_acceptors,h_bond_donor,logp]
        
        temp_df.append(descriptors)
                
        pass_ro5 = conditions.count(True) >= 3
        
        if pass_ro5 == True:
            pass_list.append(0)
        
        else:
            pass_list.append(1)
            
    
    names = 'MW HBA HBD LogP'.split(' ')
    temp_df = pd.DataFrame(temp_df, columns=names)
    
    
    df['Violation(s)_ro5'] = violation_counts_list
    df['Lipinski_ro5'] = pass_list
    
    df = pd.concat([df, temp_df], axis=1)
    
    return(df)
    
def ghose_filter(df):

    temp_df = []
    
    pass_list = []
    
    for item in df['smiles']:
    
        molecule = Chem.MolFromSmiles(item) 

        # ghose descriptors
        #     Molecular weight between 160 and 480
        #     LogP between -0.4 and +5.6
        #     Atom count between 20 and 70
        #     Molar refractivity between 40 and 130
        
        molar_refractivity = Chem.Crippen.MolMR(molecule)
        
        number_of_atoms = Chem.rdchem.Mol.GetNumAtoms(molecule)

        molecular_weight = Descriptors.MolWt(molecule)

        logp = Descriptors.MolLogP(molecule)
        
        conditions = [molecular_weight >= 160 and molecular_weight <= 480, 
              number_of_atoms >= 20 and number_of_atoms <= 70, 
              molar_refractivity >= 40 and molar_refractivity <= 130, 
              logp >= -0.4 and logp <= 5.6]
                
        pass_ghose = conditions.count(True) == 4
        
        descriptors = [number_of_atoms,molar_refractivity]
        
        temp_df.append(descriptors)
        
        if pass_ghose == True:
            pass_list.append(0)
        
        else:
            pass_list.append(1)
            
    names = 'AtomNum MolRefractivity'.split(' ')
    temp_df = pd.DataFrame(temp_df, columns=names)

    df['Ghose'] = pass_list
    
    df = pd.concat([df, temp_df], axis=1)
    
    return(df)
    
def veber_filter(df):

    temp_df = []
    
    pass_list = []
    
    for item in df['smiles']:
    
        molecule = Chem.MolFromSmiles(item) 

#         Veber filter
#         Rotatable bonds <= 10
#         Topological polar surface area <= 140

        topological_surface_area_mapping = Chem.QED.properties(molecule).PSA
    
        rotatable_bonds = Descriptors.NumRotatableBonds(molecule)

        conditions = [rotatable_bonds <= 10, 
                      topological_surface_area_mapping <= 140
        ]
                
        pass_veber = conditions.count(True) == 2
        
        descriptors = [rotatable_bonds,topological_surface_area_mapping]
        
        temp_df.append(descriptors)
        
        if pass_veber == True:
            pass_list.append(0)
        
        else:
            pass_list.append(1)
            
    names = 'RotBond TPSA'.split(' ')
    temp_df = pd.DataFrame(temp_df, columns=names)

    df['Veber'] = pass_list
    
    df = pd.concat([df, temp_df], axis=1)
    
    return(df)
    
def reos_filter(df):
    
    temp_df = []
    
    pass_list = []
    
    for item in df['smiles']:
    
        molecule = Chem.MolFromSmiles(item) 

# REOS:
#     Molecular weight between 200 and 500
#     LogP between -5.0 and +5.0
#     H-bond donor count between 0 and 5
#     H-bond acceptor count between 0 and 10
#     Formal charge between -2 and +2
#     Rotatable bond count between 0 and 8
#     Heavy atom count between 15 and 50

        molecular_weight = Descriptors.ExactMolWt(molecule)
        logp = Descriptors.MolLogP(molecule)
        h_bond_donor = Descriptors.NumHDonors(molecule)
        h_bond_acceptors = Descriptors.NumHAcceptors(molecule)
        rotatable_bonds = Descriptors.NumRotatableBonds(molecule)
        formal_charge = Chem.rdmolops.GetFormalCharge(molecule)
        heavy_atoms = Chem.rdchem.Mol.GetNumHeavyAtoms(molecule)
        
        #print(molecular_weight,logp,h_bond_donor,h_bond_acceptors,formal_charge,rotatable_bonds,heavy_atoms)
        
        conditions = [molecular_weight >= 200 and molecular_weight <= 500, 
              logp >= -0.4 and logp <= 5.6,        
              h_bond_donor >= 0 and h_bond_donor <= 5, 
              h_bond_acceptors >= 0 and h_bond_acceptors <= 10,
              formal_charge >= -2 and formal_charge <= 2,
              rotatable_bonds >= 0 and rotatable_bonds <= 8,
              heavy_atoms >= 15 and heavy_atoms <= 50]
              
        descriptors = [formal_charge,heavy_atoms]
        
        temp_df.append(descriptors)
        
        #print(conditions)
        pass_reos = conditions.count(True) == 7
        
        if pass_reos == True:
            pass_list.append(0)
        
        else:
            pass_list.append(1)
            
    names = 'Charge HeavyAtom'.split(' ')
    temp_df = pd.DataFrame(temp_df, columns=names)

    df['REOS'] = pass_list
    
    df = pd.concat([df, temp_df], axis=1)
    
    return(df)
    
def qed_filter(df):
    
    temp_df = []
    
    pass_list = []
    
    for item in df['smiles']:
    
        molecule = Chem.MolFromSmiles(item) 

# Drug-Like (QED):
#     mass < 400
#     ring count > 0
#     rotatable bond count < 5
#     h-bond donor count <= 5
#     h-bond acceptor count <= 10
#     logP < 5
    
        molecular_weight = Descriptors.ExactMolWt(molecule)
        logp = Descriptors.MolLogP(molecule)
        h_bond_donor = Descriptors.NumHDonors(molecule)
        h_bond_acceptors = Descriptors.NumHAcceptors(molecule)
        rotatable_bonds = Descriptors.NumRotatableBonds(molecule)
        num_of_rings = Chem.rdMolDescriptors.CalcNumRings(molecule)
      
        #print(molecular_weight,logp,h_bond_donor,h_bond_acceptors,formal_charge,rotatable_bonds,heavy_atoms)
        
        conditions = [molecular_weight <= 400, 
              logp <= 5,        
              h_bond_donor <= 5, 
              h_bond_acceptors <= 10,
              num_of_rings > 0,
              rotatable_bonds <= 5]
              
        #descriptors = [num_of_rings]
        
        temp_df.append(num_of_rings)
        
        #print(conditions)
        pass_qed = conditions.count(True) == 6
        
        if pass_qed == True:
            pass_list.append(0)
        
        else:
            pass_list.append(1)

    df['QED'] = pass_list
    df['RingNum'] = temp_df
    
    return(df)

def calculate_filters(df,colname_chembl):
    
    df_smiles = GetSmiles(df,colname_chembl)
    
    #df_smiles = remove_salt(df_smiles,'smiles')
    
    temp = ro5_filter(df_smiles)
    temp = ghose_filter(temp)
    temp = veber_filter(temp)
    temp = reos_filter(temp)
    temp = qed_filter(temp)
    
    df = df.loc[df.groupby('drugId')['phase'].idxmax()]
    df = df.reset_index(drop=True)
    
    df = pd.merge(temp,df[['drugId','phase']],on='drugId', how='left')
    
    df[['phase']] = df[['phase']].astype(int)
    
    df = df.drop_duplicates()
    df = df.reset_index(drop=True)
    
    return(df)

def prot2crossRefDB(prot_list):
    
    import requests
    
    final_list = []
    
    for prot in prot_list:
        
        print(prot)
        temp_list = []
        
        temp_list.append(prot)
        
        files = {
        'from': (None, 'Gene_Name'),
        'to': (None, 'UniProtKB-Swiss-Prot'),
        'ids': (None, prot),
        'taxId': (None, '9606'),
        }
    
        response = requests.post('https://rest.uniprot.org/idmapping/run', files=files)
        
        response = response.json()
        
        job_id = response['jobId']
        
#         try:
#             output = get_id_mapping_results_link(job_id)
            
#         except: 
#             output = get_id_mapping_results_link(job_id)
        
        API_URL = "https://rest.uniprot.org"
        session = requests.Session()
        url = f"{API_URL}/idmapping/details/{job_id}"
        request = session.get(url)
        check_response(request)
        print(check_response)
        print(request.json())

        output= request.json()["redirectURL"]
        

        response = requests.get(output)
        
        response = response.json()
        
        if response['results']: 
        
        #print(response)
        
            uprot = response['results'][0]['to']['primaryAccession']

            temp_list.append(uprot)

            uniProtkbId = response['results'][0]['to']['uniProtkbId'] 
            organism = response['results'][0]['to']['organism']['scientificName']

            #get all crossRef DBs
            crossRefs = response['results'][0]['to']['uniProtKBCrossReferences']

            db_list = []
            for item in crossRefs:

                db_list.append(item['database'])


            if 'OpenTargets' in db_list:
                getIndex = db_list.index('OpenTargets')
                temp_list.append(crossRefs[getIndex]['id'])

            else:
                temp_list.append('None')


            if 'ChEMBL' in db_list:
                getIndex = db_list.index('ChEMBL')
                temp_list.append(crossRefs[getIndex]['id'])

            else:
                temp_list.append('None')


            if 'HGNC' in db_list:
                getIndex = db_list.index('HGNC')
                temp_list.append(crossRefs[getIndex]['id'])

            else:
                temp_list.append('None')

            temp_list.append(uniProtkbId)
            temp_list.append(organism)

            #print(temp_list)

            final_list.append(temp_list)
        
        #print(final_list)
        
    cols = ['Protein','UniProt','ENSG','ChEMBL','HGNC','uniProtkbId','Organism']
    
    df= pd.DataFrame(final_list,columns=cols)             
        
    return(df)
    
def druglikeness_parallelCoordinates(df):
    
    df = df_final.loc[df_final['MW'] <= 2000]
    
    print('This plot excludes drugs with molecular weight > 2000')

    fig = go.Figure(data=
        go.Parcoords(
            line = dict(color = df['phase'],showscale = True,
                       colorscale = [[0,'purple'],[0.5,'lightseagreen'],[1,'gold']]),

            dimensions = list([

                dict(range = [min(df['MW']),max(df['MW'])],
                    #constraintrange = [4,8],
                    label = 'MW', values = df['MW']),

                dict(range = [min(df['TPSA']),max(df['TPSA'])],
                    #constraintrange = [4,8],
                    label = 'TPSA', values = df['TPSA']),

                dict(range = [min(df['MolRefractivity']),max(df['MolRefractivity'])],
                    #constraintrange = [4,8],
                    label = 'MolRefractivity', values = df['MolRefractivity']),

                dict(range = [min(df['LogP']),max(df['LogP'])],
                    #constraintrange = [4,8],
                    label = 'LogP', values = df['LogP']),

                dict(range = [min(df['AtomNum']),max(df['AtomNum'])],
                    #constraintrange = [4,8],
                    label = 'AtomNum', values = df['AtomNum']),

                dict(range = [min(df['RotBond']),max(df['RotBond'])],
                    #constraintrange = [4,8],
                    label = 'RotBond', values = df['RotBond']),

                dict(range = [min(df['HeavyAtom']),max(df['HeavyAtom'])],
                    #constraintrange = [4,8],
                    label = 'HeavyAtom', values = df['HeavyAtom']),

                dict(range = [min(df['RingNum']),max(df['RingNum'])],
                    #constraintrange = [4,8],
                    label = 'RingNum', values = df['RingNum']),

                dict(range = [0,4],
                    tickvals = [0,1,2,3,4],
                    label = 'Violation', values = df['Violation(s)_ro5']),

                dict(range = [-1,2],
                    tickvals = [0,1],
                    ticktext = ['Yes','No'],
                    label = 'Lipinski_ro5', values = df['Lipinski_ro5']),

                dict(range = [-1,2],
                    tickvals = [0,1],
                    ticktext = ['Yes','No'],
                    label = 'Ghose', values = df['Ghose']),

                dict(range = [-1,2],
                    tickvals = [0,1],
                    ticktext = ['Yes','No'],
                    label = 'Veber', values = df['Veber']),

                dict(range = [-1,2],
                    tickvals = [0,1],
                    ticktext = ['Yes','No'],
                    label = 'REOS', values = df['REOS']),

                dict(range = [-1,2],
                    tickvals = [0,1],
                    ticktext = ['Yes','No'],
                    label = 'QED', values = df['QED']),

                dict(range = [1,4],
                    tickvals = [1,2,3,4],
                    label = 'Phase', values = df['phase']),

            ])
        )
    )

    fig.update_layout(
        plot_bgcolor = 'white',
        paper_bgcolor = 'white'
    )

    fig.show()

def druglikeness(df,filename):
    
    #select cols with druglikness flags as 0 (Yes) and 1 (No)
    filter_cols = ['Lipinski_ro5','Ghose','Veber','REOS','QED']

    df = df[filter_cols]
    df = df.stack().groupby(level=[1]).value_counts().unstack()
    
    # Preferred order of stacked bar elements
    stack_order = [0, 1]
    df = df[stack_order]

    ax = df.plot.bar(rot=0, stacked=True,
                      figsize = [8,8],
                      title='Bar plots showing drug-likeness across various filters',
                      xlabel = 'Filters for drug-likeness',
                      ylabel = 'Total number of drugs',
                      #table = df,
                      linewidth = 4,
                      edgecolor='white', 
                      #grid = 2
                     )

    ax.set_facecolor("ghostwhite")
    #%matplotlib inline

    ax.set_axisbelow(True)
    ax.yaxis.grid(color='white', linewidth = 2)
    ax.legend(['Yes', "No"])
    
    path = '../streamlit_ITMP/images/' + filename + '.png'
    
    ax.get_figure().savefig(path, format='png', dpi = 600)
    

def removeSalt(mol):
    remover = SaltRemover()
    stripped = remover.StripMol(mol)
    smiles = Chem.MolToSmiles(stripped)

    if '.' in smiles:
        temp = smiles.split(".")
        smiles = max(temp, key=len)
        return(smiles)

    else:
        return(smiles)


    