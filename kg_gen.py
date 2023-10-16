# -*- coding: utf-8 -*-

"""Utils files with all functions relevant to generation of KG."""

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
from pybel.dsl import Protein, Abundance, Pathology, BiologicalProcess
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

def GetQuery():
    
    query = input('Please enter the disease you are interested in and we will try to find the best matches for you.' +'\n' + '\n' + 'Input: ')
    return(query)

def Generate_KG():
        
    from IPython.display import Image    
    printmd("**Welcome to the KG Generator tool. In the following steps, we will need a couple of inputs from your side.**",color = "blue")
   
    #Image(filename='data/KGG.png')
    
    #this is required to get above printed before the input is asked
    time.sleep(0.05)
    
    #dis = input('Please enter the disease you are interested in and we will try to find the best matches for you.' +'\n' + '\n' + 'Input: ')
    
    dis = GetQuery()
    
    temp = searchDisease(dis)
    #print(temp)
    #print('\n')
    if not temp.empty:
        printmd('**Here you go! Hopefully your disease of interest is in the list. If so, let\'s get started.**')
        #print('\n')
        #print(temp)
        return(temp)
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
        print('Your dataframe is ready')
        return(df)
    
    else:
        print('No drugs found in clinical trials')
    
    #chembl_list = list(set(df['drugId']))
    #return(chembl_list)
    
def KG_namespace_plot(final_kg):
    
    import matplotlib.pyplot as plt
    nspace_count = pybel.struct.summary.count_namespaces(final_kg)
    nspace_count = dict(nspace_count)

    nspace_data = {'Namespace':list(nspace_count.keys()),'Number':list(nspace_count.values())}
    nspace = pd.DataFrame(nspace_data)
    plt.figure()

    a = sns.barplot(x="Number", y="Namespace", data=nspace_data)
    a.set(xlabel='Number',ylabel='Namespace',title= 'KG Namespace in numbers')

    plt.tight_layout()
    #plt.savefig('data/export/test2.png',dpi=600)
    plt.show()    

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
    
def createKG():

    #import matplotlib.pyplot as plt
    #import matplotlib.image as mpimg
    image = mpimg.imread("data/KGG.png")
    plt.imshow(image)
    plt.axis('off')
    plt.show()
    
    doid = Generate_KG()
    #print(doid)
    
    time.sleep(0.1)
    
    efo_id = int(input('Please enter the index value of your disease of interest. Input: '))
    #print(efo_id)
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
    
    chembl_list = GetDiseaseAssociatedDrugs(doid['id'][efo_id],ct_phase)
    
    return (chembl_list)
       
    #create empty KG
    kg = pybel.BELGraph(name=kg_name, version="0.0.1")
    
    uprot_ext = ExtractFromUniProt(df_dis2prot['UniProt'])
    
    print('A total of ' + str(len(chembl_list)) + ' drugs have been identified. Now fetching relevant data')
    
    chembl2mech = RetMech(chembl_list)
    chembl2act = RetAct(chembl_list)
    
    prtn_as_chembl = Ret_chembl_protein(chembl2act) + Ret_chembl_protein(chembl2mech)
    prtn_as_chembl = set(prtn_as_chembl)
    prtn_as_chembl = list(prtn_as_chembl)
    chembl2uprot = chembl2uniprot(prtn_as_chembl)
    
    chembl2act = chembl2gene2path(chembl2uprot, chembl2act)
    chembl2mech = chembl2gene2path(chembl2uprot, chembl2mech)
    
    
    kg = uniprot_rel(uprot_ext, 'HGNC', kg)
    
    kg = chem2moa_rel(chembl2mech, 'HGNC', kg)
    kg = chem2act_rel(chembl2act, 'HGNC', kg)
    kg = gene2path_rel(chembl2uprot, 'HGNC', kg)
    
    adv_effect = GetAdverseEvents(chembl_list)
    kg = chembl2adverseEffect_rel(adv_effect,kg)
    
    print('Your KG is now generated!')
    
    return(kg)
    
    
           
