# -*- coding: utf-8 -*-

"""Utils files with all functions relevant to generation of KG."""

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
from pybel.dsl import Protein, Abundance, Pathology, BiologicalProcess, Gene, MicroRna, Population 
from tqdm.auto import tqdm
import time
import seaborn as sns
import pybel

logger = logging.getLogger("__name__")


def RetMech(chemblIds) -> dict:
    """Function to retrieve mechanism of actions and target proteins from ChEMBL

    :param chemblIds:
    :return:
    """
    getMech = new_client.mechanism

    mechList = []
    for chemblid in tqdm(chemblIds, desc='Retrieving mechanisms from ChEMBL'):
        mechs = getMech.filter(
            molecule_chembl_id=chemblid
        ).only(['mechanism_of_action', 'target_chembl_id','action_type'])
        
        mechList.append(list(mechs))

    named_mechList = dict(zip(chemblIds, mechList))
    named_mechList = {
        k: v
        for k, v in named_mechList.items()
        if v
    }
    return named_mechList


def RetDrugInd(chemblIDs) -> dict:
    """Function to retrieve associated diseases from ChEMBL

    :param chemblIDs:
    :return:
    """
    getDrugInd = new_client.drug_indication

    drugIndList = []
    for chemblid in tqdm(chemblIDs, desc='Retrieving diseases from ChEMBL'):
        drugInd = getDrugInd.filter(
            molecule_chembl_id=chemblid
        ).only('mesh_heading')
        drugIndList.append(list(drugInd))

    named_drugIndList = dict(zip(chemblIDs, drugIndList))
    named_drugIndList = {
        k: v
        for k, v in named_drugIndList.items()
        if v
    }
    return named_drugIndList
    
def RetAct_pchem_5(chemblIds) -> dict:
    """Function to retrieve associated assays from ChEMBL

    :param chemblIds:
    :return:
    """
    GetAct = new_client.activity
    getTar = new_client.target
    ActList = []
    filtered_list=['assay_chembl_id','assay_type','pchembl_value','target_chembl_id',
                   'target_organism','bao_label','target_type']

#     filtered_list = [
#         'pchembl_value',
#         'target_chembl_id',
#         'target_type',
#         'bao_label'
#     ]

    for chembl in tqdm(chemblIds, desc='Retrieving bioassays from ChEMBL'):
    #for i in range(len(chemblIds)):
        acts = GetAct.filter(
            molecule_chembl_id=chembl,
            pchembl_value__isnull=False,
            assay_type_iregex='(B|F)',
            target_organism='Homo sapiens'
        ).only(filtered_list)
        
        #print(chemblIds[i])
        data = []
        #print(acts)
        
        for d in acts:

            if float(d.get('pchembl_value')) <= 5:
                continue
                            
            if (d.get('bao_label') != 'single protein format'):
                continue
                
            tar = d.get('target_chembl_id')   
            tar_dict = getTar.get(tar)
            #print(tar_dict)
            
            try:
                if tar_dict['target_type'] in ('CELL-LINE', 'UNCHECKED'):
                    continue
            except KeyError:
                continue
                

            #uprot_id = d['target_components'][0]['accession']
            #print(uprot_id)
            #print(d)
            data.append(d)

        # acts = acts[:5]
        ActList.append(list(data))

    named_ActList = dict(zip(chemblIds, ActList))
    named_ActList = {
        k: v
        for k, v in named_ActList.items()
        if v
    }
    return named_ActList

def RetAct(chemblIds) -> dict:
    """Function to retrieve associated assays from ChEMBL

    :param chemblIds:
    :return:
    """
    GetAct = new_client.activity
    getTar = new_client.target
    ActList = []
    filtered_list=['assay_chembl_id','assay_type','pchembl_value','target_chembl_id',
                   'target_organism','bao_label','target_type']

#     filtered_list = [
#         'pchembl_value',
#         'target_chembl_id',
#         'target_type',
#         'bao_label'
#     ]

    for chembl in tqdm(chemblIds, desc='Retrieving bioassays from ChEMBL'):
    #for i in range(len(chemblIds)):
        acts = GetAct.filter(
            molecule_chembl_id=chembl,
            pchembl_value__isnull=False,
            assay_type_iregex='(B|F)',
            target_organism='Homo sapiens'
        ).only(filtered_list)
        
        #print(chemblIds[i])
        data = []
        #print(acts)
        
        for d in acts:

            if float(d.get('pchembl_value')) < 6:
                continue
                            
            if (d.get('bao_label') != 'single protein format'):
                continue
                
            tar = d.get('target_chembl_id')   
            tar_dict = getTar.get(tar)
            #print(tar_dict)
            
            try:
                if tar_dict['target_type'] in ('CELL-LINE', 'UNCHECKED'):
                    continue
            except KeyError:
                continue
                

            #uprot_id = d['target_components'][0]['accession']
            #print(uprot_id)
            #print(d)
            data.append(d)

        # acts = acts[:5]
        ActList.append(list(data))

    named_ActList = dict(zip(chemblIds, ActList))
    named_ActList = {
        k: v
        for k, v in named_ActList.items()
        if v
    }
    return named_ActList

# def RetAct(chemblIds) -> dict:
    # """Function to retrieve associated assays from ChEMBL

    # :param chemblIds:
    # :return:
    # """
    # GetAct = new_client.activity
    # ActList = []

    # filtered_list = [
        # 'pchembl_value',
        # 'target_chembl_id',
        # 'target_type'
    # ]

    # for chembl in tqdm(chemblIds, desc='Retrieving bioassays from ChEMBL'):
        # acts = GetAct.filter(
            # molecule_chembl_id=chembl,
            # pchembl_value__isnull=False,
            # assay_type_iregex='(B|F)',
            # target_organism='Homo sapiens'
        # ).only(filtered_list)

        # data = []

        # for d in acts:

            # if float(d.get('pchembl_value')) < 6:
                # continue
            # try:

                # if d['target_type'] in ('CELL-LINE', 'UNCHECKED'):
                    # continue
            # except IndexError:
                # continue

            # uprot_id = d['target_components'][0]['accession']

            # data.append(uprot_id)

        # # acts = acts[:5]
        # ActList.append(list(data))

    # named_ActList = dict(zip(chemblIds, ActList))
    # named_ActList = {
        # k: v
        # for k, v in named_ActList.items()
        # if v
    # }
    # return named_ActList


def Ret_chembl_protein(sourceList) -> list:
    """Method to retrieve ChEMBL ids which are proteins/targets

    :param sourceList:
    :return:
    """
    protein_List = []
    for item in sourceList:
        for j in range(len(sourceList[item])):
            protein_List.append(sourceList[item][j]['target_chembl_id'])

    protein_List = set(protein_List)
    protein_List = list(filter(None, protein_List))
    return protein_List


def chembl2uniprot(chemblIDs) -> dict:
    """Method to convert ChEMBL id to UNIPROT and get associated REACTOME pathways

    :param chemblIDs:
    :return:
    """
    getTarget = new_client.target
    chem2Gene2path = []
    chemHasNoPath = set()
    chemNotprotein = set()

    chem2path = defaultdict(list)

    # Loop to ensure it is a protein
    for chemblid in tqdm(chemblIDs, desc='Filtering UniProt proteins from ChEMBL'):
        chem = getTarget.filter(
            chembl_id=chemblid
        ).only('target_components')

        try:
            uprot_id = chem[0]['target_components'][0]['accession']

            if not uprot_id:
                chemHasNoPath.add(chemblid)

        except IndexError:
            chemHasNoPath.add(chemblid)

    logger.info(f'No UniProt information available for {len(chemHasNoPath)} proteins.')

    chemblIDs_filtered = [
        item
        for item in chemblIDs
        if item not in chemHasNoPath
    ]

    # Get gene symbol from ChEMBL and filtering the list for human proteins only
    for chemblid in tqdm(chemblIDs_filtered, desc='Filtering human proteins from ChEMBL'):

        chem = getTarget.filter(chembl_id=chemblid).only('target_components')
        getGene = chem[0]['target_components'][0]['target_component_synonyms']
        try:
            getGene = [item for item in getGene if item["syn_type"] == "GENE_SYMBOL"][0]

            if not getGene:
                chemNotprotein.add(chemblid)

        except IndexError:
            chemNotprotein.add(chemblid)

    chemblIDs_filtered = [
        item
        for item in chemblIDs_filtered
        if item not in chemNotprotein
    ]

    # Extracting data for valid proteins only
    for chemblid in tqdm(chemblIDs_filtered, desc='Populating ChEMBL data for human proteins'):
        chem = getTarget.filter(
            chembl_id=chemblid
        ).only('target_components')

        # UniProt data
        uprot_id = chem[0]['target_components'][0]['accession']

        # Gene symbol
        getGene = chem[0]['target_components'][0]['target_component_synonyms']
        getGene = [item for item in getGene if item["syn_type"] == "GENE_SYMBOL"][0]

        # Pathway data
        chem2path = [item for item in chem[0]['target_components'][0]['target_component_xrefs'] if
                     item["xref_src_db"] == "Reactome"]

        uprot = {'accession': uprot_id}
        chem2path.append(uprot)
        chem2path.append(getGene)
        chem2Gene2path.append(chem2path)

    named_chem2Gene2path = dict(zip(chemblIDs_filtered, chem2Gene2path))
    named_chem2Gene2path = {
        k: v
        for k, v in named_chem2Gene2path.items()
        if v
    }
    return named_chem2Gene2path


def chembl2gene2path(
    chem2geneList,
    ActList
):
    """Method for updating chembl protein nodes with gene symbol.

    :param chem2geneList:
    :param ActList:
    :return:
    """
    for item in chem2geneList:
        sizeOfitem = len(chem2geneList[item])
        gene = chem2geneList[item][sizeOfitem - 1]['component_synonym']
        #uniprot id is located in 2nd last pos
        uprot = chem2geneList[item][-2]['accession']
        for jtem in ActList:
            for i in range(len(ActList[jtem])):
                if item == ActList.get(jtem)[i]['target_chembl_id']:
                    newkey = {'Protein': gene,'Accession':uprot}
                    ActList[jtem][i].update(newkey)

    return ActList


def Ret_uprotid(chembl2uprot) -> list:
    """Method to get UniProt ids from dict of chembl2uniprot

    :param chembl2uprot:
    :return:
    """
    chemblProt = []
    for item in chembl2uprot:
        for j in range(len(chembl2uprot[item])):
            if 'accession' in chembl2uprot[item][j]:
                chemblProt.append(
                    chembl2uprot[item][j]['accession']
                )
    return chemblProt


def ExtractFromUniProt(uniprot_id) -> dict:
    """Uniprot parser to retrieve information about OMIM disease, reactome pathway, biological process,
     and molecular functions.

    :param uniprot_id:
    :return:
    """
    Uniprot_Dict = []

    mapped_uprot = []

    for id in tqdm(uniprot_id, desc='Fetching Protein-related info'):

        # Retrieve data for id in text format if found in uniprot
        ret_uprot = requests.get(
            'https://www.uniprot.org/uniprot/' + id + '.txt'
        ).text.split('\n')

        if ret_uprot == ['']:
            continue

        id_copy = id
        mapped_uprot.append(id_copy)
        
        #counters for bf and bp
        #i = 0
        #j = 0
        k = 0
        id = {}
        id['Disease'] = {}
        id['Reactome'] = {}
        id['Function'] = {}
        id['BioProcess'] = {}
        id['Gene'] = {}

        # parse each line looking for info about disease, pathway, funcn, bp and so on
        for line in ret_uprot:
            # parse lines with disease and extract disease names and omim ids
            if '-!- DISEASE:' in line:
                if ('[MIM:' in line):
                    dis = line.split(':')
                    id['Disease'].update({dis[1][1:-5]: dis[2][:-1]})

            # extract reactome ids and names
            if 'Reactome;' in line:
                ract = line.split(';')
                id['Reactome'].update({ract[2][1:-1]: ract[1][1:]}) #previously was ract[2][1:-2]

            # look for functions
            if ' F:' in line:
                #if j < 5:
                fn = line.split(';')
                id['Function'].update({fn[2][3:]: fn[1][1:]})
                #j += 1

            # look for biological processes
            if ' P:' in line and 'GO;' in line:
                #if i < 5:
                bp = line.split(';')
                # bp returns list with GO ids and names
                id['BioProcess'].update({bp[2][3:]: bp[1][1:]})
                #i += 1

            if 'GN   Name' in line:
                if k == 0:
                    gene = line.split('=')
                    gene = gene[1].split(' ')
                    if ';' in gene[0]:
                        gene = gene[0].split(';')
                        gene = {'Gene': gene[0]}
                    else:
                        gene = {'Gene': gene[0]}
                    id.update(gene)
                    k += 1

        Uniprot_Dict.append(id)

    Uniprot_Dict = dict(zip(mapped_uprot, Uniprot_Dict))

    return Uniprot_Dict


def chem2moa_rel(
    named_mechList,
    org,otp_prots,
    graph: BELGraph
) -> BELGraph:
    """Method to create the monkeypox graph

    :param named_mechList:
    :param org:
    :param graph: BEL graph of Monkeypox
    :return:
    """
    
    #identified types of chemical and protein action types
    #['INHIBITOR','NEGATIVE ALLOSTERIC MODULATOR','POSITIVE ALLOSTERIC MODULATOR','ANTAGONIST','AGONIST','MODULATOR','BLOCKER','ACTIVATOR','DISRUPTING AGENT', 'SUBSTRATE', 'OPENER','PARTIAL AGONIST','SEQUESTERING AGENT']
    # following lists are used to determine type of edge relationships
    pos = ['POSITIVE ALLOSTERIC MODULATOR','AGONIST','ACTIVATOR','PARTIAL AGONIST']
    neg = ['INHIBITOR','NEGATIVE ALLOSTERIC MODULATOR','ANTAGONIST','BLOCKER']
    misc = ['MODULATOR','DISRUPTING AGENT','SUBSTRATE','OPENER','SEQUESTERING AGENT']

    for chembl_name, chembl_entries in tqdm(named_mechList.items(), desc='Populating Chemical-MoA edges'):
        for info in chembl_entries:
            graph.add_qualified_edge(
                Abundance(namespace='ChEMBL', name=chembl_name),
                BiologicalProcess(namespace='MOA', name=info['mechanism_of_action']),
                relation='hasMechanismOfAction',
                citation='ChEMBL database',
                evidence='ChEMBL query'
            )

            if not info['target_chembl_id']:
                continue

            if 'Protein' in info and info['Accession'] in otp_prots:
                if info['action_type'] in pos:
                    graph.add_increases(
                        Abundance(namespace='ChEMBL', name=chembl_name),
                        Protein(namespace=org, name=info['Protein']),
                        citation='ChEMBL database',
                        evidence='ChEMBL query'
                    )
                if info['action_type'] in neg and info['Accession'] in otp_prots:
                    graph.add_decreases(
                        Abundance(namespace='ChEMBL', name=chembl_name),
                        Protein(namespace=org, name=info['Protein']),
                        citation='ChEMBL database',
                        evidence='ChEMBL query'
                    )

                if info['action_type'] in misc and info['Accession'] in otp_prots:
                    graph.add_qualified_edge(
                        Abundance(namespace='ChEMBL', name=chembl_name),
                        Protein(namespace=org, name=info['Protein']),
                        relation='targets',
                        citation='ChEMBL database',
                        evidence='ChEMBL query'
                    )

            # else:
                # graph.add_association(
                    # Abundance(namespace='ChEMBL', name=chembl_name),
                    # Protein(namespace=org, name=info['target_chembl_id']),
                    # citation='ChEMBL database',
                    # evidence='ChEMBL query'
                # )

    return graph


def chem2dis_rel(
    named_drugIndList,
    graph: BELGraph
) -> BELGraph:
    """Method to add drug indication edges to the KG.

    :param named_drugIndList:
    :param graph:
    :return:
    """
    for chembl_id, drug_entries in tqdm(named_drugIndList.items(), desc='Populating Drug-Indication edges'):
        for drug_data in drug_entries:
            graph.add_qualified_edge(
                Abundance(namespace='ChEMBL', name=chembl_id),
                Pathology(namespace='Disease', name=drug_data['mesh_heading']),
                relation='isUsedIn',
                citation='ChEMBL database',
                evidence='ChEMBL query'
            )
    return graph


def chem2act_rel(
    named_ActList,
    org,otp_prots,
    graph: BELGraph
) -> BELGraph:
    """Method to add bioassay edges to the KG.

    :param named_ActList:
    :param org:
    :param graph:
    :return:
    """
    for chemical, chem_entries in tqdm(named_ActList.items(), desc='Adding bioassay edges to BEL'):
        for chem_data in chem_entries:
            if chem_data['target_chembl_id']:
                if 'Protein' in chem_data and chem_data['Accession'] in otp_prots:
                    graph.add_qualified_edge(
                        Abundance(namespace='ChEMBLAssay', name=chem_data['assay_chembl_id']),
                        Protein(namespace=org, name=chem_data['Protein']),
                        relation='hasTarget',
                        citation='ChEMBL database',
                        evidence='ChEMBL query'
                    )
                # else:
                    # graph.add_association(
                        # Abundance(namespace='ChEMBLAssay', name=chem_data['assay_chembl_id']),
                        # Protein(namespace=org, name=chem_data['target_chembl_id']),
                        # citation='ChEMBL database',
                        # evidence='ChEMBL query'
                    # )

            graph.add_qualified_edge(
                Abundance(namespace='ChEMBL', name=chemical),
                Abundance(namespace='ChEMBLAssay', name=chem_data['assay_chembl_id']),
                relation='hasAssay',
                citation='ChEMBL database',
                evidence='ChEMBL query',
                annotation={
                    'assayType': chem_data['assay_type'],
                    'pChEMBL': chem_data['pchembl_value']
                }
            )

    return graph


def gene2path_rel(
    named_chem2geneList,
    org,otp_prots,
    graph
) -> BELGraph:
    """Method to add protein and reactome data to KG

    :param named_chem2geneList:
    :param org:
    :param graph:
    :return:
    """
    for item in named_chem2geneList:
        itemLen = len(named_chem2geneList[item]) - 1
        for j in range(itemLen - 1):
            
            #checks if uprot id is in otp_proteins
            if named_chem2geneList[item][itemLen - 1]['accession'] in otp_prots:
            
                graph.add_qualified_edge(
                    Protein(namespace=org, name=named_chem2geneList[item][itemLen]['component_synonym']),
                    BiologicalProcess(namespace='Reactome', name=named_chem2geneList[item][j]['xref_name']),
                    relation='hasPathway',
                    citation='ChEMBL database',
                    evidence='ChEMBL query',
                    annotation={
                        'Reactome': 'https://reactome.org/content/detail/'+named_chem2geneList[item][j]['xref_id']
                    }
                )

    return graph




# def uniprot_rel(
    # named_uprotList,
    # org,
    # graph
# ) -> BELGraph:
    # """Method to add UniProt related edges

    # :param named_uprotList:
    # :param org:
    # :param graph:
    # :return:
    # """
    # for item in named_uprotList:
        # fun = list(named_uprotList[item]['Function'].keys())
        # bp = list(named_uprotList[item]['BioProcess'].keys())
        # for f in fun:
            # if str(named_uprotList[item]['Gene']) != 'nan' and not isinstance(named_uprotList[item]['Gene'], dict):
                # graph.add_association(
                    # Protein(namespace=org, name=named_uprotList[item]['Gene']),
                    # BiologicalProcess(namespace='GOMF', name=f),
                    # citation='UniProt database',
                    # evidence='UniProt query'
                # )
            # else:
                # graph.add_association(
                    # Protein(namespace=org, name=item),
                    # BiologicalProcess(namespace='GOMF', name=f),
                    # citation='UniProt database',
                    # evidence='UniProt query'
                # )

        # for b in bp:
            # if str(named_uprotList[item]['Gene']) != 'nan' and not isinstance(named_uprotList[item]['Gene'], dict):
                # graph.add_association(
                    # Protein(namespace=org, name=named_uprotList[item]['Gene']),
                    # BiologicalProcess(namespace='GOBP', name=b),
                    # citation='UniProt database',
                    # evidence='UniProt query'
                # )
            # else:
                # graph.add_association(
                    # Protein(namespace=org, name=item),
                    # BiologicalProcess(namespace='GOBP', name=b),
                    # citation='UniProt database',
                    # evidence='UniProt query'
                # )

    # return graph
    
def uniprot_rel(
    named_uprotList,
    org,
    graph
) -> BELGraph:
    """Method to add UniProt related edges

    :param named_uprotList:
    :param org:
    :param graph:
    :return:
    """
    for item in tqdm(named_uprotList,desc='Populating Uniprot edges'):
        fun = list(named_uprotList[item]['Function'].keys())
        bp = list(named_uprotList[item]['BioProcess'].keys())
        pathway = list(named_uprotList[item]['Reactome'].keys())
        for f in fun:
            if str(named_uprotList[item]['Gene']) != 'nan' and not isinstance(named_uprotList[item]['Gene'], dict):
                graph.add_qualified_edge(
                    Protein(namespace=org, name=named_uprotList[item]['Gene']),
                    BiologicalProcess(namespace='GOMF', name=f),
                    relation='hasMolecularFunction',
                    citation='UniProt database',
                    evidence='UniProt query'
                )
            else:
                graph.add_qualified_edge(
                    Protein(namespace=org, name=item),
                    BiologicalProcess(namespace='GOMF', name=f),
                    relation='hasMolecularFunction',
                    citation='UniProt database',
                    evidence='UniProt query'
                )

        for b in bp:
            if str(named_uprotList[item]['Gene']) != 'nan' and not isinstance(named_uprotList[item]['Gene'], dict):
                graph.add_qualified_edge(
                    Protein(namespace=org, name=named_uprotList[item]['Gene']),
                    BiologicalProcess(namespace='GOBP', name=b),
                    relation='hasBiologicalProcess',
                    citation='UniProt database',
                    evidence='UniProt query'
                )
            else:
                graph.add_qualified_edge(
                    Protein(namespace=org, name=item),
                    BiologicalProcess(namespace='GOBP', name=b),
                    relation='hasBiologicalProcess',
                    citation='UniProt database',
                    evidence='UniProt query'
                )
               
        
        if str(named_uprotList[item]['Gene']) != 'nan' and not isinstance(named_uprotList[item]['Gene'], dict):
            nx.set_node_attributes(graph,{Protein(namespace=org, name=named_uprotList[item]['Gene']):'https://3dbionotes.cnb.csic.es/?queryId='+item},'3Dbio')
            
            nx.set_node_attributes(graph,{Protein(namespace=org, name=named_uprotList[item]['Gene']):'https://www.uniprot.org/uniprotkb/'+item},'UniProt')
        
        else:
            nx.set_node_attributes(graph,{Protein(namespace=org, name=item):'https://3dbionotes.cnb.csic.es/?queryId='+item},'3Dbio')
            
            nx.set_node_attributes(graph,{Protein(namespace=org, name=item):'https://www.uniprot.org/uniprotkb/'+item},'UniProt')
        
    return graph


def _get_target_data(
    protein_list: list,
    organism: str
) -> pd.DataFrame:
    """Get chemical for target data from ChEMBL.

    :param protein_list:
    :param organism:
    :return:
    """
    df_data = []

    target = new_client.target
    activity = new_client.activity

    for protein in tqdm(protein_list, desc='Retrieving chemicals for proteins'):
        if pd.isna(protein):
            continue
        try:
            prot_data = [target.search(protein)[0]]

            # Search for protein with same synonym
            if prot_data == [None]:
                prot_data = target.filter(
                    target_synonym__icontains=protein, target_organism__istartswith=organism
                ).only(['target_chembl_id', 'target_pref_name', 'molecule_chembl_id', 'molecule_pref_name'])
        except HttpBadRequest:
            print(f'Cannot search for {protein} due to chembl error')
            continue

        # No results found
        if not prot_data:
            continue

        for prot in prot_data:
            # Absence of chembl id
            if not prot['target_chembl_id']:
                continue

            prot_activity_data = activity.filter(
                target_chembl_id=prot['target_chembl_id'],
                assay_type_iregex='(B|F)',
            ).only([
                'pchembl_value', 'molecule_chembl_id', 'activity_id', 'target_pref_name', 'molecule_pref_name'
            ])

            if len(prot_activity_data) < 1:
                continue

            for i in prot_activity_data:
                tmp = {}

                try:
                    if i['pchembl_value'] is None:
                        continue
                except HttpApplicationError:
                    continue

                pchembl_val = i['pchembl_value']

                if float(pchembl_val) < 6:
                    tmp['activity'] = 'active'
                else:
                    tmp['activity'] = 'inactive'

                tmp.update({
                    'protein_symbol': protein,
                    'protein_name': i['target_pref_name'],
                    'aid': str(i['activity_id']),
                    'chembl_id': i['molecule_chembl_id'],
                    'compound_name': i['molecule_pref_name'].capitalize() if i['molecule_pref_name'] else ''
                })
                df_data.append(tmp)
            
            #time.sleep(1)
            #print('1 seconds of break')

    # Merge duplicated protein-chemical entries into one
    df = pd.DataFrame()

    for idx, row in tqdm(enumerate(df_data), total=len(df_data), desc='Preparing data'):
        if idx == 0:
            df = df.append(row, ignore_index=True)
        else:
            _in_df = df.loc[
                (df['protein_symbol'] == row['protein_symbol']) & (df['chembl_id'] == row['chembl_id'])
                ]

            if _in_df.empty:
                df = df.append(row, ignore_index=True)
            else:
                row_index = _in_df.index

                # Check existing citations
                existing_assays = set(df.loc[row_index, 'aid'].values[0].split(' | '))
                old_count = len(existing_assays)
                existing_assays.add(row['aid'])
                new_count = len(existing_assays)

                # Check if new citation added, if yes - add respective data
                if old_count < new_count:
                    df.loc[row_index, 'aid'] = ' | '.join(existing_assays)
    df = df[['activity', 'protein_symbol', 'protein_name', 'aid', 'chembl_id', 'compound_name']]
    return df


def target_list_to_chemical(
    proteins: list,
    organism: str = 'Homo sapiens',
) -> pd.DataFrame:
    """Extract chemical information on list of targets
    Usage:
    >> target_list_to_chemical(proteins=['RIPK'])
    """

    df = _get_target_data(protein_list=proteins, organism=organism)
    return df

def getNodeList(nodeName,graph):
    #import pybel
    node_list = []
    for node in graph.nodes():
        if isinstance(node,pybel.dsl.Abundance):
            if node.namespace == nodeName:
                node_list.append(node.name)
    return(node_list)
    
def chembl_annotation(graph):
    chemblids = getNodeList('ChEMBL',graph)
    for item in tqdm(chemblids,desc='Adding ChEMBL URLs'):
        nx.set_node_attributes(graph,{Abundance(namespace='ChEMBL',
                                               name=item):'https://www.ebi.ac.uk/chembl/compound_report_card/'+item},'ChEMBL')
    return graph
    

def chembl_name_annotation(graph, drugs_df):
    # create a dict {'Metformin':'CHEMBL1431'} for adding node attributes
    drugName_chembl_dict = {}
    for v, k in drugs_df[['prefName', 'drugId']].values:
        # print(k,v)
        drugName_chembl_dict.update({k: v})

    molecule = new_client.molecule

    # get all chemicals/drugs from kg
    chemblids = getNodeList('ChEMBL', graph)

    for chem in tqdm(chemblids, desc='Adding preferred names and trade names'):

        #print(chem)

        trade_names = []

        # fetch prefName and synonym from chembl
        getNames = molecule.filter(molecule_chembl_id=chem).only(['pref_name', 'molecule_synonyms'])

        # get preferred name #not preferred because it can be empty sometimes
        # pref_name = getNames[0]['pref_name']

        # get preferred name of a drug from openTarget diseaseAssociatedDrug file #always has a preferredName
        nx.set_node_attributes(graph, {Abundance(namespace='ChEMBL', name=chem): drugName_chembl_dict[chem]},
                               'PreferredName')

        # get trade names and append to a list
        try:
            for item in getNames[0]['molecule_synonyms']:
                if item['syn_type'] == 'TRADE_NAME':
                    trade_names.append(item['molecule_synonym'])

            nx.set_node_attributes(graph, {Abundance(namespace='ChEMBL', name=chem): trade_names}, 'TradeName')

        except:
            continue

    return graph

# def chembl_name_annotation(graph):
    
    # molecule = new_client.molecule
    
    # #get all chemicals/drugs from kg
    # chemblids = getNodeList('ChEMBL',graph)

    # for chem in tqdm(chemblids,desc='Adding preferred names and trade names'):
        
        # trade_names= []
        
        # #fetch prefName and synonym from chembl
        # getNames = molecule.filter(molecule_chembl_id=chem).only(['pref_name','molecule_synonyms'])
        
        # #get preferred name
        # pref_name = getNames[0]['pref_name']
        
        # #get trade names and append to a list
        # for item in getNames[0]['molecule_synonyms']:
            # if item['syn_type'] == 'TRADE_NAME':
                # trade_names.append(item['molecule_synonym'])
                
        # nx.set_node_attributes(graph,{Abundance(namespace='ChEMBL', name=chem):pref_name},'PreferredName')

        # nx.set_node_attributes(graph,{Abundance(namespace='ChEMBL',
                                               # name=chem):trade_names},'TradeName')
    # return graph
    
def getGeneOntolgyNodes(nodeName,graph):
    #import pybel
    node_list = []
    for node in graph.nodes():
        if isinstance(node,pybel.dsl.BiologicalProcess):
            #print(node)
            if node.namespace == nodeName:
                node_list.append(node.name)
    return(node_list)
    
def gene_ontology_annotation(graph,uprotDict):
    
    gobp_dict = {}
    gomf_dict = {}
    
    #create a merged dict of {bioprocess: bp_id} for mapping
    for prot in uprotDict:
        gobp_dict.update(uprotDict[prot]['BioProcess'])
        
    for prot in uprotDict:
        gomf_dict.update(uprotDict[prot]['Function'])
        
    
    #extract bps and mfs from kg
    bp_kg = getGeneOntolgyNodes('GOBP',graph)
    mf_kg = getGeneOntolgyNodes('GOMF',graph)
    
    for item in tqdm(bp_kg,desc='adding biological process annotations'):
        
        gobp_id = gobp_dict[item]
        #print(gobp_id)
        #print('https://www.ebi.ac.uk/QuickGO/term/'+gobp_id)
        nx.set_node_attributes(graph,{BiologicalProcess(namespace='GOBP', name=item):'https://www.ebi.ac.uk/QuickGO/term/'+gobp_id},'QuickGO')
        nx.set_node_attributes(graph,{BiologicalProcess(namespace='GOBP', name=item):gobp_id},'Gene Ontology identifier')
    
    for item in tqdm(mf_kg,desc='adding molecular function annotations'):
        
        gomf_id = gomf_dict[item]
        #print(gobp_id)
        #print('https://www.ebi.ac.uk/QuickGO/term/'+gobp_id)
        nx.set_node_attributes(graph,{BiologicalProcess(namespace='GOMF', name=item):'https://www.ebi.ac.uk/QuickGO/term/'+gomf_id},'QuickGO')
        nx.set_node_attributes(graph,{BiologicalProcess(namespace='GOMF', name=item):gomf_id},'Gene Ontology identifier')
        
    return(graph)
    
#Method 2, using OT dump i.e paraquet files
# function can directly use druggability.to_csv file created above. Much faster

def protein_annotation_druggability(graph):
    
    df = pd.read_csv('../data/misc/DruggableProtein_annotation_OT.csv')
    
    #print(df.head(5))
                                      
    protein_list = []
    for node in graph.nodes():
        if isinstance(node,pybel.dsl.Protein):
            if node.namespace == 'HGNC':
                protein_list.append(node.name)
                
    mapping_dict_id_symbol = {}
    for k,v in df[['approvedSymbol','ENSG']].values:
        mapping_dict_id_symbol.update({k:v})
        
    mapping_dict_prot_drugability = {}
    for k,v in df[['approvedSymbol','Druggable Family']].values:
        mapping_dict_prot_drugability.update({k:v})
        
    for protein in tqdm(protein_list,desc='Adding druggability annotation'):
        
        #print(protein)
        
        if protein in mapping_dict_id_symbol.keys():
                     
            nx.set_node_attributes(graph,{Protein(namespace="HGNC",name=protein):f"https://platform.opentargets.org/target/{mapping_dict_id_symbol[protein]}"},'OpenTargets')
            nx.set_node_attributes(graph,{Protein(namespace="HGNC",name=protein):mapping_dict_prot_drugability[protein]},'Druggability')
          
        else: 
            
            #print(protein)
            nx.set_node_attributes(graph,{Protein(namespace="HGNC",name=protein):'No'},'Druggability')
          
    return(graph)
    

def chembl2rxn_rel(
    chemblid_list,
    graph: BELGraph
) -> BELGraph:
    """

    :param chemblid_list:
    :param graph:
    :return:
    """
    infile = open('data/opentargets/drugReactions.pkl', 'rb')
    rxn_df = pickle.load(infile)
    infile.close()

    chembl_id_rxn = rxn_df[rxn_df['chembl_id'].isin(chemblid_list)]
    chembl_id_rxn = chembl_id_rxn.reset_index(drop=True)
    for i in range(len(chembl_id_rxn)):
        graph.add_qualified_edge(
            Abundance(namespace='ChEMBL', name=chembl_id_rxn['chembl_id'][i]),
            Pathology(namespace='SideEffect', name=chembl_id_rxn['event'][i]),  
            relation='hasSideEffect',
            citation="OpenTargets Platform",
            evidence='DrugReactions'
        )

    return graph

def cas2cid(cas):
    import re
    cid_clean = []
    #cid = [pcp.get_compounds(id,'name') for id in cas]
    for id in tqdm(cas):
        #print(id)
        #retrieves cid for given cas 
        cid = pcp.get_compounds(id,'name')
        #print(cid)
        #strip unwanted strings to get cid 
        filter_cid = re.findall('\((.+?)\)',str(cid))
        #print(filter_cid[0])
        if len(filter_cid) > 0:
            cid_clean.append(filter_cid[0])
        else:
            cid_clean.append('None')

    return(cid_clean)

def cid2chembl(cidList) -> list:
    """Method to convert Pubchem CIDs to ChEMBL ids

    :param cidList:
    :return:
    """

    cid2chembl_list = []

    for id in tqdm(cidList, desc='Converting PubChem ids to ChEMBL ids'):
        c = pcp.Compound.from_cid(id)

        for synonym in c.synonyms:
            if synonym.startswith('CHEMBL'):
                print(synonym)
                cid2chembl_list.append(synonym)

    return cid2chembl_list
  
#function to create sub-graph  
def filter_graph(mainGraph, vprotList):
    nsp_list = []
    chem_list = []
    for u, v, data in mainGraph.edges(data=True):
        if u.name in vprotList or v.name in vprotList:
            #print(u.name)
            nsp_list.append(u)
            #print(u)
            #print(v.name)
            nsp_list.append(v)

    nsp_graph = mainGraph.subgraph(nsp_list)
    #nsp_graph = pybel.struct.mutation.induction_expansion.get_subgraph_by_second_neighbors(mpox_graph, nsp_list, filter_pathologies=False)
    return(nsp_graph)
    
# def chembl2rxn_rel(itmpGraph):
    
    # infile = open('data/normalized_data/drugReactions.pkl','rb')
    # rxn_df = pickle.load(infile)
    # infile.close()
    
    # chembl_id = []
    # for node in itmpGraph.nodes():
        # if isinstance(node,pybel.dsl.Abundance):
            # if node.namespace == 'ChEMBL':
                # chembl_id.append(node.name)
            
    # chembl_id_rxn = rxn_df[rxn_df['chembl_id'].isin(chembl_id)]
    # chembl_id_rxn = chembl_id_rxn.reset_index(drop=True)
    # for i in range(len(chembl_id_rxn)):
        # itmpGraph.add_association(Abundance(namespace='ChEMBL',name = chembl_id_rxn['chembl_id'][i]),
                                  # Pathology(namespace='SideEffect',name = chembl_id_rxn['event'][i]),
                                  # citation = "OpenTargets Platform",evidence = 'DrugReactions')
        
    # return(itmpGraph)
                        
    
# def chembl2rxn_rel(itmpGraph):
    
    # infile = open('data/normalized_data/drugReactions.pkl','rb')
    # rxn_df = pickle.load(infile)
    # infile.close()
    
    # chembl_id = []
    # for node in itmpGraph.nodes():
        # if isinstance(node,pybel.dsl.Abundance):
            # if node.namespace == 'ChEMBL':
                # chembl_id.append(node.name)
            
    # chembl_id_rxn = rxn_df[rxn_df['chembl_id'].isin(chembl_id)]
    # chembl_id_rxn = chembl_id_rxn.reset_index(drop=True)
    # for i in range(len(chembl_id_rxn)):
        # itmpGraph.add_association(Abundance(namespace='ChEMBL',name = chembl_id_rxn['chembl_id'][i]),
                                  # Pathology(namespace='SideEffect',name = chembl_id_rxn['event'][i]),
                                  # citation = "OpenTargets Platform",evidence = 'DrugReactions')
        
    # return(itmpGraph)
                              

#import time
# def _get_target_data(
#     protein_list: list,
#     organism: str
# ) -> pd.DataFrame:
#     """Get chemical for target data from ChEMBL.

#     :param protein_list:
#     :param organism:
#     :return:
#     """
#     df_data = []

#     target = new_client.target
#     activity = new_client.activity

#     for protein in tqdm(protein_list, desc='Retrieving chemicals for proteins'):
#         if pd.isna(protein):
#             continue
#         try:
#             prot_data = [target.search(protein)[0]]

#             # Search for protein with same synonym
#             if prot_data == [None]:
#                 prot_data = target.filter(
#                     target_synonym__icontains=protein, target_organism__istartswith=organism
#                 ).only(['target_chembl_id', 'target_pref_name', 'molecule_chembl_id', 'molecule_pref_name'])
#         except HttpBadRequest:
#             print(f'Cannot search for {protein} due to chembl error')
#             continue

#         # No results found
#         if not prot_data:
#             continue

#         for prot in prot_data:
#             # Absence of chembl id
#             if not prot['target_chembl_id']:
#                 continue

#             prot_activity_data = activity.filter(
#                 target_chembl_id=prot['target_chembl_id'],
#                 assay_type_iregex='(B|F)',
#             ).only([
#                 'pchembl_value', 'molecule_chembl_id', 'activity_id', 'target_pref_name', 'molecule_pref_name'
#             ])

#             if len(prot_activity_data) < 1:
#                 continue

#             for i in prot_activity_data:
#                 tmp = {}

#                 try:
#                     if i['pchembl_value'] is None:
#                         continue
#                 except HttpApplicationError:
#                     continue

#                 pchembl_val = i['pchembl_value']

#                 if float(pchembl_val) < 6:
#                     tmp['activity'] = 'inhibitor'
#                 else:
#                     tmp['activity'] = 'activator'

#                 tmp.update({
#                     'protein_symbol': protein,
#                     'protein_name': i['target_pref_name'],
#                     'aid': str(i['activity_id']),
#                     'chembl_id': i['molecule_chembl_id'],
#                     'compound_name': i['molecule_pref_name'].capitalize() if i['molecule_pref_name'] else ''
#                 })
#                 df_data.append(tmp)
            
#             time.sleep(5)
#             print('5 seconds of break')

#     # Merge duplicated protein-chemical entries into one
#     df = pd.DataFrame()
#     if len(df_data)==0:
#         return df
        

#     for idx, row in tqdm(enumerate(df_data), total=len(df_data), desc='Preparing data'):
#         if idx == 0:
#             df = df.append(row, ignore_index=True)
#         else:
#             _in_df = df.loc[
#                 (df['protein_symbol'] == row['protein_symbol']) & (df['chembl_id'] == row['chembl_id'])
#                 ]

#             if _in_df.empty:
#                 df = df.append(row, ignore_index=True)
#             else:
#                 row_index = _in_df.index

#                 # Check existing citations
#                 existing_assays = set(df.loc[row_index, 'aid'].values[0].split(' | '))
#                 old_count = len(existing_assays)
#                 existing_assays.add(row['aid'])
#                 new_count = len(existing_assays)

#                 # Check if new citation added, if yes - add respective data
#                 if old_count < new_count:
#                     df.loc[row_index, 'aid'] = ' | '.join(existing_assays)
#     print(df)
#     print(df_data)
#     df = df[['activity', 'protein_symbol', 'protein_name', 'aid', 'chembl_id', 'compound_name']]
    
#     return df


# def target_list_to_chemical(
#     proteins: list,
#     organism: str = 'Homo sapiens',
# ) -> pd.DataFrame:
#     """Extract chemical information on list of targets
#     Usage:
#     >> target_list_to_chemical(proteins=['RIPK'])
#     """

#     df = _get_target_data(protein_list=proteins, organism=organism)
#     return df

#def _get_target_data(protein_list: list, organism: str):
#     """Get chemical for target data from ChEMBL"""
#     df_data = []

#     target = new_client.target
#     activity = new_client.activity

#     for protein in protein_list:
#         if pd.isna(protein):
#             continue
#         try:
#             prot_data = [target.search(protein)[0]]

#             # Search for protein with same synonym
#             if prot_data == [None]:
#                 prot_data = target.filter(
#                     target_synonym__icontains=protein, target_organism__istartswith=organism
#                 ).only(['target_chembl_id', 'target_pref_name', 'molecule_chembl_id', 'molecule_pref_name'])
#         except HttpBadRequest:
#             print(f'Cannot search for {protein} due to chembl error')
#             continue

#         # No results found
#         if not prot_data:
#             continue

#         for prot in tqdm(prot_data, f'Analying data for {protein}'):
#             # Absence of chembl id
#             if not prot['target_chembl_id']:
#                 continue

#             prot_activity_data = activity.filter(
#                 target_chembl_id=prot['target_chembl_id'],
#                 assay_type_iregex='(B|F)',
#             ).only([
#                 'pchembl_value', 'molecule_chembl_id', 'activity_id', 'target_pref_name', 'molecule_pref_name'
#             ])

#             if len(prot_activity_data) < 1:
#                 continue

#             for i in prot_activity_data:
#                 tmp = {}

#                 if i['pchembl_value'] is None:
#                     continue

#                 pchembl_val = i['pchembl_value']

#                 if float(pchembl_val) < 6:
#                     tmp['activity'] = 'inhibitor'
#                 else:
#                     tmp['activity'] = 'activator'

#                 tmp['protein_symbol'] = protein
#                 tmp['protein_name'] = i['target_pref_name']
#                 tmp['aid'] = str(i['activity_id'])
#                 tmp['chembl_id'] = i['molecule_chembl_id']
#                 tmp['compound_name'] = i['molecule_pref_name'].capitalize() if i['molecule_pref_name'] else ''
#                 df_data.append(tmp)
                
#             time.sleep(10)

#     # Merge duplicated protein-chemical entries into one
#     df = pd.DataFrame()

#     for idx, row in tqdm(enumerate(df_data), total=len(df_data), desc='Preparing data'):
#         if idx == 0:
#             df = df.append(row, ignore_index=True)
#         else:
#             _in_df = df.loc[
#                 (df['protein_symbol'] == row['protein_symbol']) & (df['chembl_id'] == row['chembl_id'])
#             ]

#             if _in_df.empty:
#                 df = df.append(row, ignore_index=True)
#             else:
#                 row_index = _in_df.index

#                 # Check existing citations
#                 existing_assays = set(df.loc[row_index, 'aid'].values[0].split(' | '))
#                 old_count = len(existing_assays)
#                 existing_assays.add(row['aid'])
#                 new_count = len(existing_assays)

#                 # Check if new citation added, if yes - add respective data
#                 if old_count < new_count:
#                     df.loc[row_index, 'aid'] = ' | '.join(existing_assays)
#     df = df[['activity', 'protein_symbol', 'protein_name', 'aid', 'chembl_id', 'compound_name']]
#     return df

# def target_list_to_chemical(
#     proteins: list,
#     organism: str = 'Homo sapiens',
#     output_dir: str = ''
# ) -> None:
#     """Extract chemical information on list of targets
#     Usage:
#     >> target_list_to_chemical(proteins=['RIPK'])
#     """

#     df = _get_target_data(protein_list=proteins, organism=organism)
#     #os.makedirs(output_dir, exist_ok=True)
#     #df.to_csv(os.path.join(output_dir, 'chemical_annotated.csv'), sep='\t', index=False)
#     return(df)
   
#function to annotate protein nodes with uniprot/3dbionotes
   
# def uniprot_rel_test(
    # named_uprotList,
    # org,
    # graph
# ) -> BELGraph:
    # """Method to add UniProt related edges

    # :param named_uprotList:
    # :param org:
    # :param graph:
    # :return:
    # """
    # for item in named_uprotList:
        # fun = list(named_uprotList[item]['Function'].keys())
        # bp = list(named_uprotList[item]['BioProcess'].keys())
        # for f in fun:
            # if str(named_uprotList[item]['Gene']) != 'nan' and not isinstance(named_uprotList[item]['Gene'], dict):
                # graph.add_association(
                    # Protein(namespace=org, name=named_uprotList[item]['Gene']),
                    # BiologicalProcess(namespace='GOMF', name=f),
                    # citation='UniProt database',
                    # evidence='UniProt query'
                # )
            # else:
                # graph.add_association(
                    # Protein(namespace=org, name=item),
                    # BiologicalProcess(namespace='GOMF', name=f),
                    # citation='UniProt database',
                    # evidence='UniProt query'
                # )

        # for b in bp:
            # if str(named_uprotList[item]['Gene']) != 'nan' and not isinstance(named_uprotList[item]['Gene'], dict):
                # graph.add_association(
                    # Protein(namespace=org, name=named_uprotList[item]['Gene']),
                    # BiologicalProcess(namespace='GOBP', name=b),
                    # citation='UniProt database',
                    # evidence='UniProt query'
                # )
            # else:
                # graph.add_association(
                    # Protein(namespace=org, name=item),
                    # BiologicalProcess(namespace='GOBP', name=b),
                    # citation='UniProt database',
                    # evidence='UniProt query'
                # )
               
        
        # if str(named_uprotList[item]['Gene']) != 'nan' and not isinstance(named_uprotList[item]['Gene'], dict):
            # nx.set_node_attributes(graph,{Protein(namespace=org, name=named_uprotList[item]['Gene']):'https://3dbionotes.cnb.csic.es/?queryId='+item},'3Dbio')
        
        # else:
            # nx.set_node_attributes(graph,{Protein(namespace=org, name=item):'https://3dbionotes.cnb.csic.es/?queryId='+item},'3Dbio')
            
            
        
    # return graph
    
    
# get types of action types for chemicals and proteins    
# actionType_list = []
# for item in chem2mech:
    # for i in range(len(chem2mech[item])):
        # if chem2mech[item][i]['action_type'] != None:
            # if chem2mech[item][i]['action_type'] not in actionType_list:
                # actionType_list.append(chem2mech[item][i]['action_type'])
                
# actionType_list



# def chem2moa_rel(named_mechList, org, graph: BELGraph) -> BELGraph:
    # """Method to create the graph"""
    # pos = ["POSITIVE ALLOSTERIC MODULATOR", "AGONIST", "ACTIVATOR", "PARTIAL AGONIST"]
    # neg = ["INHIBITOR", "NEGATIVE ALLOSTERIC MODULATOR", "ANTAGONIST", "BLOCKER"]
    # misc = [
        # "MODULATOR",
        # "DISRUPTING AGENT",
        # "SUBSTRATE",
        # "OPENER",
        # "SEQUESTERING AGENT",
    # ]

    # for chembl_name, chembl_entries in named_mechList.items():
        # for info in chembl_entries:
            # graph.add_association(
                # Abundance(namespace="ChEMBL", name=chembl_name),
                # BiologicalProcess(namespace="MOA", name=info["mechanism_of_action"]),
                # citation="ChEMBL database",
                # evidence="ChEMBL query",
            # )

            # if not info["target_chembl_id"]:
                # continue

            # if "Protein" in info:
                # if info["action_type"] in pos:
                    # graph.add_increases(
                        # Abundance(namespace="ChEMBL", name=chembl_name),
                        # Protein(namespace=org, name=info["Protein"]),
                        # citation="ChEMBL database",
                        # evidence="ChEMBL query",
                    # )
                # if info["action_type"] in neg:
                    # graph.add_decreases(
                        # Abundance(namespace="ChEMBL", name=chembl_name),
                        # Protein(namespace=org, name=info["Protein"]),
                        # citation="ChEMBL database",
                        # evidence="ChEMBL query",
                    # )

                # if info["action_type"] in misc:
                    # graph.add_association(
                        # Abundance(namespace="ChEMBL", name=chembl_name),
                        # Protein(namespace=org, name=info["Protein"]),
                        # citation="ChEMBL database",
                        # evidence="ChEMBL query",
                    # )

    # return graph
    
# def chem2act_rel(named_ActList, org, graph: BELGraph) -> BELGraph:
    # """Method to add bioassay edges to the KG."""
    # for chemical, chem_entries in named_ActList.items():
        # for chem_data in chem_entries:
            # if chem_data["target_chembl_id"]:
                # if "Protein" in chem_data:
                    # graph.add_association(
                        # Abundance(
                            # namespace="ChEMBLAssay", name=chem_data["assay_chembl_id"]
                        # ),
                        # Protein(namespace=org, name=chem_data["Protein"]),
                        # citation="ChEMBL database",
                        # evidence="ChEMBL query",
                    # )

            # graph.add_association(
                # Abundance(namespace="ChEMBL", name=chemical),
                # Abundance(namespace="ChEMBLAssay", name=chem_data["assay_chembl_id"]),
                # citation="ChEMBL database",
                # evidence="ChEMBL query",
                # annotation={
                    # "assayType": chem_data["assay_type"],
                    # "pChEMBL": chem_data["pchembl_value"],
                # },
            # )

    # return graph