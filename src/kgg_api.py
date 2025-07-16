from utils_v2 import *

def createKG(queryDisease,disease_df,phase):

    ct_phase = phase
    doid=disease_df

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
    
    # if vir_prot:
    #     vir_uprot_ext = ExtractFromUniProt(vir_prot)

    if not drugs_df.empty:
    
        #print('A total of ' + str(len(list(set(drugs_df['drugId'])))) + ' drugs have been identified. Now fetching relevant data')
        
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
    
    # if vir_prot:
    #     kg = uniprot_rel(vir_uprot_ext,'VP',kg)
    
    #snp_dgnet = GetDiseaseSNPs(doid['id'][efo_id])  
    
    # if snp_dgnet != None:
        # kg = snp2gene_rel(snp_dgnet,kg)
    
    if not dis2snp_df.empty:
        kg = snp2gene_rel(dis2snp_df,kg)
    
    #print('Your KG is now generated!','\n')
    
    saveFiles(queryDisease, uprot_df, adv_effect, kg, drugs_df, dis2snp_df,uprot_ext)
    
    return(kg)