import numpy as np
import basicFunctions as function

#-------------SUBMODEL CREATION--------------
def submodel(model,efm):
    subS = model.matrice[:,function.find_efm(efm)] #stochiometric submatrix containing the reactions contained in the efm that are non-zero only

    metab_idx = [] #list of the efm metabolite indices 
    metab_in_submat = [] #list of the efm metabolite 

    ind_metab_ext = [] #list of the efm external metabolite indices 
    ind_metab_int = [] #list of the efm internal metabolite indices 
    list_int_met = [] #list of the efm internal metabolite 
    list_ext_met = [] #list of the efm external metabolite 
    
    ind_in_submat = 0

    for idx in range(subS.shape[0]): 

        if np.any(subS[idx,:]) : #if the metabolite is contained in at least one reaction of the efm
            metab_idx.append(idx) #storage of the metabolite index 
            metab_in_submat.append(model.metabolites[idx]) #storage of the metabolite

            if "ext" in model.metabolites[idx].name or "BIOMASS" in model.metabolites[idx].name:
                ind_metab_ext.append(ind_in_submat)
                list_ext_met.append(model.metabolites[idx])

            else :
                ind_metab_int.append(ind_in_submat)
                list_int_met.append(model.metabolites[idx])

            ind_in_submat += 1 #indice of the metabolite in the submatrix
    
    submat_metab_e = model.matrice[metab_idx,:] #creation of the submatrix containing the metabolites contained in the efm only
    
    #Storage of all information in a dictionary 
    dico_e_infos = {"submatrix":submat_metab_e, "submatrix_metab_idx":metab_idx, "submatrix_metab":metab_in_submat,"submatrix_metabExt_idx":ind_metab_ext,"submatrix_metabInt_idx":ind_metab_int,"submatrix_metab_sorted":[list_int_met,list_ext_met]}

    return dico_e_infos

#-------------INDICES OF CERTAIN METABOLITES IN THE SUBMODEL--------------
def metabolites_indexes(dico):

    k_H=function.indice_find(dico["submatrix_metab"],'H',"name","")
    k_Hext=function.indice_find(dico["submatrix_metab"],'H-ext',"name","")
    k_ATP=function.indice_find(dico["submatrix_metab"],'ATP',"name","")
    k_ADP=function.indice_find(dico["submatrix_metab"],'ADP',"name","")
    k_NAD=function.indice_find(dico["submatrix_metab"],'NAD',"name","")
    k_NADH=function.indice_find(dico["submatrix_metab"],'NADH',"name","")
    k_NADP=function.indice_find(dico["submatrix_metab"],'NADP',"name","")
    k_NADPH=function.indice_find(dico["submatrix_metab"],'NADPH',"name","")
    k_Q8=function.indice_find(dico["submatrix_metab"],'Q8',"name","")
    k_Q8H2=function.indice_find(dico["submatrix_metab"],'Q8H2',"name","")
    k_ACOA=function.indice_find(dico["submatrix_metab"],'ACETYL-CoA',"name","")
    k_COA=function.indice_find(dico["submatrix_metab"],'CoASH',"name","")
    k_SUCC_COA=function.indice_find(dico["submatrix_metab"],'SUCC-CoA',"name","")
    k_O2=function.indice_find(dico["submatrix_metab"],'O2',"name","")

    dico = {"k_H":k_H, "k_Hext":k_Hext, "k_ATP":k_ATP, "k_ADP":k_ADP, "k_NAD":k_NAD, "k_NADH":k_NADH, "k_NADP":k_NADP, "k_NADPH":k_NADPH,"k_Q8":k_Q8, "k_Q8H2":k_Q8H2, "k_ACOA":k_ACOA, "k_COA":k_COA, "k_SUCC_COA":k_SUCC_COA, "k_O2":k_O2}
   
    return dico
