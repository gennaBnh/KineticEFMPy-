import numpy as np
from math import *
import basicFunctions as function
import cvxEfm 
import submodelCreation
import xlwt
import asp_results 


 #------------------EFM SOLVING----------------
def efm_solving(model,efm,dico_params):
    
    #CREATION OF THE SUBMODEL ASSOCIATED TO THE EFM (CONTAINS ONLY THE METABOLITES AND REACTIONS IN THE EFM)
    submodel = submodelCreation.submodel(model,efm)

    #W CONTAINING ONLY THE REACTIONS OF THE EFM
    W_efm = dico_params["W"].iloc[:,function.find_efm(efm)]

    #DICTIONARY CONTAINING THE INDICES OF CERTAIN METABOLITES IN THE MODEL 
    dico_metabolites_indexes = submodelCreation.metabolites_indexes(submodel)

    #EFM'S RESOLUTION WITH CVXYPY
    result = cvxEfm.cvx_problem(submodel,efm,dico_params,model,W_efm,dico_metabolites_indexes) 

    return result


 #------------------STORAGE OF RESULTS IN AN EXCEL FILE-----------------
def result_file_redaction(model,efm,result_cvx,efm_number):

    #EXTRACTION OF CVX RESULTS
    flux  = []
    for value in efm :
        flux.append(float(result_cvx["alph"])*float(value)*3600)
    
    enzyme = np.transpose(result_cvx["Enz"])
    concentration = np.transpose(result_cvx["X"])
    metabolite = [model.metabolites]
    reactions = [model.react_name]
    

    #WRITING THE RESULTS FILE
    workbook = xlwt.Workbook()
    function.excel_file(workbook,"Metabolites",["Metabolite", "Concentration"],metabolite,concentration,"")
    function.excel_file(workbook,"Enzymes",["Reaction","Flux (mmol.gDW-1.h-1)","Enzyme"],reactions,flux,enzyme)
    workbook.save('results/results_efm_'+str(efm_number)+'.xls')

    return 0


#------------------APPLICATION OF "efm_solving" AND "result_file" ON THE EFM(s) OF INTEREST-----------------
def efm_definition(efm_number,model,dico_param):

    #INDICES OF THE EFM(s) TO SOLVE AMONG THOSE HAVING A BIOMASS DIFFERENT FROM 0
    k_opt = [] 
    efm_name = []

    for i in range(len(dico_param["lobj"])):
        if efm_number == "all":
            efm_name.append(model.efms.axes[0][dico_param["lobj"][i]]) 
            k_opt.append(i) 
        else: #If only a particular efm to solve, storage of its number and indexe 
            if model.efms.axes[0][dico_param["lobj"][i]] == int(efm_number) :
                efm_name.append(model.efms.axes[0][dico_param["lobj"][i]])
                k_opt.append(i)

    #RESOLUTION OF THE ENZYMATIC PROFILE OF THESE EFM(s)
    indice = 0

    for efm_indexe in k_opt: 
        efm = dico_param["EM"].iloc[dico_param["lobj"][efm_indexe],:] #storage of the efm to solve 
        for i in efm :
            i = float(i) 
        print("efm :",efm_name[indice]) #display of the efm solved 
        cvx_result = efm_solving(model,efm,dico_param) #efms solving with CVXPY
        result_file_redaction(model,efm,cvx_result,efm_name[indice]) #redaction of the results file 
        indice += 1

    return 0





