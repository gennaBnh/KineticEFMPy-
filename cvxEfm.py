import numpy as np
import cvxpy as cp
from math import *
import basicFunctions as function

#--------------RESOLUTION OF THE ENZYMATIC PROFILE-------------- 
def cvx_problem(submodel_dico,efm,dico_params,model,W_efm,dico_ind):

    ##DEFINITIONS OF THE PROBLEM VARIABLES
    x = cp.Variable(len(submodel_dico["submatrix_metab"]),pos=True) #CVX variable containing the concentration of all the efm metabolites 
    alph = cp.Variable(pos=True) #CVX variable to resolve the flux rate of each reaction
    E = cp.Variable((len(function.find_efm(efm)),1),pos=True) #CVX variable containing the enzymatic concentration of the flows 


    #TURN C AND e INTO VECTORS TO MULTIPLY THEM IN THE OBJECTIVE 
    c_vector = function.create_vector(dico_params["C"].loc[1,:]) #Turns the column of C containing the numerical values into a vector
    e_vector = function.create_vector(efm) #Turns e into a vector


    #DEFINITIONS OF THE PROBLEM OBJECTIVE
    objective = cp.Maximize(alph*(c_vector@e_vector))


    #DEFINITIONS OF THE PROBLEM CONSTRAINTS
    constraints = []
    
    i = 0 #Index of the enzymes of E
    j = 0 #Index of efm reactions

    while j < len(efm) : #For each reaction of the efm 

        if efm[j] != 0 : #if the reaction value is different from 0

            lsub_npArray = function.find(submodel_dico["submatrix"][:,j]<0) #Index of reaction substrates in the submodel matrix 
            lprod_npArray = function.find(submodel_dico["submatrix"][:,j]>0) #Index of reaction products in the submodel matrix 
            lsubstrat_npArray = function.find(model.matrice[:,j]<0) #Index of reaction substrates in the model matrix 
            lproduit_npArray = function.find(model.matrice[:,j]>0) #Index of reaction products in the model matrix 

            #If it is a diffusion reaction, recuperation of the associated diffusion coefficient
            idxDiff = function.dData(model,model.diffData,j) 

            #Applications of diffusion constraints
            if idxDiff != []:
                if len(lsub_npArray) != 1 or len(lprod_npArray) != 1 :
                    print('  In cvxCall: diffusion reaction with more than one substrate or product.')
                elif efm[j] > 0 :   
                    constraints.append((alph * abs(efm[j]) + model.diffData[1][idxDiff] * x[lprod_npArray[0]]) / (model.diffData[1][idxDiff] * x[lsub_npArray[0]]) <= 1)
                else : 
                    constraints.append((alph * abs(efm[j]) + model.diffData[1][idxDiff] * x[lsub_npArray[0]]) / (model.diffData[1][idxDiff]* x[lprod_npArray[0]]) <= 1)
                i+=1
                j+=1
            

            #If it is not an diffusion reaction 
            else : 

                if model.react_name[j] == 'nadh-dh':
                    sH = 3
                    We = (x[dico_ind["k_Hext"]]/model.constantes['H_int'])**(-sH +model.constantes['a']*model.constantes['F']/(log(10)*model.constantes['R']*model.constantes['T'])) * exp(-model.constantes['b']*model.constantes['F']/(model.constantes['R']*model.constantes['T']))
                    Wie = We
                    T1w = (x[dico_ind["k_NADH"]]/model.km[function.indice_find(model.metabolites,'NADH',"name",""),j]) * (x[dico_ind["k_Q8"]]/model.km[function.indice_find(model.metabolites,'Q8',"name",""),j])* (x[dico_ind["k_H"]]/model.km[function.indice_find(model.metabolites,'H',"name",""),j]) * (x[dico_ind["k_Hext"]]/model.km[function.indice_find(model.metabolites,'H-ext',"name",""),j])
                
                elif model.react_name[j] == 'cytbd':
                    sH = 2
                    We = (x[dico_ind["k_Hext"]]/model.constantes['H_int'])**(-sH +model.constantes['a']*model.constantes['F']/(log(10)*model.constantes['R']*model.constantes['T'])) * exp(-model.constantes['b']*model.constantes['F']/(model.constantes['R']*model.constantes['T']))
                    Wie =We
                    T1w =(x[dico_ind["k_O2"]]/model.km[function.indice_find(model.metabolites,'O2',"name",""),j])**0.5 * (x[dico_ind["k_Q8H2"]]/model.km[function.indice_find(model.metabolites,'Q8H2',"name",""),j]) * (x[dico_ind["k_Hext"]]/model.km[function.indice_find(model.metabolites,'H-ext',"name",""),j])**2
                    
                elif model.react_name[j] == 'atp-synth' : 
                    sH = 4
                    Wi = (model.constantes['H_int']/x[dico_ind["k_Hext"]])**(sH +model.constantes['a']*model.constantes['F']/(log(10)*model.constantes['R']*model.constantes['T'])) * exp(-model.constantes['b']*model.constantes['F']/(model.constantes['R']*model.constantes['T']))
                    Wie = Wi
                    T1w = (x[dico_ind["k_ADP"]]/model.km[function.indice_find(model.metabolites,'ADP',"name",""),j]) * (x[dico_ind["k_H"]]/model.km[function.indice_find(model.metabolites,'H',"name",""),j])**4 
                else : 
                    sH = 0
                

                #If it is the biomass reaction
                if j == function.find(dico_params["C"].iloc[1,:])[0] : 
                    T1=function.produitBiom(x,lsub_npArray,lsubstrat_npArray,model.km,j)
                    #T2=function.produitBiom(x,lprod_npArray,lproduit_npArray,model.km,j)
                    T2 == 0
                    constraints.append((((1+T1+T2)*alph*efm[j]) + (model.kcatm.iloc[j,1]*T2*E[i]))/(model.kcatp.iloc[j,1]*T1*E[i]) <=1) 
                    i+=1
                    j+=1

                #For the other reactions
                else :

                    T1=function.produit(x,lsub_npArray,lsubstrat_npArray,model.km,submodel_dico["submatrix"],j)
                    T2=function.produit(x,lprod_npArray,lproduit_npArray,model.km,submodel_dico["submatrix"],j)

                    if float(efm[j]) > 0 :
                        if sH != 0 :
                            constraints.append((((1+T1+T2)*alph*efm[j]) + (model.kcatm.iloc[j,1]*T2*E[i]))/(model.kcatp.iloc[j,1]*T1w*Wie*E[i])<=1)
                        else :
                            constraints.append((((1+T1+T2)*alph*efm[j]) + (model.kcatm.iloc[j,1]*T2*E[i]))/(model.kcatp.iloc[j,1]*T1*E[i])<=1)
                        i+=1
                        j+=1

                    elif float(efm[j]) < 0:
                        if sH != 0 :
                            constraints.append((((1+T1+T2)*alph*abs(efm[j])) + (model.kcatp.iloc[j,1]*T1w*Wie*E[i]))/(model.kcatm.iloc[j,1]*T2*E[i])<=1)
                        else :
                            constraints.append((((1+T1+T2)*alph*abs(efm[j])) + (model.kcatp.iloc[j,1]*T1*E[i]))/(model.kcatm.iloc[j,1]*T2*E[i])<=1)
                        i+=1
                        j+=1  
        else: 
            j+=1  


    #CONSTRAINTS ON METABOLITES 

    #TURN W AND E INTO VECTORS TO MULTIPLY THEM IN A CONSTRAINT
    w_vector = function.create_vector( W_efm.iloc[1,:])
    E_vector = function.create_vector(E)

    constraints.append(w_vector@E_vector<= model.constantes['Etot'])
    constraints.append((x[dico_ind["k_H"]] == model.constantes['H_int']))
    constraints.append(x[dico_ind["k_Hext"]]*1e-3 <= exp(-log(10)*model.constantes['pHMin']))
    constraints.append(exp(-log(10)*model.constantes['pHMax']) <= x[dico_ind["k_Hext"]]*1e-3 )
    idx=len(submodel_dico["submatrix_metab_sorted"][0])+1

    #CONSTRAINTS ON INTERNAL METABOLITES
    sumX=0
    for k in submodel_dico["submatrix_metabInt_idx"]:
        if submodel_dico["submatrix_metab"][k].name != 'H':
            sumX=sumX+x[k]
    constraints.append(sumX <= 300)


   #CONSTRAINTS ON EXTERNAL METABOLITES
   
   #STORAGE OF THE SUBMODEL EXTERNAL METABOLITE INDICES IN THE MODEL
    listofmetabext = []
    ind = 0 
    for met in model.ext_met :
        for metSubmat in submodel_dico["submatrix_metab_sorted"][1]:
            if met == metSubmat:
                listofmetabext.append(ind)
        ind += 1

    #CREATION OF A LIST l CONTAINING THE VALUES OF EXTERNAL METABOLITES SORTED ACCORDING TO THE POSITION OF THE EXTERNAL METABOLITES IN THE MODEL
    listmetextsorted = []
    for i in model.xext.iloc[:,0]:
        listmetextsorted.append(0)

    metInd = 0
    for i in model.xext.iloc[:,0]:
        ind = 0
        for j in model.ext_met:
            if i == j.name :
                listmetextsorted[ind] = model.xext.iloc[metInd,1]
            ind += 1
        metInd += 1


    ind = 0   
    for k in range(len(submodel_dico["submatrix_metab_sorted"][1])):
        constraints.append(x[k] <= listmetextsorted[listofmetabext[ind]])
        ind += 1


    #CONSERVATION LAW
    if (dico_ind["k_ATP"] != "None") and (dico_ind["k_ADP"] != "None"):
        if dico_ind["k_ATP"] > 0  :
            constraints.append((x[dico_ind["k_ATP"]]+x[dico_ind["k_ADP"]]) <= 10.1600)

    if (dico_ind["k_NAD"] != "None") and (dico_ind["k_NADH"] != "None") :
        if dico_ind["k_NAD"] > 0 :
            constraints.append((x[dico_ind["k_NAD"]]+x[dico_ind["k_NADH"]]) <= 2.6830)

    if (dico_ind["k_NADP"] != "None") and (dico_ind["k_NADPH"] != "None"):
        if dico_ind["k_NADP"] > 0 :
            constraints.append((x[dico_ind["k_NADP"]]+x[dico_ind["k_NADPH"]]) <= 0.1221)

    if (dico_ind["k_Q8"] != "None") and (dico_ind["k_Q8H2"] != "None"):
        if dico_ind["k_Q8"] > 0 :
            constraints.append((x[dico_ind["k_Q8"]]+x[dico_ind["k_Q8H2"]]) <= 1.35)

    if (dico_ind["k_ACOA"] != None ) and (dico_ind["k_COA"] != None) and (dico_ind["k_SUCC_COA"] != None) :
        if (dico_ind["k_ACOA"] > 0) and (dico_ind["k_COA"] > 0) and (dico_ind["k_SUCC_COA"] > 0) :
            constraints.append(x[dico_ind["k_ACOA"]]+x[dico_ind["k_COA"]]+x[dico_ind["k_SUCC_COA"]] <= 2.24)
        
        elif dico_ind["k_ACOA"] > 0 and dico_ind["k_COA"] > 0 :
            constraints.append(x[dico_ind["k_ACOA"]]+x[dico_ind["k_COA"]] <= 2.24)

        elif dico_ind["k_SUCC_COA"] > 0 and dico_ind["k_COA"] > 0 :
            constraints.append(x[dico_ind["k_COA"]]+x[dico_ind["k_SUCC_COA"]] <= 2.24)


    kpts = function.indice_find(model.reactions,'pts',"name","")
    kbiom = function.indice_find(model.reactions,'biomass',"name","")

    #CONSTRAINT ON ALPHA 
    constraints.append((alph*efm[int(kpts)]*3600) <= 8.8)

    #CVX PROBLEM DEFINITION AND RESOLUTION
    problem = cp.Problem(objective, constraints)
    problem.solve(gp=True)

    #PROBLEM SATUTS
    print("=================")
    print(problem.status)
    print("=================")

    opt = problem.value * 3600
    cvxStatus = problem.status
    
    #STORAGE OF THE CVX VARIABLES VALUES
    lk = function.find_efm(efm)
    Enz = function.remplissage_val(efm,lk,E)
    X = function.remplissage_val(model.metabolites,submodel_dico["submatrix_metab_idx"],x)

    dico_result = {"opt":opt,"X":X,"Enz":Enz,"alph":alph.value,"cvxStatus":cvxStatus}
    
    return dico_result


