
import numpy as np
import cvxpy as cp
import pandas as pd
import math
import xlwt


def create_vector(elem):
    list_elem = []

    for i in elem :
        list_elem.append(i) #Storage of the values in a list

    #Turn the list into a vector
    elem_vector = np.array(list_elem)
    
    return elem_vector

#-------------------------------------------------------------------

def km_sorting(excelfile,model):

    km_matrix =  np.full((len(model.metabolites),len(model.reactions)),0.5) ##creation of a model matrix filled by default to 0.5
    km = pd.read_excel(excelfile,sheet_name="km",header = None) #reading the km values from the excel file

    firstCol = km[0] #storage of km names 

    dicMetab = {}
    dicReact = {}

    a = 0
    for met in model.metabolites : #storage of model metabolite names and their position in the matrix 
        dicMetab[met.name] = a
        a += 1

    b = 0
    for react in model.react_name : #storage of model reaction names and their position in the matrix  
        dicReact[react] = b
        b += 1
    
    a = 0
    for i in firstCol: #for each km  
        km_name = i.split("_")
        react = km_name[1] #storage of the reaction to which the km is associated
        metab = km_name[2] #and the metabolite
        if metab in dicMetab.keys() and react in dicReact.keys():
            km_matrix[dicMetab[metab],dicReact[react]] = km[1][a] #filling of the matrix with the values of km according to the reaction and the metabolite
        a+=1
    return km_matrix

#-------------------------------------------------------------------

def reactions_homogenity(model,element):

    if element == "W":
        dataframe = np.transpose(model.W) 
    else :
        dataframe = np.transpose(model.obj)

    list=[]

    for i in dataframe.iloc[0,:]: #retrieves the names of the reactions in the parameter
        list.append(i)
    dataframe.columns = list #each column has the name of its reaction

    pres = []
    abst = []
    
    ##CHECK IF THE EFM REACTIONS ARE IN THE PARAMETERS, IF NOT, REMOVE THE REACTION 
    if element == "W":
        for i in model.efms.axes[1]: 
            if i in list:
                pres.append(i)
            else :
                abst.append(i)
                model.efms.pop(i)
    else :
        for i in list:
            if i in model.efms.axes[1]:
                pres.append(i)
            else :
                abst.append(i)

    
    #CHECKS IF THE REACTIONS OF THE PARAMETER ARE IN THE EFM, OTHERWISE REMOVES IT
    for i in dataframe.axes[1]:
        if i not in model.efms.axes[1]:
            dataframe.pop(i)
    
    return dataframe

#-------------------------------------------------------------------

def indice_find(elem,A,A_type,type_elem): #return index of the desired metabolite/reaction (A) in a list/matrix (elem)
    ind = 0
    if type_elem == "matrice":
        for i in elem.index :
            if i == A:
                return ind
            ind += 1
    else :
        for i in elem :
            if A_type == "name":
                if i.name == A :
                    return ind
            else :
                if i.index == A :
                    return ind
            ind += 1

def find_efm(M): #return index of the reactions different from 0 in the efm
    ind = 0
    l = []

    for i in M:
        if i != 0:
            l.append(ind)
        ind += 1

    if len(l) == 1: #If only one index in the list returns it 
        return l[0]
    else:
        return l


def find(X): #return array of index of the desired element 
    L=np.array([])
    for i in range(len(X)):
        if X[i] != 0:
            L = np.append(L,i)
    return L


#-------------------------------------------------------------------

def remplissage_val(elem,elem_Submodel,variable): #function to fill an array with the values of the variable filled by cxv 
    #If variable = x : elem = the metabolites of the model, elem_Submodel = metabolites of x 
    #If variable = E : elem = reactions of the model (model.reactions or efm), elem_Submodel = reactions of E (so different from 0 in the efm)
    array = np.zeros((1,len(elem)))
    ind = 0
    for i in elem_Submodel:
        array[0,i] = variable[ind].value
        ind += 1
    return array

#-------------------------------------------------------------------

def dData(model,diffusionReactionNames,reaction_ind): #search for the diffusion data
    ind = 0
    for name in diffusionReactionNames[0]:
        if name == model.react_name[reaction_ind]:
            return ind
        ind += 1
    return []

#-------------------------------------------------------------------
    # x = decision variable
    # l = list of reactant indices in the submodel
    # lmat = list of reactant indices in model.S
    # km = vector of affinity constants
    # st = stoichiometry matrix
    # j = index of the current reaction


def produit(x,l,lmat,km,st,j) :
    res = 1
    if len(l) != 0 : 
        for i in range(len(l)):
            if km[int(lmat[i]),j] != 0:
                res = res*((x[l[i]]/km[int(lmat[i]),j])**abs(st[int(l[i]), j]))
            else :
                raise TypeError('Should no longer contain any null value.')
    return res


def produitBiom(x,l,lmat,km,j):
    res = 1
    if len(l) !=0 :
        for i in range(len(l)):
            if km[int(lmat[i]),j] != 0 : 
                res = res*((x[l[i]]/km[int(lmat[i]),j]))
            else : 
                raise TypeError('Should no longer contain any null value.')
    return res

#-------------------------------------------------------------------

def excel_file(workbook,sheetName,colnames,names,variable1,variable2):

    #TEXT STYLE DEFINITION
    style_titres_colonnes = xlwt.easyxf('font: bold 1 ; alignment: horizontal center, vertical center')
    style_texte = xlwt.easyxf('alignment: horizontal center, vertical center')

    #CREATION OF THE SHEET
    sheet = workbook.add_sheet(sheetName)

    #FILL THE FIRST ROW WITH THE NAME OF THE COLUMNS
    ind = 0
    for name in colnames:
        sheet.write(0,ind,name,style_titres_colonnes)
        ind += 1


    names=names[0] 
    lg_max_word = 0

    if sheetName == "Metabolites":
        for i in range(len(names)):
            sheet.write(i+1, 0, str(names[i]),style_titres_colonnes)
            sheet.write(i+1, 1, str(variable1[i][0]),style_texte)

            #Get the size of the longest value displayed
            lg_max_A_B = max(len(str(names[i])),len(str(variable1[i][0])))
            if lg_max_A_B > lg_max_word:
                lg_max_word = lg_max_A_B
            

    if sheetName == "Enzymes":
        names = names
        variable2=variable2[:,0]

        for i in range(len(names)):
            sheet.write(i+1, 0, str(names[i]),style_titres_colonnes)
            sheet.write(i+1, 1, str(variable1[i]),style_texte)
            sheet.write(i+1, 2, str(variable2[i]),style_texte)

            #Get the size of the longest value displayed
            lg_max_A_B_C = max(len(str(names[i])),len(str(variable1[i])),len(str(variable2[i])))
            if lg_max_A_B_C > lg_max_word:
                lg_max_word = lg_max_A_B_C


    #DEFINITION OF THE WIDTH OF THE COLUMNS ACCORDING TO THE LONGEST VALUE TO DISPLAY 
    for i in range(len(colnames)):
        sheet.col(i).width = (1 + lg_max_word)*256
