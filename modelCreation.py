import cobra
import pandas as pd
import numpy as np
import basicFunctions as function
import pickle
import asp_results 

#------------------ CREATION OF THE "MODEL" CLASS ---------------------
# Contains all the informations about the model 

class Model():

    def __init__(self, model, metabolites, reactions, react_name, km, kcatp, kcatm, obj, keq, W, xext, diffData, lobjectif, efms, constantes, int_met, ext_met, matrice):
        self.model = model
        self.metabolites = metabolites
        self.reactions = reactions
        self.react_name = react_name
        self.km = km
        self.kcatp = kcatp
        self.kcatm = kcatm
        self.obj = obj
        self.keq = keq
        self.W = W
        self.xext = xext
        self.diffData = diffData
        self.lobjectif = lobjectif
        self.efms = efms
        self.constantes = constantes
        self.int_met = int_met
        self.ext_met = ext_met
        self.matrice = matrice

    #Display functions used for objects of the Mode classl
    def __str__(self):
        return f"Model: {self.model} | Métabolites: {self.metabolites} | Réactions: {self.reactions} | Réactions name: {self.react_name} | KM: {self.km} | KCATP(s-1): {self.kcatp} | KCATM: {self.kcatm} | Obj: {self.obj} | Keq: {self.keq} | W: {self.W} | XEXT: {self.xext} | diffData(s-1): {self.diffData} | lObjectif: {self.lobjectif} | EFMS: {self.efms} | Constantes: {self.constantes} | Met_int: {self.int_met} | Met_ext {self.ext_met} \n"

    def __repr__(self):
        return str(self)


#----------- DECLARATION OF THE MODEL CONSTANTS -----------
constant_dico = {
"z_H" : 1, 
"a" : 33.33e-3,
"b" : -143.33e-3, 
"H_int" : 1e-4, 
"R" : 8.314e-3, 
"T" : 310, 
"F" : 96485e-3, 
"Etot" : 5.64, 
"pHMin" : 5, 
"pHMax" : 9 }


#------------------------ READING OF PARAMETERS ------------------------
def parameters_extraction(xlsFile,pklFile,efm_number):

    #READING OF PARAMETERS FORM THE XLS FILE (objective,w,ext,diffusion,keq,kcatp,kcatm)
    list_param = ["objective","w","ext","diffusion","keq"]
    dico_param = {}

    for param in list_param:
        dico_param[param]=pd.read_excel(xlsFile, sheet_name=param, header = None)
    dico_param["kcatp"] = pd.read_excel(xlsFile,sheet_name="kcat", header = None, index_col=None, usecols = "A:B")
    dico_param["kcatm"] = pd.read_excel(xlsFile,sheet_name="kcat", header = None, index_col=None, usecols = "A,C")
    
    #READING OF EFMS FROM THE PKL FILE
    if efm_number == "all" :
        dico_param["efms"] = asp_results.recup_efms(pklFile)

    else :
        with open(pklFile, 'rb') as f:
            efms = pickle.load(f)
        dico_param["efms"] = efms

    return dico_param


#------------------ CRÉATION DU MODEL ---------------------
def cobra_model(model_sbml,params,xlsFile):
    
    #IDENTIFICATION OF INTERNAL AND EXTERNAL METABOLITES
    int_met = []
    ext_met = []

    for met in model_sbml.metabolites:
        if "ext" in met.name or 'BIOMASS' in met.name:
            ext_met.append(met)
        else:
            int_met.append(met)
    
    #REMOVAL OF REACTIONS ADDED TO THE MODEL BY COBRA 
    nb_react_added_by_cobra = 0
    for react in model_sbml.reactions:
        if "EX_" in react.name :
           nb_react_added_by_cobra += 1

    model_sbml.reactions = model_sbml.reactions[nb_react_added_by_cobra:]

    #NAME OF THE MODEL'S REACTIONS
    names=[]
    for react in model_sbml.reactions:
        names.append(react.name)
    
    #CREATION OF THE MODEL'S STOICHIOMETRIC MATRIX
    stoichiometric_matrix = cobra.util.create_stoichiometric_matrix(model_sbml,'dense',None)
    
    #CREATION OF THE OBJECT "complete_model" THAT CONTAINS ALL INFORMATION 
    complete_model = Model(model_sbml, model_sbml.metabolites, model_sbml.reactions, names, None, params["kcatp"], params["kcatm"],params["objective"],params["keq"] ,params["w"], params["ext"], params["diffusion"], None, params["efms"],constant_dico, int_met, ext_met,stoichiometric_matrix)
    complete_model.km = function.km_sorting(xlsFile,complete_model) 

    return complete_model


#------------------ DEFINITION OF OTHER PARAMETERS -----------------
def parameters_definition(model):

    #KEEPING THE MODEL REACTIONS ONLY IN THE PARAMETERS 
    W_dataframe= function.reactions_homogenity(model,"W")
    C_dataframe = function.reactions_homogenity(model,"obj")

    #INDEX OF THE REACTION TO MAXIMIZE
    kobj_val = function.find(C_dataframe.iloc[1,:])
    col_reaction_to_maximize = model.efms.iloc[:,kobj_val] #Column corresponding to this reaction for efms 

    #INDICES OF THE EFMS FOR WHICH THE REACTION TO BE MAXIMIZED HAS A VALUE DIFFERENT FROM 0
    lobj = []
    for i in range(model.efms.iloc[:,kobj_val].shape[0]):
        if float(col_reaction_to_maximize.iloc[i,0]) != 0:
            lobj.append(i) 

    EM = model.efms 

    dico = {"C":C_dataframe,"W":W_dataframe,"kobj":kobj_val,"lobj":lobj,"EM":EM} #storage of variables in a dictionary
    return dico

#------------------ CREATION OF A DICTIONARY CONTAINING THE INDICES OF CERTAIN METABOLITES/REACTIONS IN THE MODEL -----------------
def met_indices(model):

    k_H = function.indice_find(model.metabolites,"H","name","")
    k_Hext = function.indice_find(model.metabolites,'H-ext',"name","")

    katpase= function.indice_find(model.reactions,'atp-synth',"name","")
    kpfl= function.indice_find(model.reactions,'pfl',"name","")
    knadhDH= function.indice_find(model.reactions,'nadh-dh',"name","")
    kcytbd= function.indice_find(model.reactions,'cytbd',"name","")
    ksdh= function.indice_find(model.reactions,'sdh',"name","")
    kadh= function.indice_find(model.reactions,'adh',"name","")
    
    dico = {"k_H":k_H, "k_Hext":k_Hext, "katpase":katpase, "kpfl":kpfl, "knadhDH":knadhDH, "kcytbd":kcytbd, "ksdh":ksdh, "kadh":kadh}
    return dico
