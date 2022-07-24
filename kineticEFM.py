import modelCreation
import cobra
import optimEFM
import sys


def efmAnalysis(excel_file,pkl_file,xml_file,efm_to_solve="all"):

    #MODEL CREATION
    model_parameters = modelCreation.parameters_extraction(excel_file, pkl_file,efm_to_solve) #Reading the model parameters
    model_sbml = cobra.io.read_sbml_model(xml_file) #Model formalization with cobrapy
    complete_model = modelCreation.cobra_model(model_sbml, model_parameters,excel_file) #Complete model (cobrapy model + parameters)
    dico_parameters = modelCreation.parameters_definition(complete_model) #Dictionary containing new parameters C,W,kobj,lobj,EM
    
    #EFM(s) RESOLUTION
    optimEFM.efm_definition(efm_to_solve,complete_model,dico_parameters) #Identification of the efm(s) of interest and resolution with CVXPY


#--------------MAIN FUNCTION-------------- 
if __name__ == '__main__':

    if len(sys.argv) == 4 : #Calls efmAnalysis with the default setting "efm_to_solve": all efms identified with pareto will be solved
        efmAnalysis(sys.argv[1],sys.argv[2],sys.argv[3]) 

    elif len(sys.argv) == 5 : #Calls efmAnalysis with the number of the specific efm to solve.  
        efmAnalysis(sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4]) 

    else :
        raise TypeError ("You must provide all the requiered files. See instructions.txt for more informations")

