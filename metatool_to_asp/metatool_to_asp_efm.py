all_efms = {}
dict_all_efms = {}

#Stockage des efms du fichier pkl dans un dictionnaire 
i = 0
for _,element in efms.iterrows():
    last_elem = element.loc[element != 0].index
    last_elem = list(last_elem)
    all_efms[efms.axes[0][i]] = last_elem
    last_elem_dict = {x: bool(x in last_elem) for x in col}
    dict_all_efms[efms.axes[0][i]] = last_elem_dict
    i += 1
    
#Stockage des réactions de l'efm de metatool dans une liste  
def read_meta (file) :
    f = open(file,"r")
    for readline in f :
        list_elem = readline.split(" ")
    for elem in list_elem :
        if '0' in elem or '1' in elem or '2' in elem or '3' in elem or '4' in elem or '5' in elem or '6' in elem or '7' in elem or '8' in elem or '9' in elem or 'irreversible' in elem:
            list_elem.remove(elem)
    return list_elem

#Identification l'efm correspondant dans le pkl
def comparaison(l_metatool, dict_pkl):
    l_metatool.sort()
    efm_same_length = {}
    for efm in dict_pkl: #stockage des efms ayant le même nombre de réactions que celui de metatool dans un dico
        if len(dict_pkl[efm]) == len(l_metatool):
            efm_same_length[efm] = dict_pkl[efm]
    for efm in efm_same_length: #recherche de l'efm ayant les mêmes réactions parmi les efms ayant le même nombre de réactions
        if(set(efm_same_length[efm]) == set(l_metatool)):
            return efm