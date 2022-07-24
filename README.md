# KineticEFM
## Description
Python library to solve enzymatic and metabolic profiles of elementary flow modes

## Installation
### Cloning the repository
`git clone git@github.com:gennaBnh/KineticEFM.git`

### Prerequisites
- [Python 3](https://www.python.org/downloads/) 
- pip install cvxpy
- pip install numpy
- pip install pkl
- pip install cobra
- pip install pandas
- pip install math 
- pip install xlwt
- pip install pickle 

### Files format
In order to run the library, 3 files must be provided : 

- A file contaning the model parameters (km, kcat ...) :  
    Must be in an excel file (example : param_ecoli.xls)  
    NB : the reactions of every sheet must be sorted alphabetically

- A file containing the model to study :  
    Must be in SBML format (example : ccm.xml)


- A file containing the efms of the model :   
    Must be in pkl format (example : output_ccm.pkl)


An example of these files is available in the "data_example" folder.


### Model constants 
The model constants can be modified in the file modelCreation.py, line 42


## Usage: 

- To study every efm identified with the pareto filter in asp-efm:  
`python3 kineticEFM.py [excel_file] [pkl_file] [xml_file]`  
Example with data_example files:  
`python3 kineticEFM.py data_example/param_ecoli.xls data_example/output_ccm.pkl data_example/ccm.xml`  

- To study only an efm among all the efms identified with asp-efm:  
`python3 kineticEFM.py [excel_file] [pkl_file] [xml_file] [efm number]`  
Example with data_example files:  
`python3 kineticEFM.py data_example/param_ecoli.xls data_example/output_ccm.pkl data_example/ccm.xml 4230`

Results are available an the excel file : results_efm_[efm number].xls
