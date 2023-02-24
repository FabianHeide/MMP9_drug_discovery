"""
This set of functions is used to generate and process drug-like parameter features from small molecules in SMILES format convention.

"""

import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski

#function that calculates lipinski values based on a SMILES line code
def lipinski(smiles, verbose=False):

    try:
        molecule_data= []
        for element in smiles:
            mol=Chem.MolFromSmiles(element) 
            molecule_data.append(mol)

        baseData= np.arange(1,1)
        i=0  
        for mol in molecule_data:        

            paramater_MolWt = Descriptors.MolWt(mol)
            paramater_MolLogP = Descriptors.MolLogP(mol)
            paramater_NumHDonors = Lipinski.NumHDonors(mol)
            paramater_NumHAcceptors = Lipinski.NumHAcceptors(mol)

            row = np.array([paramater_MolWt,
                            paramater_MolLogP,
                            paramater_NumHDonors,
                            paramater_NumHAcceptors])   

            if(i==0):
                baseData=row
            else:
                baseData=np.vstack([baseData, row])
            i=i+1      

        columnNames=["MW","LogP","NumHDonors","NumHAcceptors"]   
        parameters = pd.DataFrame(data=baseData,columns=columnNames)
        
    except:
        print('Check for NaN values in the SMILES column')

    return parameters


#normalize IC50 values to 100000000 or below
def normalize(input):
    norm = []

    for i in input['standard_value']:
        if i > 100000000:
          i = 100000000
        norm.append(i)

    input['standard_value_norm'] = norm
    x = input
        
    return x
     

#function that converts IC50 value to pIC50 for better data handling
def pIC50(input):
    pIC50 = []

    for i in input['standard_value_norm']:
        molar = i*(10**-9) #convert nM to M
        pIC50.append(-np.log10(molar))

    input['pIC50'] = pIC50
    x = input.drop('standard_value_norm', 1)
        
    return x
