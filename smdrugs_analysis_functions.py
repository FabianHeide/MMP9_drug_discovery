"""
This set of functions is used to generate and process drug-like parameter features from small molecules in SMILES format convention.

"""

import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski, rdMolDescriptors

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

            parameter_MolWt = Descriptors.MolWt(mol)
            parameter_MolLogP = Descriptors.MolLogP(mol)
            parameter_NumHDonors = Lipinski.NumHDonors(mol)
            parameter_NumHAcceptors = Lipinski.NumHAcceptors(mol)
            parameter_TPSA = rdMolDescriptors.CalcTPSA(mol)

            row = np.array([parameter_MolWt,
                            parameter_MolLogP,
                            parameter_NumHDonors,
                            parameter_NumHAcceptors,
                            parameter_TPSA])   

            if(i==0):
                baseData=row
            else:
                baseData=np.vstack([baseData, row])
            i=i+1      

        columnNames=['MW','LogP','NumHDonors','NumHAcceptors', 'TPSA']   
        parameters = pd.DataFrame(data=baseData,columns=columnNames)
        
    except:
        print('Check for NaN values in the smiles column')

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


#function that calculates Mann-Whitney statistics and then puts out a summary table
def mannwhitney(descriptor, verbose=False): 
    from numpy.random import seed
    from numpy.random import randn
    from scipy.stats import mannwhitneyu

    # seed the random number generator
    seed(1)

    # actives and inactives
    selection = [descriptor, 'activity_class']
    df = df1[selection]
    active = df[df['activity_class'] == 'active']
    active = active[descriptor]

    selection = [descriptor, 'activity_class']
    df = df1[selection]
    inactive = df[df['activity_class'] == 'inactive']
    inactive = inactive[descriptor]

    # compare samples
    stat, p = mannwhitneyu(active, inactive)
    #print('Statistics=%.3f, p=%.3f' % (stat, p))
    
    # interpret
    alpha = 0.05
    if p > alpha:
        interpretation = 'Same distribution (fail to reject H0)'
    else:
        interpretation = 'Different distribution (reject H0)'
    
    results = pd.DataFrame({'Parameter':descriptor,
                            'Statistics':stat,
                            'p':p,
                            'alpha':alpha,
                            'Interpretation':interpretation}, index=[0])
    filename = 'mannwhitneyu_' + descriptor + '.csv'
    results.to_csv(f'data/{filename}', index=False)

    return results
