import os

#function is used to extract the affinities from ADV output files
#the ADV output files are single log files for each compound that is docked where the affinities can be found in lines 27-31
def extract():
    
    affinity_values = []
    
    for i in range(0,len(os.listdir('mmp9_vina_run1/output_logs'))):
        try:
            with open(f'mmp9_vina_run1/output_logs/sm_{i}.txt', 'r') as file:
                lines = file.readlines()
                line27= lines[26]
                affinity_values.append(line27.split()[1])
        
        except:
            affinity_values.append(0)      

    with open('mmp9_vina_run1/affinity_list.txt', 'w') as file:
            for item in affinity_values:
                file.write(str(item) + '\n')

    return 'Extraction finished'

extract()