import re, os
import numpy as np
import pandas as pd

# Make phyml commands for all AA analyses

dat = pd.read_csv('model_data.csv')
pbs_template = open('example.pbs', 'r').read()

counter_var = 0
for c in dat.columns[1:]:
    model_command = list()
#    model_command.append('PHYML_PATH')
    model_command.append('-i '+re.sub('\n', '', c))
    best_model = dat[c][0]
    if 'LG' in best_model:
        model_command.append('-m LG')
    elif 'Rt' in best_model:
        model_command.append('-m RtREV')

    if '+I' in best_model:
        model_command.append('-v e')

    if '+G' in best_model:
        model_command.append('-a e')

    if '+F' in best_model:
        model_command.append('-f e')
    
    model_command.append('-d aa')

    model_str =  ' '.join(model_command)
    temp_name = c[:5]+str(counter_var)

    pbs_script = re.sub('FILE_NAME', model_str,  re.sub('RUN_NAME', temp_name, pbs_template))

    counter_var = counter_var + 1
    
    open(temp_name+'.pbs', 'w').write(pbs_script)

    print 'script made for '+c+' in '+temp_name+'pbs' 
