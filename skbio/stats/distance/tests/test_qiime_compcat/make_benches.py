import numpy as np
from skbio import DistanceMatrix
import pandas as pd

sizes = [10, 50, 100, 150, 200, 300, 500, 1000]

for i in sizes:
    mf_file_name = 'rand_mf' + str(i) + '.txt'
    gr = []
    col = ['#SampleID', 'Treatment']
    for i in range(i):
        if i%2==0:
            gr.append([str(i), 'foo'])
        else:
            gr.append([str(i), 'bar'])
        
    df = pd.DataFrame(gr, columns=col)
    df.set_index('#SampleID', inplace=True)
    df.to_csv(mf_file_name, sep='\t')
