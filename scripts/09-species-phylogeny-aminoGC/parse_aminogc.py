import sys
import os
import pandas as pd

data = {'orthogroup':[], 'aminogc_dev':[]}
for path in sys.argv[1:-3]:
    orthogroup = os.path.split(path)[1]
    orthogroup = os.path.splitext(orthogroup)[0]
    print(orthogroup)
    with open(path) as f:
        for line in f:
            line = line.strip().split()
            if line[0] == 'STD.DEV.':
                data['orthogroup'].append(orthogroup)
                data['aminogc_dev'].append(line[1])


df = pd.DataFrame.from_dict(data)
df.sort_values(by='aminogc_dev', inplace=True)
df.sort_values(by='aminogc_dev', inplace=True)
df['orthogroup'].head(25).to_csv(sys.argv[-3], header=False, index=False)
df['orthogroup'].head(25).to_csv(sys.argv[-2], header=False, index=False)
df.to_csv(sys.argv[-1], index=False, sep='\t')
