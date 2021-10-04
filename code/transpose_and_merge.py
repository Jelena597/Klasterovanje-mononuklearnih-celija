import glob

import pandas as pd

GENI_DIR = '../data/raw/std/'
SPISAK = '../data/raw/redukovan_spisak_gena_faza_2.csv'

spisak = ['hg38_' + g for g in pd.read_csv(SPISAK)['Gene name'].values.tolist()]

# print(spisak)


with open('data/interim/geni.csv', 'w') as output_file:
    for i, file in enumerate(glob.glob(GENI_DIR + '*.csv')):
        print(i)
        file_name = file.split('/')[-1][:-4]
        print(file_name)
        df = (
            pd
            .read_csv(file, index_col=0)
            .T
            .filter(spisak, axis=1)
        )
        df.insert(0, 'file', file_name)

        print(len(df.columns))
        # print(df)

        df.to_csv(output_file, header=~bool(i))
        # quit()
