import pandas as pd 


df = pd.read_csv('../../data/raw/moses_dataset.csv')
test_smis = df[df["SPLIT"] == 'test']["SMILES"]
train_smis = df[df["SPLIT"] == 'train']["SMILES"]
eval_smis = train_smis.sample(frac=0.1, random_state=42)
train_smis = train_smis[~train_smis.isin(eval_smis)]

def write_smis_to_newline_delimited_smi_file(smis ,fname):
    with open(fname, 'w') as f:
        for smi in smis:
            f.write(smi+'\n')

write_smis_to_newline_delimited_smi_file(train_smis, '../../data/processed/train_moses.smi')
write_smis_to_newline_delimited_smi_file(eval_smis, '../../data/processed/eval_moses.smi')
write_smis_to_newline_delimited_smi_file(test_smis, '../../data/processed/test_moses.smi')