from tqdm import tqdm 

from tokenizer import encoder, decoder, fragment_augmentation
from moses_train_test_split import write_smis_to_newline_delimited_smi_file


def read_smi_file(fname):
    with open(fname, 'r') as f:
        smis = [i.strip() for i in f.readlines()]
    return smis


for mode in ["train", "test", "eval"]:
    smis = read_smi_file(f'../../data/processed/{mode}_moses.smi')

    failed = set()
    successful = set()
    for smi in tqdm(smis):
        for frags in fragment_augmentation(smi, 5):
            try:
                if decoder(encoder(frags)) == smi:
                    successful.update([encoder(frags)])
                else:
                    failed.update([smi])
            except:
                failed.update([smi])

    write_smis_to_newline_delimited_smi_file(failed, f'../../data/final/{mode}_failed_to_encode_moses.smi')
    write_smis_to_newline_delimited_smi_file(successful, f'../../data/final/{mode}_encoded_moses.smi')
