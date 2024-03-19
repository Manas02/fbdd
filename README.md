# Masked Language Models are Fragment Based Drug Designers

![](./idea.png)

 
## Download, Process and Train FragmentBERT
```bash
bash run.sh
```

## Steps to reproduce : 
- `fragment_tokenizer/tokenizer.py` Defines Fragment Augmentation, Encoder and Decoder.
- `moses_train_test_split.py` Takes in `moses_dataset.csv` and produces train, test moses data.
- `encode_moses.py` Produces failed and augmented encoded moses smiles files.
- `model/tokenizer.py` Creates FragmentBERT Tokenizer.
- `model.py` Train the model and log to WandB.


## Conda

```bash
conda env create -f environment.yml
conda activate fbdd
```
