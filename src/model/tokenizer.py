#!/usr/bin/env python
# coding: utf-8

import os
from tokenizers import Tokenizer
from tokenizers.models import WordLevel
from tokenizers.trainers import WordLevelTrainer
from tokenizers.pre_tokenizers import WhitespaceSplit
from tokenizers.processors import TemplateProcessing

# File paths
files = ["../../data/final/train_encoded_moses.smi", 
         "../../data/final/test_encoded_moses.smi",
         "../../data/final/eval_encoded_moses.smi"]
model_path = '../FragmentBERT'
os.makedirs(model_path, exist_ok=True)

# Tokenizer training
tokenizer = Tokenizer(WordLevel(unk_token="[UNK]"))
trainer = WordLevelTrainer(special_tokens=["[UNK]", "[CLS]", "[SEP]", "[PAD]", "[MASK]"])
tokenizer.pre_tokenizer = WhitespaceSplit()
tokenizer.train(files, trainer)

# Tokenizer post-processing
tokenizer.post_processor = TemplateProcessing(
    single="[CLS] $A [SEP]",
    pair="[CLS] $A [SEP] $B:1 [SEP]:1",
    special_tokens=[
        ("[CLS]", tokenizer.token_to_id("[CLS]")),
        ("[SEP]", tokenizer.token_to_id("[SEP]")),
    ],
)

max_length = 128

# Save the tokenizer
tokenizer.enable_truncation(max_length=max_length)
tokenizer.save(f'{model_path}/tokenizer.json')
