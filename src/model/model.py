#!/usr/bin/env python
# coding: utf-8

import wandb
from transformers import (
    PreTrainedTokenizerFast,
    BertConfig,
    BertForMaskedLM,
    LineByLineTextDataset,
    DataCollatorForLanguageModeling,
    TrainingArguments,
    Trainer,
)


wandb.init(project="FragmentBERT")
model_path = '../FragmentBERT'

tokenizer = PreTrainedTokenizerFast(tokenizer_file=f"{model_path}/tokenizer.json")

# Set special tokens
tokenizer.mask_token = "[MASK]"
tokenizer.unk_token = "[UNK]"
tokenizer.pad_token = "[PAD]"
tokenizer.sep_token = "[SEP]"
tokenizer.cls_token = "[CLS]"

max_length = 128

# Load datasets
train_dataset = LineByLineTextDataset(
    tokenizer=tokenizer,
    file_path="../../data/final/sampled/sampled_2M_train.smi",
    block_size=128,
)

eval_dataset = LineByLineTextDataset(
    tokenizer=tokenizer,
    file_path="../../data/final/sampled/sampled_100K_eval.smi",
    block_size=128,
)

# Data collator for language modeling
data_collator = DataCollatorForLanguageModeling(
    tokenizer=tokenizer, mlm=True, mlm_probability=0.15
)

# Model configuration
model_config = BertConfig(vocab_size=30000, max_position_embeddings=max_length)
model = BertForMaskedLM(config=model_config)

# Model training arguments
training_args = TrainingArguments(
    output_dir=model_path,
    evaluation_strategy="steps",
    overwrite_output_dir=True,
    num_train_epochs=5,
    per_device_train_batch_size =128,
    gradient_accumulation_steps=8,
    per_device_eval_batch_size = 128,
    logging_steps=1,
    save_steps=1,
    load_best_model_at_end=True,
    resume_from_checkpoint=model_path,
    disable_tqdm=False,
    save_total_limit=1,
    report_to='wandb',
)

# Trainer initialization
trainer = Trainer(
    model=model,
    args=training_args,
    data_collator=data_collator,
    train_dataset=train_dataset,
    eval_dataset=eval_dataset
)

# Model training with mixed-precision
trainer.train()
