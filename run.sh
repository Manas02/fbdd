#!/bin/bash

mkdir data/raw data/final data/processed

wget -c https://media.githubusercontent.com/media/molecularsets/moses/master/data/dataset_v1.csv
mv dataset_v1.csv ./data/raw/moses_dataset.csv

# List of Python file names in the desired order
python_files=("./src/fragment_tokenizer/moses_train_test_split.py", 
              "./src/fragment_tokenizer/encode_moses.py",
              "./src/model/tokenizer.py",
              "./src/model/model.py"
)

# Iterate over each Python file and run it
for file in "${python_files[@]}"; do
    echo "Running $file..."
    python "$file"
    echo "Finished running $file."
done
