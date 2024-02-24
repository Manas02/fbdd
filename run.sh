#!/bin/bash

# Function to get directory of a file
get_directory() {
    echo "$(cd "$(dirname "$1")" && pwd)"
}

# Create directories
mkdir -p data data/raw data/final data/processed

# Download dataset
wget -c https://media.githubusercontent.com/media/molecularsets/moses/master/data/dataset_v1.csv
mv dataset_v1.csv ./data/raw/moses_dataset.csv

# List of Python file names in the desired order
python_files=("src/fragment_tokenizer/moses_train_test_split.py" 
              "src/fragment_tokenizer/encode_moses.py"
              "src/model/tokenizer.py"
              "src/model/model.py"
)

# Iterate over each Python file and run it
for file in "${python_files[@]}"; do
    echo "Running $file..."
    dir=$(get_directory "$file")
    cd "$dir" || exit
    python "$(basename "$file")"
    cd - || exit
    echo "Finished running $file."
done
