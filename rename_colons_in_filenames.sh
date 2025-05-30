#!/bin/bash
# filepath: /home/onur/repos/Allergen-Prediction/rename_files.sh

# Check if a directory argument is provided
if [ $# -eq 0 ]; then
    echo "Usage: $0 <directory>"
    echo "Example: $0 /path/to/folder"
    exit 1
fi

# Get the directory from the first argument
DIR="$1"

# Check if the directory exists
if [ ! -d "$DIR" ]; then
    echo "Error: Directory '$DIR' does not exist"
    exit 1
fi

# Counter for renamed files
count=0

# Find all files with colons in their names
find "$DIR" -maxdepth 1 -type f -name '*:*' | while read -r file; do
    # Get the directory and filename separately
    dir=$(dirname "$file")
    filename=$(basename "$file")
    
    # Replace colons with underscores
    newname="${filename//:/_}"
    
    # Construct the full new path
    newpath="$dir/$newname"
    
    # Check if the target filename already exists
    if [ -e "$newpath" ]; then
        echo "Warning: Cannot rename '$filename' - '$newname' already exists"
    else
        # Rename the file
        mv "$file" "$newpath"
        echo "Renamed: '$filename' -> '$newname'"
        ((count++))
    fi
done

echo "Renaming complete. $count files renamed."