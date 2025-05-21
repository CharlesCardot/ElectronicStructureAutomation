#!/bin/bash

# Define your list of files and folders to add
FILES_TO_ADD=(
    "automation_template/"
    "backups/"
	"utils/"
	"materials_examples"
    "CHANGELOG"
    "TODO"
	"README.md"
)

# Add each file/folder
for item in "${FILES_TO_ADD[@]}"; do
    git add "$item"
done

echo "Custom git add complete!"

