#!/bin/bash
materials=$(ls materials)

# Function to recursively rename files
rename_files() {
    local dir=$1
    local entry=$2

    # Iterate through files and directories in the given directory
    for file in "$dir"/*; do
        if [ -f "$file" ]; then
            # Check if the filename contains "NAME"
            if [[ $file == *NAME* ]]; then
                # Construct the new filename by replacing "NAME" with $entry
                new_name="${file//NAME/$entry}"

                # Rename the file
                mv "$file" "$new_name"
            fi
        elif [ -d "$file" ]; then
            # Recursively call the function for subdirectories
            rename_files "$file" "$entry"
        fi
    done
}

mkdir "output"
mkdir "output_assets"
for entry in $materials; do
    echo $entry
    mkdir "output/$entry"

    workdir="output/${entry}"

    # Copy all wrkdir_template files to current working directory
	cp -r wrkdir_template/* "${workdir}/"
	rename_files "${workdir}" "$entry"

    # Copy all python processing files to working directory
	find "${workdir}" -maxdepth 1 -type d -name '*FEFF*' -exec cp write_feffinp.py {}/ \;

    # Edit qsub file to include material name
    sed -i -e "s/QSUBNAME/${entry}/g" "${workdir}/qsub.script"

    cd "${workdir}"
    qsub qsub.script
    cd ../..
    
done
