#!/bin/bash
output=$(ls output)

workdir=$(pwd)
for entry in $output; do
    echo "Making plots for ${entry}"

    cd "output/${entry}/"
	python make_job_plots.py
	mv "${entry}_plots" "${workdir}/output_assets/"
	mv "${entry}_summary.pdf" "${workdir}/output_assets/"
	cd $workdir

done
