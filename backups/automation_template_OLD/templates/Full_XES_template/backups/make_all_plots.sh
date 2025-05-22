#!/bin/bash
output=$(ls output)

# This is just an example, it will need to be edited depending on the template
workdir=$(pwd)
for entry in $output; do
    echo "Making plots for ${entry}"

    workdir="output/${entry}/LFMcalc/QueueStuff/"
    cd "${workdir}/10elec_0.0001_testUD/"
    python plot_temp.py &> plot_temp.out
    cp "${entry}_CT_param_explore.png" "/home/ccardot3/QuantyScripts/automation/output_images/"
    cd ../10elec_0.0001_testUD_extreme
    python plot_temp.py &> plot_temp.out
    cp "${entry}_CT_param_explore_extreme.png" "/home/ccardot3/QuantyScripts/automation/output_images/"

    cd $workdir
    
done
