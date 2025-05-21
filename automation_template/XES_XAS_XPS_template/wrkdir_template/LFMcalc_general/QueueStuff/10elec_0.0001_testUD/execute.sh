#!/bin/bash
Uvals="1.0 2.0 3.0 4.0 5.0 6.0 7.0 8.0 9.0 10.0"
Dvals="1.0 2.0 3.0 4.0 5.0 6.0 7.0 8.0 9.0 10.0"
for U in $Uvals
do
for D in $Dvals
do
        mkdir D_$D+U_$U
        cd D_$D+U_$U
        # Make the input_file
        sed -e"s/UUUU/$U/g" ../template.quanty > testing
        sed -e"s/DDDD/$D/g" testing > input_file.quanty
        rm testing
        cp ../diag_mat.py ./
        echo '#!/bin/bash -l' >> qsub.script
        echo '#PBS -l nodes=1:ppn=4' >> qsub.script
        echo "#PBS -N D_$D+U_$U" >> qsub.script
        echo '#PBS -o quanty.sout' >> qsub.script
        echo '#PBS -e quanty.serr' >> qsub.script
        echo '#PBS -V' >> qsub.script
        echo 'cd $PBS_O_WORKDIR' >> qsub.script
        echo 'quanty input_file.quanty &> testing.out' >> qsub.script

        qsub -q max_24 qsub.script 
        cd ..
done
done


