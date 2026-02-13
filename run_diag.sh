#!/bin/bash

model_name="2D_Ising"
symm=True
Lmax=2000
iteration=0
n=1
num_collect=120
hermitian=False

mkdir -p tmp

for dcut in 10
do
    for unit_size in 1
    do
        for iT in {0..8}
	    do
            for iL in {0..10}
	        do

                job_name="diag_D${dcut}_US${unit_size}_${model_name}_iter${iteration}_iT${iT}_iL${iL}"

                cp code_main_diag.py                           code_tmp_${job_name}.py
                sed -i "s@#model_name#@${model_name}@g"        code_tmp_${job_name}.py
                sed -i "s@#symm#@${symm}@g"                    code_tmp_${job_name}.py
                sed -i "s@#dcut#@${dcut}@g"                    code_tmp_${job_name}.py
                sed -i "s@#Lmax#@${Lmax}@g"                    code_tmp_${job_name}.py
                sed -i "s@#unit_size#@${unit_size}@g"          code_tmp_${job_name}.py
                sed -i "s@#iteration#@${iteration}@g"          code_tmp_${job_name}.py
                sed -i "s@#iT#@${iT}@g"                        code_tmp_${job_name}.py
                sed -i "s@#iL#@${iL}@g"                        code_tmp_${job_name}.py
                sed -i "s@#n#@${n}@g"                          code_tmp_${job_name}.py
                sed -i "s@#num_collect#@${num_collect}@g"      code_tmp_${job_name}.py
                sed -i "s@#hermitian#@${hermitian}@g"          code_tmp_${job_name}.py

                cp job_temp.sh                                 tmp/job_tmp_${job_name}.sh
                sed -i "s@#job_name#@${job_name}@g"            tmp/job_tmp_${job_name}.sh
                sbatch                                         tmp/job_tmp_${job_name}.sh

            done
        done
    done
done
