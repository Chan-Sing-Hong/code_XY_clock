#!/bin/bash

model_name="2D_Ising"
symm=True
Lmax=$((2000))

mkdir -p tmp

for dcut in 10
do
    for unit_size in 1
    do
        for T in $(seq 0.8 0.02 0.96)
	    do

            job_name="HOTRG_D${dcut}_US${unit_size}_${model_name}_T${T}"

            cp code_main_HOTRG.py                          code_tmp_${job_name}.py
            sed -i "s@#model_name#@${model_name}@g"        code_tmp_${job_name}.py
            sed -i "s@#symm#@${symm}@g"                    code_tmp_${job_name}.py
            sed -i "s@#T#@${T}@g"                          code_tmp_${job_name}.py
            sed -i "s@#dcut#@${dcut}@g"                    code_tmp_${job_name}.py
            sed -i "s@#Lmax#@${Lmax}@g"                    code_tmp_${job_name}.py
            sed -i "s@#unit_size#@${unit_size}@g"          code_tmp_${job_name}.py
            
            cp job_temp.sh                                 tmp/job_tmp_${job_name}.sh
            sed -i "s@#job_name#@${job_name}@g"            tmp/job_tmp_${job_name}.sh
            sbatch                                         tmp/job_tmp_${job_name}.sh

        done
    done
done
