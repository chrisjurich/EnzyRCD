#!/usr/bin/env bash 
python main.py             work_dir/enzyme.pdb      \
                --ligand_1 work_dir/ETI.mol2        \
                --ligand_1_name ETI                 \
                --ligand_2 work_dir/SAH.mol2        \
                --ligand_2_name SAH                 \
                --constraints work_dir/RDock.cst    \
                --n_struct 1 

