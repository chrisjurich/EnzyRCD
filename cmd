#!/usr/bin/env bash 
python main.py work_dir/enzyme.pdb  --ligand_1 work_dir/ETI_H.sdf \
                                    --ligand_1_name ETI           \
                                    --ligand_2 work_dir/SAH_H.sdf \
                                    --ligand_2_name SAH           \
                                    --n_struct 5

