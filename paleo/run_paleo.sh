PISM_PREFIX=$HOME/local-rl8/pism python run_glacial_cycle.py -w 72:00:00 -e ../uncertainty_quantification/ensemble_gris_paleo_climate-calving_w_posterior_lhs_20.csv -g 12000 -s chinook-rl8 -q t2small --exstep 250 --calving hybrid_calving -n 20 --stress_balance ssa+sia --age --bed_def lc --o_dir 2023_10_climate_calving_w_posterior --o_format netcdf4_parallel --compression_level 2 -i ../../best_v1/g4500m_const_ctrl_e_1.25_ppq_0.6_tefo_0.02_ssa_n_3.25_philow_5.0_hydro_null_100a.nc




