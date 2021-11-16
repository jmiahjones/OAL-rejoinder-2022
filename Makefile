
run_all_sims: R/call_sims.R
	nohup Rscript ./R/call_sims.R > logs/call_sims.log 2>&1

results/full_results.qs: R/call_sims.R
	run_all_sims

R/call_sims.R: R/simulate_oal.R

R/simulate_oal.R: R/oal_funs.R

test: tests/test-glmnet.R
	Rscript ./tests/test-glmnet.R

