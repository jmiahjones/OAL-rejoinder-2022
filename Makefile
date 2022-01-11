
run_all_sims: R/call_sims.R
	nohup Rscript ./R/call_sims.R > logs/call_sims.log 2>&1

test: tests/test-glmnet.R
	Rscript ./tests/test-glmnet.R

plot: R/visualize.R
	nohup Rscript ./R/visualize.R  > logs/visualize.log 2>&1

plots/positivity_violations.png: plot

plots/ate_bias.png: plot
