# Simulations for the Outcome-Adaptive Lasso Rejoinder Paper

This is the code for simulations in our rejoinder paper to a reader reaction. The main files are:

- `R/call_sims.R`: This does the "expand grid" logic to create several different simulation scenarios.
Once the scenarios are created, it calls the `simulate_oal()` function.

- `R/simulate_oal.R`: This is where the eponymous function is defined. It simulates the data and calls the OAL function.

- `R/oal_funs.R`: This a big workhorse file, where the (nG)OAL logic and the grid search is described. 
This deserves its own section below.

- `R/ate_funs.R`: Another big workhorse file, where all the ATE estimation logic is contained.

- `R/visualize.R`: This creates the ggplot2 plots for our viewing pleasure.

- `tests/test-oal-utilities.R`: The ATE utility functions are tested to verify they are working as intended.

- `tests/test-glmnet.R`: This is where a lot of the logic around scaling the lambda1 values is demonstrated. See below.

## ATE Estimation Procedures

The `R/ate_funs.R` file contains a lot of logic around estimating the ATE.
First, it defines a big function `ATE_est()` that estimates the ATE using
whatever method is passed in. It either does IPW, AIPW, or TMLE-based estimation.
The TMLE code should roughly match what is used in the package.

Other utility functions involve the addition of truncation in the creation of IPW weights, and the creation of the overlap weights. All of the ATE estimation methods-IPW, AIPW, and TMLE-involve a generic weight function which would usually be used for IPW weights. By making this more general, it is simple to compare the use of truncation or overlap weights across the different estimators.

Finally, a `cross_fitting()` function is also defined that uses cross-fitting with `lm()` to get the outcome regression. The same estimates are used for all of the ATE
estimators. This is based on the `mlr3` package for machine learning, so more complex estimators can be specified by passing in a `mlr3`-compatible name in the `learner` argument.

## OAL Fitting using Glmnet

The function `oal_fitting()` in the `R/oal_funs.R` file performs the fitting of 
the OAL penalization function. When lambda2=0 (default), the original OAL is used. 
If lambda1 is NULL (default), the lambda1 values used are those determined internally by glmnet. 
Otherwise if a lambda1 values is passed in, as the `simulate_oal()` function does, that values is passed into glmnet.
The supplied value for lambda1 will result in a coefficient which minimizes the objective function in beta:

-ln(beta) + lambda1 * sum( pen * abs(beta) ),

where ln() is the log-likelihood, beta are the coefficients, and the function is written in R notation.

To do this, we have to scale lambda1 by `mean(pen)`, as is demonstrated in tests/test-glmnet.R. We also must
account for the difference in row-sizes of the naive GOAL method, if lambda2 != 0.

The `grid_search_oal_fit()` function repeatedly calls `oal_fitting()`, but uses a grid search based on the wAMD criterion.

Other utility functions are also provided in this file, although they are all toward the top and rather straightforward
adaptations of the original code.

## Running the Code

First thing's first, you will need to install the necessary packages. I'm using renv to keep track of my package list.

The first thing you will need to do, if you do not have renv installed, is run `install.packages('renv')` to get this functionality. Then, either open this folder in Rstudio, or at a terminal located in this folder, run R and type the following commands:

```{R}
renv::restore()
```

This will instruct renv to install all the packages at the versions recorded in the `renv.lock` file. It may take several minutes to finish this step.

If you are on a UNIX-like system with GNU Make, you can simply use the commands in the Makefile to run the code.
To run all the simulations, you may run `make` or `make`.
To create the plots afterwards, you may run `make plot`.
Otherwise, you may simply run the `./R/call_sims.R` and `./R/visualize.R` scripts, respectively.

You may wish to make sure you have a `results/` directory beforehand, for the results to be saved in, if it did not get created when you downloaded this code.
