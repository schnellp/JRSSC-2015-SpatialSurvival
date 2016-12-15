A marginal cure-rate proportional hazards model for spatial survival data

Description of Data Set

mcguire.csv is a text file containing comma-separated data. The first row consists of column headers.
Of interest in the paper are the following:
- id: subject ID
- tooth: tooth number
- smoke: smoking status (1 = smoker)
- age: age (years) at baseline
- hygiene1: hygiene rating ([G]ood, [F]air, [P]oor)
- cr0: crown-to-root ratio at baseline
- probe0: baseline probing depth
- mobile0: baseline mobility
- time: time of last observation (event or censoring time)
- censor: right-censoring indicator (1 = censored)

Programs and Use

- load_data.R loads the raw data from mcguire.csv, produces objects useful for the method presented in the paper, and saves the result as NunnData.RData, which is used directly in subsequent scripts.
- MCMCfx_PS_factor_PCURE.R fits the positive stable model via MCMC.
- MCMCfx_Gauss_PCURE_NUfixed.R fits the competing Gaussian random effects model via MCMC. This is a modification of MCMCfx_PS_factor_PCURE.R.
- MCMCfx_Ind_factor_PCURE.R fits the competing independence model via MCMC. This is a modification of MCMCfx_PS_factor_PCURE.R.
- FitPSmodel.R loads the example data set and calls MCMCfx_PS_factor_PCURE.R to fit the positive stable model.
- summarize.R produces summary statistics and various plots from a model fit.
- callSim.R simulates data sets and fits three competing models via the MCMCfx[...] scripts. Results are output to the "results/" folder. The first set of three fits is saved to results/fit.RObject, and a table of results is updated after each set of three fits is completed and saved as results/partialResult.#.RObject.

For questions, please contact Brian Reich at brian_reich@ncsu.edu.