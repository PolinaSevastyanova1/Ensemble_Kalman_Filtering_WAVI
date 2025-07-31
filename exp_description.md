ensemble_0:
    observe_sinusoid, 5 iterations, G=1, 

ensemble_1:
    observe_sinusoid, 10 iterations, G=0.05,

ensemble_2:
    observe_sinusoid, 20 iterations, G=0.05,

ensemble_3:
    observe_sinusoid, 20 iterations, G=0.001

ensemble_4:
    exp_PIG_driver_s3, reached 35 full iterations, G=0.001
	actually used the old run_computer_model and priors : (

ensemble_onerun:
	exp_driver, only 1 iteration, all my starting bits 
	#1 is 87653 #2 is 10101 #3 is 25252 #4 is 36363 #5 is 484 #6 is 8976543
	1: some off center (melt) and with outliers (glen, weert)
	2: mostly decent, (bump, glen) skewed left, (weert) skewed right, pct melt good.
	3:good, (glen, melt) good, (bump) outliers but interesting, (pct, weert) skewed but all have good range.
	4: most a bit off-center with large spikes
	5: actually quite gaussian with large spikes, can compare positive, zero, and negative pct. 
	6: (pct) bad, large spikes
Best: I think 5 will give me an ideal experiment due to gaussain shape and large range, to see interaction between parameters. Changing intialisation in main ensemble.yaml to 484.
2 also makes a good study case, as the known melt prefactor is narrow, and others wider. 

ensemble__onerun b1-b5: new initialisations of the truth values. I like 5. results stored in the truth plots, G1_1 being the original forcing. I want to keep the vague gradient similar, but more randomised, so I choose G1_5 forcing.
1 is 235487 #2 is 10101 #3 is 43 #4 is 3934 #5 is 999999 #6 is 1234
I will make my value 999999.

First, I want to see my results without my new parameter, with this new G, and with my own files, for 20 iterations, using exp_driver.jl. Later, I also want to vary the pucnocline forcing random seed. 

ensemble_5:
    exp_PIG_driver_s**, 50 iterations, G=0.001

ensemble_onerunr:
    exp_PIG_driver_s3, 1 iteration
My experiments: have fast bin initial_conditions files fed in, and everything in experiments folder

Need: more gaussian initialization, less iterations, have it work with original G. 

Ensemble takes in whatever file is called exp_PIG_driver.jl. Create a new ensemble file and run that for different things. 

