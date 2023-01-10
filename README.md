# GammaCygni
Analysis on GammaCygni region

The starting files are saved locally in a parent folder

* dataset_generator
    - You can use this notebook to create all the datasets you need for each "window"
* sliding_window_nullhypothesis
    - You can use this notebook to run the fit on each "window", using one of the datasets created using dataset_generator. Here you can only use the simulated model, but you can remove the "diffuse" component.
* sliding_window_faster
    - Same as sliding_window_nullhypothesis, but here you can remove the "diffuse" component and add a different model (spectral: PWL, constant; spatial: DISK, GAUSSIAN, HESS-LIKE).
* check_results
    - In this notebook you can check the parameters of your fit, the TS and the residuals for each window
* modules
    - useful functions to inspect the results