# Rcomputing_covid19
try out to estimate time-varying reproduction numbers (R factor inspired from Anne Cori works) with covid19API numbers

## 1-API_covid19_requests notebook
Principal goal is to request fresh data from the API and write them in workdata folder, in CSV format : 
    'Confirmed','Deaths','Recovered','Active' and 'Date'

* request according to country list
* filter responses (drop useless data) and format dates
* compute incidence data (new cases per time unit)
* compute some average column in order to smooth incidence

    7 days smoothing for moving average
    0.1 alpha for exponential moving average

* contain some plot tools to visualise data in raw form

## 2-unit_compute_R notebook

Notebook destined to compute on a unique dataframe
It is mostly use to test/build python files : **utils_R.py, TEST_utils_R.py**

>the MAIN Rcomputation cell is copied into **R_main.py** which is used in loop notebook : *modification/improvement on this cell must be copied in R_main functions*

* load tables writed down in workdata folder by Notebook1
* execute some tests on utils_R functions
* select incidence serie on which compute R
* select variables values and test incorrect variables:

    1. **mean and standard deviation** of epidemic transmission rate and, if uncertainty calculation mode, variations authorized around those values : (std, min, max)

    2. **samplesize** : if uncertainty mode, will compute several couples of mean/std around authorized distribution of main mean/std and compute R on these values as well, then will average all results.
    Careful : increase computing time significantly as this give a matrix of [uncertaintySampleSize (number of mean/std couples) * posteriorSamplesize] more computations for EACH time step !!

    3. **computation steps and length**: frequency of computation in time (steps =1, every day), and on which timeframe (length=7, take values over a week to compute R)
    steps effect just reduce computation time as its coarse results
    length effect smooth out Rcurve as it average on incidence variation during time frame

    4. **CVThreshold**(coefficient of variation) : helps determine Minimum number of incident cases to start R computation (first data of epidemic are inconsistent).
    
    5. **meanprior and stdprior** : initial values to start the first computation (as it is a recursive computation). might decide tu use main mean and std values but priors parameter have gret influence over CVthreshold calculations. Need to be sure before any simplification.


Computation result in a dataframe in this form :
timestep    incidence   aPosterior  bPosterior  MeanR   StdR    RQuantile025    RQuantile05     RQuantile25     Rmedian     RQuantile75     RQuantile95     RQuantile975

## 3-loop_R_on_countrylist
Notebook destined to compute/plot on a multiple dataframe and compare countries using the same structure as describe in notebook2

1. contains a computing loop over several countries dataframes. This loop uses functions from **R_main.py** to compute.
Results are written down into CSV files : workdata/resultframe_ **countryname** _Table.csv.


2. contains a second part dedicated to plotting results (from results CSV, so can be used without the computing loop)

* Rrate = f(time) to visualize one country or compare several
* Rrate computed on different incidence types : confirmed, active, death, recovery.
