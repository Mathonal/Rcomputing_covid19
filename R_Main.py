import scipy.stats as stats
import math
import pandas as pd
import numpy as np
import utils_R

def execmain(): print('meuh')

def Rloop(IncidenceSerie,MeanPrior,StdPrior,SI_mean, SI_stdev,
    CVThreshold,ComputationLength,ComputationStep):
    """
        return dataframe['timestep''aPosterior''bPosterior''MeanR''StdR'
        'RQuantile025''RQuantile05''RQuantile25':RQuantile25,
        'Rmedian''RQuantile75''RQuantile95''RQuantile975']
        for each determined timeframe in incidence
    """
    # =============================
    # R ESTIMATION 
    # =============================
    # local VARS 
    aPosterior=[]
    bPosterior=[]
    MeanR=[]
    StdR=[]
    RQuantile025=[]
    RQuantile05=[]
    RQuantile25=[]
    Rmedian=[]
    RQuantile75=[]
    RQuantile95=[]
    RQuantile975=[]

    # STEP 1
    aPrior,bPrior = utils_R.get_abPrior(MeanPrior,StdPrior)
    CumulIncThreshold = 1 / (CVThreshold * CVThreshold) - aPrior
    # Time steps
    startTimes,endTimes = utils_R.get_TimeStepSlices(IncidenceSerie,CumulIncThreshold,
                                     ComputationLength,ComputationStep)
    # VERIF 
    #  If (endTime(TimePeriodNb - 1) < StartEstimDate) : 
    #    "Warning: you are trying to estimate R too early in the epidemic to get the desired posterior CV. Estimation will be performed anyway."
    #  If any startTime[x] > endTime[x]) : "Time period X has its starting date after its ending date. Estimation aborted."
    #  if any endTime[x] > TimeMaxnb : "Time period X ends after the end of the epidemic.  Estimation aborted.  

    #computation visualization variable
    compviz = len(startTimes)/10
    i = 0
    j = 1

    # LOOP ON TIME PERIODS
    for curStart,curEnd in zip(startTimes,endTimes):
        #Compute A & B posterior
        TimePeriodRange = range(curStart,curEnd)
        #print(TimePeriodRange)
        result = utils_R.CalculatePosterior(IncidenceSerie,aPrior, bPrior, SI_mean, SI_stdev, TimePeriodRange)
        aPosterior.append(result[0])
        bPosterior.append(result[1])

        # compute mean and std using last computed a and b 
        MeanR.append(aPosterior[-1] * bPosterior[-1])
        StdR.append(math.sqrt(aPosterior[-1]) * bPosterior[-1]) 

        RQuantile025.append(stats.gamma.ppf(0.025, aPosterior[-1], scale=bPosterior[-1]))
        RQuantile05.append(stats.gamma.ppf(0.05, aPosterior[-1], scale=bPosterior[-1]))
        RQuantile25.append(stats.gamma.ppf(0.25, aPosterior[-1], scale=bPosterior[-1]))
        Rmedian.append(stats.gamma.ppf(0.5, aPosterior[-1], scale=bPosterior[-1]))
        RQuantile75.append(stats.gamma.ppf(0.75, aPosterior[-1], scale=bPosterior[-1]))
        RQuantile95.append(stats.gamma.ppf(0.95, aPosterior[-1], scale=bPosterior[-1]))
        RQuantile975.append(stats.gamma.ppf(0.975, aPosterior[-1], scale=bPosterior[-1]))

        # Keep track of computation
        if i > compviz :
            compviz = compviz + len(startTimes)/10
            print('{}0% complete'.format(j))
            j=j+1
        i=i+1

    incidencelist = IncidenceSerie.loc[(IncidenceSerie.index >= endTimes[0]-1)]
    # DATAFRAMING ALL RESULTS
    testdict = {'timestep':endTimes,
                'incidence':incidencelist,
                'aPosterior':aPosterior,
                'bPosterior':bPosterior,
                'MeanR':MeanR,
                'StdR':StdR,
                'RQuantile025':RQuantile025,
                'RQuantile05':RQuantile05,
                'RQuantile25':RQuantile25,
                'Rmedian':Rmedian,
                'RQuantile75':RQuantile75,
                'RQuantile95':RQuantile95,
                'RQuantile975':RQuantile975}

    rsltdf = pd.DataFrame(testdict)
    print('computation complete')
    # RETURN
    return rsltdf

def Rloop_I(IncidenceSerie,MeanPrior,StdPrior,
    SI_mean, SI_stdev,mean_vars,stdev_vars,
    CVThreshold,ComputationLength,ComputationStep,
    uncertaintySampleSize,posteriorSamplesize):
    
    # =============================
    # R ESTIMATION 
    # =============================
    # local VARS 
    aPosterior=[]
    bPosterior=[]
    MeanR=[]
    StdR=[]
    RQuantile025=[]
    RQuantile05=[]
    RQuantile25=[]
    Rmedian=[]
    RQuantile75=[]
    RQuantile95=[]
    RQuantile975=[]

    # STEP 1
    aPrior,bPrior = utils_R.get_abPrior(MeanPrior,StdPrior)
    CumulIncThreshold = 1 / (CVThreshold * CVThreshold) - aPrior
    # Time steps
    startTimes,endTimes = utils_R.get_TimeStepSlices(IncidenceSerie,CumulIncThreshold,
                                     ComputationLength,ComputationStep)

    # VERIF 
    #  If (endTime(TimePeriodNb - 1) < StartEstimDate) : 
    #    "Warning: you are trying to estimate R too early in the epidemic to get the desired posterior CV. Estimation will be performed anyway."
    #  If any startTime[x] > endTime[x]) : "Time period X has its starting date after its ending date. Estimation aborted."
    #  if any endTime[x] > TimeMaxnb : "Time period X ends after the end of the epidemic.  Estimation aborted.  

    #computation visualization variable
    compviz = len(startTimes)/10
    i = 0
    j = 1

    # =============================
    # R ESTIMATION WITH UNCERTAINTY
    # =============================
    # COMPUTE UNCERTAINTY SAMPLE
    mu,sigma = utils_R.compute_mu_sigma_distribution(uncertaintySampleSize,SI_mean,SI_stdev,mean_vars,stdev_vars)

    # LOOP ON TIME PERIODS
    for curStart,curEnd in zip(startTimes,endTimes) :
        #Compute A & B posterior
        TimePeriodRange = range(curStart,curEnd)
        #print(TimePeriodRange)   

        SampleR = []
        aPosterior_cur = []
        bPosterior_cur = []
        # LOOP ON SAMPLESIZE
        for k in range(len(mu)):
            result = utils_R.CalculatePosterior(IncidenceSerie,aPrior, bPrior, mu[k], sigma[k], TimePeriodRange)    
            aPosterior_cur.append(result[0])
            bPosterior_cur.append(result[1])

            SampleR_curloop = []
            for l in range(posteriorSamplesize):
                SampleR_curloop.append(stats.gamma.ppf(np.random.random(), 
                                                       aPosterior_cur[-1],
                                                       scale=bPosterior_cur[-1]))
            SampleR.append(SampleR_curloop)
        sampleRprint = np.array(SampleR)

        # compute mean and std using last computed a and b 
        MeanR.append(np.mean(sampleRprint))
        StdR.append(np.std(sampleRprint))
        aPosterior.append(np.mean(aPosterior_cur))
        bPosterior.append(np.mean(bPosterior_cur))

        percentiles = np.percentile(sampleRprint,[2.5,5,25,50,75,95,97.5])

        RQuantile025.append(percentiles[0])
        RQuantile05.append(percentiles[1])
        RQuantile25.append(percentiles[2])
        Rmedian.append(percentiles[3])
        RQuantile75.append(percentiles[4])
        RQuantile95.append(percentiles[5])
        RQuantile975.append(percentiles[6])

        if i > compviz :
            compviz = compviz + len(startTimes)/10
            print('{}0% complete'.format(j))
            j=j+1
        i=i+1

    incidencelist = IncidenceSerie.loc[(IncidenceSerie.index >= endTimes[0]-1)]
    # DATAFRAMING ALL RESULTS
    testdict = {'timestep':endTimes,
                'incidence':incidencelist,
                'aPosterior':aPosterior,
                'bPosterior':bPosterior,
                'MeanR':MeanR,
                'StdR':StdR,
                'RQuantile025':RQuantile025,
                'RQuantile05':RQuantile05,
                'RQuantile25':RQuantile25,
                'Rmedian':Rmedian,
                'RQuantile75':RQuantile75,
                'RQuantile95':RQuantile95,
                'RQuantile975':RQuantile975}

    rsltdf = pd.DataFrame(testdict)
    print('computation complete')
    # RETURN
    return rsltdf

