import scipy.stats as stats
import math
import numpy as np


# FUNCTION - 1  : PRIOR SI DATA
def get_abPrior(MeanPrior,StdPrior):
    """
        Compute starting simulation value of Alpha, Beta from
        input prior values of MEAN,STD 
    """
    aPrior = MeanPrior * MeanPrior / (StdPrior * StdPrior)
    bPrior = (StdPrior * StdPrior) / MeanPrior
    return aPrior,bPrior

# FUNCTION - 2  : STARTING TIME FRAME
def search_StartingTimeStep(IncidenceSerie,CumulIncThreshold,Verbose=False):
    """
        Compute starting time step for future R computation.
        Total incidence from start of the epidemy is added until 
        input treshold is topped.
    """
    # INIT VARS
    TimeMaxnb = IncidenceSerie.shape[0]
    CumulInc = 0
    t = 0
    # compute first possible calculation (FIRST STEP)
    for t in range(0,TimeMaxnb):
        CumulInc = CumulInc + IncidenceSerie[t]
        if Verbose and (IncidenceSerie[t] > 0) : 
            print("{} : (+{}) {} for {}".format(t,IncidenceSerie[t],CumulInc,round(CumulIncThreshold,1)))
        if (CumulInc >= CumulIncThreshold) : break
        t = t + 1
        
    # VERIF ON TIME PERIODS
    if t == TimeMaxnb: 
        print("ERROR : The epidemic is too small to ever get the desired posterior CV.\
 Estimation aborted. Try a higher value for the aimed posterior CV.")
        return None
    else : 
        if Verbose: print("CumulIncThreshold for proper computations topped off at the TimePeriods {} ".format(t))

    StartEstimDate = t+1 # t+1 as it will be used as excluded upper limit : [t-x, t[
    return StartEstimDate

# FUNCTION - 3 : TIME STEPS DEFINITION
def get_TimeStepSlices(IncidenceSerie,CumulIncThreshold,
    RcomputingRange,RcomputingFreq,Verbose=False):
    """
        Find the best starting point in IncidenceSerie for computing R
        and define TimeSteps Domains for R Computing according to :
            RcomputingRange : size of computing slice
            RcomputingFreq : frequency of calculation

        DEPENDENCIES : search_StartingTimeStep
    """
    # INIT VARS
    endTime=[]
    startTime=[]
    TimeMaxnb = IncidenceSerie.shape[0]
    
    StartEstimDate = search_StartingTimeStep(IncidenceSerie,CumulIncThreshold)
    if not StartEstimDate : return startTime,endTime
    
    StartEstimDate = max(StartEstimDate,
                         RcomputingRange)
    t = StartEstimDate
    TimePeriodNb = 1
    endTime=[]
    startTime=[]
    # FIRST VALS      
    endTime.append(t)
    startTime.append(t - RcomputingRange)
    
    #Counting loop
    while t<= (TimeMaxnb - RcomputingFreq):
        t = t + RcomputingFreq
        TimePeriodNb = TimePeriodNb + 1
        endTime.append(t)
        startTime.append(t - RcomputingRange)
    
    # details
    if Verbose :
        print('Estimation Starting Step : {}'.format(StartEstimDate))
        print('Number of Time periods to study : {}'.format(len(startTime)))
        print('Example of period [x,y[ :')
        for i in range (0,5):
            print('{}-{}'.format(startTime[i],endTime[i]))
    
    return startTime,endTime

# FUNCTION - 4 : DiscreteShifted Gamma SIDistr
def DiscreteShiftedGammaSIDistr(k, mean, sd):
    '''
        Compute a shifted cumulative distribution function of an exponential
        distribution on values [k-2 ; k]

        DEPENDENCIES : scipy.stats LIB
    '''
    a = (mean - 1) * (mean - 1) / (sd * sd) # ALPHA param
    b = sd * sd / (mean - 1) # beta param
    
    if k == 0:
        DiscreteShiftedGammaSIDistr = 0
    if k == 1 : 
        DiscreteShiftedGammaSIDistr = k * stats.gamma.cdf(k, a, scale=b) \
            - a * b * stats.gamma.cdf(k, a + 1, scale=b)
    else :
        DiscreteShiftedGammaSIDistr = k * stats.gamma.cdf(k, a, scale=b) \
            + (k - 2) * stats.gamma.cdf(k - 2, a, scale=b) \
            - 2 * (k - 1) * stats.gamma.cdf(k - 1, a, scale=b) \
            + a * b * (2 * stats.gamma.cdf(k - 1, a + 1, scale=b) \
                       - stats.gamma.cdf(k - 2, a + 1, scale=b) \
                       - stats.gamma.cdf(k, a + 1, scale=b))
        
    return DiscreteShiftedGammaSIDistr

# FUNCTION - 5 : serial interval definition
def FinalSIDistributionWithoutUncertainty(SImean,SIstdev,TimeMaxnb):
    """
        Compute the serial interval distribution from the SIMPLE CASE PARAMETER
        WITHOUT uncertainty  : MEAN and STD only.
        TimeMaxnb is the length of the INCIDENCE SERIE

        return : 
        - 2 recomputed values of MEAN and STD
        - a list of SIditribution values from [0; incidenceserie size] for plotting purpose

        DEPENDENCIES : DiscreteShiftedGammaSIDistr , math lib
    """
    SumPi = 0
    SumPiXi = 0
    SumPiXi2 = 0
    SIDistr = []
    
    for i in range(0,TimeMaxnb):
        current_SIDistr = max(0, DiscreteShiftedGammaSIDistr(i, SImean, SIstdev))
        SumPi = SumPi + current_SIDistr
        SumPiXi = SumPiXi + i * current_SIDistr
        SumPiXi2 = SumPiXi2 + i * i * current_SIDistr
        SIDistr.append(current_SIDistr)
        
    MeanSIFinal = SumPiXi
    sdSIFinal = math.sqrt(SumPiXi2 - SumPiXi * SumPiXi)
    #VERIFICATION if inconsistent data
    if SumPi < 0.99 :
        print("The epidemic is too short compared to the distribution of the SI. Estimation aborted.")
        return False
    else :
        return MeanSIFinal,sdSIFinal,SIDistr

# FUNCTION - 6 & 7 : LAMBDA & CALCULATE POSTERIOR
def LambdaCalc(t, MeanSI, stdSI, Incidence):
    """
        Compute the sum of DiscreteGamma distribution for ALL points from [0, t[
        This function can be quite heavy, so there is a break of loop
        once new computed values have passed the main distribution and does not add
        significant values.

        DEPENDENCIES : DiscreteShiftedGammaSIDistr
    """
    Lambda = 0
    for s in range(0 ,t):
        newvalue = Incidence[t - s] * DiscreteShiftedGammaSIDistr(s, MeanSI, stdSI)
        if (s > (MeanSI+2*stdSI) and newvalue < 0.01):
            # print('break LambdaCalc at {}/{}'.format(s,t))
            break
        Lambda = Lambda + Incidence[t - s] * DiscreteShiftedGammaSIDistr(s, MeanSI, stdSI)
    return Lambda

def CalculatePosterior(Incidence,aPrior,bPrior, MeanSI, stdSI, TimePeriodRange):
    """
        compute posterior Alpha and beta parameter for the current steptime according
        to past steps.

        DEPENDENCIES : LambdaCalc
    """
    # with TimePeriodRange = range(startTime = startTime(TimePeriodNb),endTime(TimePeriodNb)
    SumI = 0
    SumLambda = 0
    
    for t in TimePeriodRange:
        SumI = SumI + Incidence[t]
        SumLambda = SumLambda + LambdaCalc(t, MeanSI, stdSI, Incidence)
    
    aPosterior = aPrior + SumI #aPosterior
    bPosterior = 1 / (1 / bPrior + SumLambda) #bPosterior
    
    return [aPosterior,bPosterior]

# FUNCTION - 8 : COMPUTE MU & SIGMA SAMPLE ON UNCERTAINTY PARAMETERS
def compute_mu_sigma_distribution(SampleSizeSI,SI_mean,SI_stdev,mean_vars,stdev_vars):
    """
        Compute and return 2 list of random sample of MEAN,STD according to :
        - SI_mean : main value of MEAN
        - SI_stdev : main value of STD
        - mean_vars{'std': standard dev around MEAN main value,
                    'MAX': maximum value authorised,
                    'min': minimum value authorised}
        - stdev_vars{'std': standard dev around STD main value,
                    'MAX': maximum value authorised,
                    'min': minimum value authorised}

        DEPENDENCIES : Numpy lib, scipy.stats lib
    """
    mu=[]
    sigma=[]
    for k in range(SampleSizeSI):
        mu_cur = stats.norm.ppf(np.random.random(), loc=SI_mean, scale=mean_vars['std'])     
        while (mu_cur < mean_vars['min'] or mu_cur > mean_vars['MAX']):
            mu_cur = stats.norm.ppf(np.random.random(), loc=SI_mean, scale=mean_vars['std'])
            
        sigma_cur = stats.norm.ppf(np.random.random(), loc=SI_stdev, scale=stdev_vars['std'])
        while (sigma_cur < stdev_vars['min'] or sigma_cur > stdev_vars['MAX']):
            sigma_cur = stats.norm.ppf(np.random.random(), loc=SI_stdev, scale=stdev_vars['std'])
        
        mu.append(mu_cur)
        sigma.append(sigma_cur)
    return mu,sigma