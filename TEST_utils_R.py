# FUNCTION - 1  : PRIOR SI DATA
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt

from utils_R import get_abPrior
from utils_R import search_StartingTimeStep
from utils_R import get_TimeStepSlices
from utils_R import DiscreteShiftedGammaSIDistr
from utils_R import FinalSIDistributionWithoutUncertainty

# TEST INCIDENCE DATA FROM ORIGINAL PAPER
TEST_INC_LIST = pd.Series([1, 0, 0, 1, 0, 0, 2, 0, 2, 2, 1, 1,
    1, 0, 0, 0, 4, 1, 2, 4, 13, 23, 35, 26, 12, 17, 19, 17, 28,
    23, 27, 27, 11, 21, 21, 25, 31, 103, 96, 69, 58, 48, 33, 25,
    43, 37, 30, 29, 28, 34, 34, 32, 24, 17, 16, 23, 18, 15, 11,
    19, 16, 9, 17, 14, 6, 4, 8, 9, 11, 11, 11, 6, 13, 3, 8, 2, 4,
    4, 6, 4, 5, 4, 8, 2, 3, 2, 5, 3, 3, 4, 3, 1, 2, 0, 1, 3, 3, 1,
    1, 1, 1, 1, 1, 1, 0, 2])

def TEST_get_abPrior():
    flag = True
    testval = get_abPrior(5,5)
    if (not testval[0] == 1 or not testval[1] == 5) :
        print('TESTERROR on get_abPrior function value test')
        flag = False
    if flag  : print("TEST get_abPrior : OK")
    return flag

def TEST_search_starting_timestep(Verbose=False):
    flag = True
    print('TEST search_starting_timestep : starting')
    verifval = search_StartingTimeStep(TEST_INC_LIST,100,Verbose)
    if not verifval == 24 : 
        print('TESTERROR in search_starting_timestep function value test')
        flag =False
    verifval = search_StartingTimeStep(TEST_INC_LIST,100000,Verbose)
    if verifval : 
        print('TESTERROR in search_starting_timestep function test : provoquemax')
        flag =False
    if flag  : print('TEST search_starting_timestep : OK')
    return flag


def TEST_get_TimeStepSlices(Verbose=False):
    flag = True
    print('TEST get_TimeStepSlices : starting')
    start,end = get_TimeStepSlices(TEST_INC_LIST,10,7,1,Verbose)
    if not len(start) == 95 or not len(end) == 95 : 
        print('TESTERROR in get_TimeStepSlices : ouput size (verify loop) ')
        flag =False
    if not start[0]==5 and not start[-1] == 99 : 
        print('TESTERROR in get_TimeStepSlices : start points ')
        flag =False
    if not end[0]==12 and not end[-1] == 106 : 
        print('TESTERROR in get_TimeStepSlices : end points ')
        flag =False
    start,end = get_TimeStepSlices(TEST_INC_LIST,10,20,1,Verbose)
    if not start[0]==0 and not end[-1] == 20 : 
        print('TESTERROR in get_TimeStepSlices : shift start large calclength')
        flag =False
    start,end = get_TimeStepSlices(TEST_INC_LIST,100000,7,1,Verbose)
    if start or end : 
        print ('TESTERROR in get_TimeStepSlices : outsize sample')
        flag =False
    if flag  : print('TEST get_TimeStepSlices : OK')
    return flag

def TEST_DiscreteShiftedGammaSIDistr(Verbose=False):
    flag = True
    value = DiscreteShiftedGammaSIDistr(0, 10, 2)
    if not value == 0.0 : 
        print('TESTERROR on DiscreteShiftedGammaSIDistr : is not 0 when k=0; verify conditions')
        flag =False
    value = DiscreteShiftedGammaSIDistr(1, 2, 1)
    if not value == 0.3678794411714424 :
        print('TESTERROR on DiscreteShiftedGammaSIDistr : when k=1; \
    verify conditions and calculation')
        flag =False
    value = DiscreteShiftedGammaSIDistr(4, 2, 1)
    if not value == 0.05407678538961869 : 
        print('TESTERROR on DiscreteShiftedGammaSIDistr : when k>2; \
    verify conditions and calculation')
        flag =False

    if Verbose :
        # VISUALISATION
        print('Visualisation of DiscreteShiftedGammaSIDistr function :')
        listtest=[]
        for k in range(0,TEST_INC_LIST.shape[0]):
            listtest.append(DiscreteShiftedGammaSIDistr(k,20,5))
        y1 = np.array(listtest)
        x = np.array(range(0,TEST_INC_LIST.shape[0]))
        plt.plot(x, y1, "y-", color='blue')

    if flag  : print('TEST DiscreteShiftedGammaSIDistr : OK')
    return flag

def TEST_FinalSIDistributionWithoutUncertainty():
    flag = True

    comparSIDistr = [0, 0.0006733229161071663, 0.012173576884085532,
        0.0426690088074797, 0.07784691240841637, 0.10405725229122348,
        0.11651350251941256, 0.11644877791280207, 0.10755556091840318,
        0.09375227486846986, 0.07818005452912091, 0.06295462868020994,
        0.04928027310537493, 0.03768540879359242, 0.028258881607430397,
        0.020839313032990067, 0.015148217455912619]

    MeanSIFinal,sdSIFinal,SIDistr = FinalSIDistributionWithoutUncertainty(8.4,3.8,TEST_INC_LIST.shape[0])
    if not MeanSIFinal == 8.40000000001321 or not sdSIFinal == 3.8218661059491086 :
        print('TESTERROR on FinalSIDistributionWithoutUncertainty : computation of mean and std does not gives good results')
        print(MeanSIFinal,sdSIFinal)
        flag = False

    resulttest = SIDistr[:17]
    if not resulttest == comparSIDistr :
        print('TESTERROR on FinalSIDistributionWithoutUncertainty : computation of SIDistr[] does not gives good results')
        print(resulttest)
        flag = False

    if flag  : print('TEST TEST_FinalSIDistributionWithoutUncertainty : OK')
    return flag
# ==================================================
def launch_utils_R_TESTS(Verbose=False):
    print('STARTING TESTS ON UTILS_R FUNCTION')
    TEST_get_abPrior()
    TEST_search_starting_timestep(Verbose)
    TEST_get_TimeStepSlices(Verbose)
    TEST_DiscreteShiftedGammaSIDistr(Verbose)
    TEST_FinalSIDistributionWithoutUncertainty()