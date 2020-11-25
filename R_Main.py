import scipy.stats as stats
import math
import numpy as np
import utils_R

def execmain(): print('meuh')

def Rloop():



    # DATAFRAMING ALL RESULTS
    testdict = {'timestep':endTime,
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

    # RETURN
    return rsltdf