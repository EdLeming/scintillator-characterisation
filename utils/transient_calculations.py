import ROOT
import numpy as np

colors = [ROOT.kBlack, ROOT.kRed-4, 
          ROOT.kGreen+2, ROOT.kCyan+2,
          ROOT.kBlue, ROOT.kMagenta+1,
          ROOT.kGray, ROOT.kOrange+7]

def rms(alist):
    '''Calc rms of 1d array'''
    if len(alist) > 1:
        listsum = sum((i - np.mean(alist))**2 for i in alist)
        return np.sqrt(listsum/(len(alist) - 1.0))
    else:
        #print "More than one item needed to calculate RMS, returning 0"
        return 0.
    
def positiveCheck(y):
    if np.abs(max(y[0,:])) > np.abs(min(y[0,:])):
        return True
    else:
        return False

def interpolateThreshold(x, y, thresh, rise=True, start=0):
    """Calculate the threshold crossing using a linear interpolation"""
    if rise == True:
        index_high = np.where( y > thresh )[0][start]
    else:
        index_high = np.where( y < thresh )[0][start]
    index_low = index_high - 1
    dydx = (y[index_high] - y[index_low])/(x[index_high]-x[index_low])
    time = x[index_low] + (thresh - y[index_low]) / dydx
    return time

def calcArea(x,y):
    """Calc area of pulses"""
    integral = np.zeros( len(y[:,0]) )
    for i in range(len(y[:,0])):
        integral[i] = np.trapz(y[i,:],x)
    return np.mean(integral), rms(integral), integral

def calcPulseRatios(x, y1, y2):
    """Calc the ratio of the pulse areas"""
    ratio = np.zeros( len(y1[:,0]) )
    for i in range( len(y1[:,0]) ):
        area1 = np.trapz(y1[i,:], x)
        area2 = np.trapz(y2[i,:], x)
        ratio[i] = area1 / area2
    return np.mean(ratio), rms(ratio), ratio
        
def calcRise(x,y):
    """Calc rise time of pulses"""
    rise = np.zeros( len(y[:,0]) )
    rising = True
    f = positiveCheck(y)
    if f:
        rising = False

    for i, event in enumerate(y):
        if f:
            m = max(event)
        else:
            m = min(event)
        m_index = np.where(event == m)[0][0]
        lo_thresh = m*0.1
        hi_thresh = m*0.9
        y_reverse = event[m_index+1:0:-1]
        try:
            low = interpolateThreshold(x[:m_index+1], y_reverse, lo_thresh, rise=rising)
            high = interpolateThreshold(x[:m_index+1], y_reverse, hi_thresh, rise=rising)
        except:
            continue
        rise[i] = low - high
    return np.mean(rise), rms(rise), rise

def calcFall(x,y):
    """Calc fall time of pulses"""
    fall = np.zeros( len(y[:,0]) )
    rising = True
    f = positiveCheck(y)
    if f:
        rising = False

    for i, event in enumerate(y):
        if f:
            m = max(event)
        else:
            m = min(event)
        m_index = np.where(event == m)[0][0]
        lo_thresh = m*0.1
        hi_thresh = m*0.9
        try:
            low = interpolateThreshold(x[m_index-1:], event[m_index-1:], lo_thresh, rise=rising)
            high = interpolateThreshold(x[m_index-1:], event[m_index-1:], hi_thresh, rise=rising)
        except:
            continue
        fall[i] = low - high
    return np.mean(fall), rms(fall), fall
        
def calcWidth(x,y):
    """Calc width of pulses"""
    width = np.zeros( len(y[:,0]) )
    f = positiveCheck(y)
    if f == True:
        flag_first = True
        flag_second = False
    else:
        flag_first = False
        flag_second = True

    for i, event in enumerate(y):
        if f:
            m = max(event)
        else:
            m = min(event)
        m_index = np.where(event == m)[0][0]
        thresh = m*0.5
        try:
            first = interpolateThreshold(x[:m_index+1], event[:m_index+1], thresh, rise=flag_first)
            second = interpolateThreshold(x[m_index-1:], event[m_index-1:], thresh, rise=flag_second)
        except:
            continue
        width[i] = second - first
    return np.mean(width), rms(width), width

def calcPeak(x,y):
    """Calc min amplitude of pulses"""
    peak = np.zeros( len(y[:,0]) )
    f = positiveCheck(y)
    if f == True:
        for i in range(len(y[:,0])):
            peak[i] = max(y[i,:])
    else:
        for i in range(len(y[:,0])):
            peak[i] = min(y[i,:])
    return np.mean(peak), rms(peak), peak

def calcSNR(x,y,nSamples=50):
    """Calc the signal-to-noise ratio of a set of pulses"""
    snr = np.zeros(len(y[:,0]))
    f = positiveCheck(y)
    if f == True:
        for i in range(len(y[:,0])):
            snr[i] = np.max(y[i,:]) / rms(y[i,:nSamples]) 
    else:
        for i in range(len(y[:,0])):
            snr[i] = np.abs( np.min(y[i,:]) / rms(y[i,:nSamples]) )
    return np.mean(snr), snr

def calcSinglePeak(pos_check, y_arr):
    """Calculate peak values for single trace inputs can be positive or negative."""
    if pos_check == True:
        m = max(y_arr)
    else:
        m = min(y_arr)
    return m

def calcJitter(x1, y1, x2, y2, threshold=0.1):
    """Calc jitter between trig and signal using CFD"""
    p1 = positiveCheck(y1)
    p2 = positiveCheck(y2)
    times = np.zeros(len(y1[:,0]))
    for i in range(len(y1[:,0])):
        m1 = calcSinglePeak(p1, y1[i,:])
        m2 = calcSinglePeak(p2, y2[i,:])
        time_1 = interpolateThreshold(x1, y1[i,:], threshold*m1, rise=p1)
        time_2 = interpolateThreshold(x2, y2[i,:], threshold*m2, rise=p2)
        times[i] = time_1 - time_2
    return np.mean(times), np.std(times), np.std(times)/np.sqrt(2*len(y1[:,0]))

def calcLeadingEdgeTimestamp(x, y, thresh, peak_indicies, rise=False):
    '''
    Step back from peak to get the most noise-free timestamp
    '''
    if rise:
        m = max(y)
    else:
        m = min(y)
    m_index = np.where(y == m)[0][0]
    y_reverse = y[m_index+1:0:-1]
    try:
        timestamp = pd.interpolate_threshold(x[:m_index+1],
                                             y_reverse,
                                             thresh,
                                             rise=rise)
    except:
        timestamp = -999
    return timestamp


def rootifyXY(x, y, scaling=1, name="", title=""):
    '''Turn xy array into a TH1D'''
    hist = ROOT.TH1D(name, title, len(x), x)
    for i, entry in enumerate(y):
        hist.SetBinContent(i, entry*scaling)
    return hist

def peakFinder(x, y, thresh=-0.07, positive=False, min_deltaT=8., plot=False):
    '''
    A peak finding algorithm
    '''
    # Find all points above / below threashold
    if positive:
        above_thresh = np.where( y > thresh, 1., 0. )
    else:
        above_thresh = np.where( y < thresh, 1., 0. )
    # Differentiate to find first point of each pulse to cross thresh
    diff = np.diff(above_thresh)
    start_indicies = np.where( diff > 0.5 )[0]
    if not start_indicies.any() or (min(y) < -1):
        return []
    stop_indicies = np.where( diff < -0.5)[0]
    peak_indicies = np.zeros( len(start_indicies) )
    # Find the peak at each threshold crossing
    for i in range( len(start_indicies) ):
        if positive:
            peak = max(y[start_indicies[i]:stop_indicies[i]])
        else:
            peak = min(y[start_indicies[i]:stop_indicies[i]])
        peak_indicies[i] = (np.where(y[start_indicies[i]:stop_indicies[i]] == peak )[0][0]
                            + start_indicies[i])
    # Loop through indicies and remove any within min_deltaT of each other
    dx = x[1] - x[0]
    peak_check = True
    while peak_check:
        too_close = np.where( np.diff(peak_indicies)*dx < min_deltaT )[0]
        if too_close.size > 0:
            peak_indicies = np.delete(peak_indicies, too_close[0]+1)
        else:
            peak_check = False

    # Some plotting stuff - probably delete after debugging
    if plot and peak_indicies.size > 1:
        x_select = [x[int(i)] for i in peak_indicies]
        y_select = [y[int(i)] for i in peak_indicies]
        plt.plot(x, y)
        plt.plot(x_select, y_select, 'x')
        plt.axhline(thresh, x[0], x[-1], linestyle="--", linewidth=0.5, alpha=0.5)
        plt.show()
    return peak_indicies
