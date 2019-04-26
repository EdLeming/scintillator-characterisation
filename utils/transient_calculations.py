import ROOT
import numpy as np
import matplotlib.pyplot as plt

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

def calcArea(x, y):
    """Calc area of pulses"""
    integral = np.zeros( len(y[:,0]) )
    for i in range(len(y[:,0])):
        integral[i] = np.trapz(y[i,:],x)
    return np.mean(integral), rms(integral), integral
 
def calcPartialArea(x,y,window=20,threshold=-0.05,areaCut=1,early=True):
    '''Calc area between two times, can be used for either late or early light'''
    integral = []
    f = positiveCheck(y)
    valid_events = 0
    #Find area threshold by sorting pulses in data set and choosing only a certain fraction
    no_events=len(y[:,0])
    area = np.zeros(no_events)
    for i, event in enumerate(y):
        area[i] = np.trapz(event,x)
    area.sort()
    if f:
        index = int(no_events*(1-areaCut))-1
    else:
        index = int(no_events*areaCut)-1
    areaThresh = area[index]
    #Event loop
    for i, event in enumerate(y):
        eventArea = np.trapz(event,x)
        #find peak and cut on area threshold
        if f:
            m = max(event)
            if eventArea < areaThresh:
                continue
        else: 
            m = min(event)
            if eventArea > areaThresh:
                continue
        #Place the time window for integration and save integral
        time_start = interpolateThreshold(x, event, threshold, rise=f)
        if early:
            time_1 = time_start
            time_2 = time_start + window
        else:
            time_1 = time_start + window
            time_2 = time_start + 200
        if time_2 > 500:
            time_2 = 500
        if time_1 >500:
            continue
        index_1 = np.where(x > time_1)[0][0]
        index_2 = np.where(x > time_2)[0][0]
        y_section = event[index_1:index_2]
        x_section = x[index_1:index_2]
        integral.append(np.trapz(y_section,x_section))
        valid_events = valid_events + 1
        if valid_events % 2000 == 0:
            print "{0} events processed".format(valid_events)
    
    return integral    

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
        except Exception as e:
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


def rootifyXY(x, y, scaling=1, name="", title=""):
    '''Turn xy array into a TH1D'''
    hist = ROOT.TH1D(name, title, len(x), x)
    for i, entry in enumerate(y):
        hist.SetBinContent(i, entry*scaling)
    return hist

def peakFinder(x, y, thresh=-0.075, positive=False, min_deltaT=10., plot=False):
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
    if not start_indicies.any():
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
    
    # Convert peaks to ints
    peak_indicies = [int(peak) for peak in peak_indicies]
    
    # Some plotting stuff - probably delete after debugging
    if plot:
        x_select = [x[i] for i in peak_indicies]
        y_select = [y[i] for i in peak_indicies]
        plt.plot(x, y)
        plt.plot(x_select, y_select, 'x')
        plt.axhline(thresh, 0, 1, linestyle="--", linewidth=0.5, alpha=0.5)
        plt.show()
    return peak_indicies

def calcLeadingEdgeTimestamp(x, y, peak_index, thresh, positive=False, plot=False):
    '''
    Step back from peak to minimise noise contributions
    '''
    rise=True
    if positive:
        rise=False
    x_select = x[:peak_index+1]
    y_reverse = y[peak_index+1:0:-1]
    pick_off = interpolateThreshold(x_select,
                                    y_reverse,
                                    thresh,
                                    rise=rise)

    timestamp = x[peak_index+1] - pick_off
    if plot:
        plt.plot(x, y)
        plt.axvline(timestamp, 0, 1, linestyle="--", linewidth=0.5, alpha=0.5)
        plt.axhline(thresh, 0, 1, linestyle="--", linewidth=0.5, alpha=0.5)
        plt.show()
    return timestamp

def calcEarlyLateRatio(x, event, threshold=-0.05, window=15, full_width=200, rise=False):
    """Calc early late ratio for single transient"""
    start_time = interpolateThreshold(x, event, threshold, rise=rise)
    window_edge = start_time + window
    end_time = start_time + full_width
    if end_time > 500:
        end_time = 500
        
    start_index = np.where(x > start_time)[0][0]
    window_edge_index = np.where(x > window_edge)[0][0]
    end_index = np.where(x > end_time)[0][0]
    
    y_early = event[start_index:window_edge_index]
    y_late = event[window_edge_index:end_index]
    
    x_early = x[start_index:window_edge_index]
    x_late = x[window_edge_index:end_index]
    
    early_integral = np.trapz(y_early,x_early)
    late_integral = np.trapz(y_late, x_late)
    
    return float(early_integral / late_integral)

def thresholdFinder(x,y,early_window=[0, 150]):
    f = positiveCheck(y)
    noise_h = ROOT.TH1D("noise_h","noise",50,-1,1)
    #Get time stamps for first section of trace, should contain only noise
    index_1 = np.where(x>early_window[0])[0][0]
    index_2 = np.where(x>early_window[1])[0][0]
    #Create histogram of the maximum noise
    for i, event in enumerate(y):
        y_section = event[index_1:index_2]
        noise = min(y_section)
        noise_h.Fill(noise)
    #Set threshold at 3 sigma above the mean noise
    threshold = noise_h.GetMean()
    std = noise_h.GetStdDev()
    if f:
        std = std
    else:
        std = -std
    threshold = threshold + std*3
    return threshold
