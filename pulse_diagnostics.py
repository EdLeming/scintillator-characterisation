import ROOT
import numpy as np

colors = [ROOT.kBlack, ROOT.kRed-4, 
          ROOT.kGreen+2, ROOT.kCyan+2,
          ROOT.kBlue, ROOT.kMagenta+1,
          ROOT.kGray, ROOT.kOrange+7]

def positive_check(y):
    if np.abs(max(y[0,:])) > np.abs(min(y[0,:])):
        return True
    else:
        return False

def rms(alist):
    '''Calc rms of 1d array'''
    if len(alist) > 1:
        listsum = sum((i - np.mean(alist))**2 for i in alist)
        return np.sqrt(listsum/(len(alist) - 1.0))
    else:
        #print "More than one item needed to calculate RMS, thus returning 0"
        return 0.

def interpolate_threshold(x, y, thresh, rise=True, start=0):
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

def calcRise(x,y):
    """Calc rise time of pulses"""
    rise = np.zeros( len(y[:,0]) )
    rising = True
    f = positive_check(y)
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
            low = interpolate_threshold(x[:m_index+1], y_reverse, lo_thresh, rise=rising)
            high = interpolate_threshold(x[:m_index+1], y_reverse, hi_thresh, rise=rising)
        except:
            continue
        rise[i] = low - high
    return np.mean(rise), rms(rise), rise

def calcFall(x,y):
    """Calc fall time of pulses"""
    fall = np.zeros( len(y[:,0]) )
    rising = True
    f = positive_check(y)
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
            low = interpolate_threshold(x[m_index-1:], event[m_index-1:], lo_thresh, rise=rising)
            high = interpolate_threshold(x[m_index-1:], event[m_index-1:], hi_thresh, rise=rising)
        except:
            continue
        fall[i] = low - high
    return np.mean(fall), rms(fall), fall
        
def calcWidth(x,y):
    """Calc width of pulses"""
    width = np.zeros( len(y[:,0]) )
    f = positive_check(y)
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
            first = interpolate_threshold(x[:m_index+1], event[:m_index+1], thresh, rise=flag_first)
            second = interpolate_threshold(x[m_index-1:], event[m_index-1:], thresh, rise=flag_second)
        except:
            continue
        width[i] = second - first
    return np.mean(width), rms(width), width

def calcPeak(x,y):
    """Calc min amplitude of pulses"""
    peak = np.zeros( len(y[:,0]) )
    f = positive_check(y)
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
    f = positive_check(y)
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
    p1 = positive_check(y1)
    p2 = positive_check(y2)
    times = np.zeros(len(y1[:,0]))
    for i in range(len(y1[:,0])):
        m1 = calcSinglePeak(p1, y1[i,:])
        m2 = calcSinglePeak(p2, y2[i,:])
        time_1 = interpolate_threshold(x1, y1[i,:], threshold*m1, rise=p1)
        time_2 = interpolate_threshold(x2, y2[i,:], threshold*m2, rise=p2)
        times[i] = time_1 - time_2
    return np.mean(times), np.std(times), np.std(times)/np.sqrt(2*len(y1[:,0]))

def dataCleaning(y):
    """Check for any transients that are clearly saturated
    """
    indices = []
    sig_to_noise, saturated = 0, 0
    f = positive_check(y)
    for i, event in enumerate(y):
        if f == True:
            sn = np.abs( np.min(event) / rms(event[:20]) )
            if sn <= 4:
                sig_to_noise = sig_to_noise + 1
                continue
            if np.where( event == max(event))[0].size > 10:
                saturated = saturated + 1
                continue
        else:
            sn = np.abs( np.min(event) / rms(event[:20]) )
            if sn <= 4:
                sig_to_noise = sig_to_noise + 1
                continue
            if np.where( event == min(event))[0].size > 10:
                saturated = saturated + 1
                continue
        indices.append(i)

    print "# Data cleaning"
    print "# SNR < 4:\t{0:d}".format(sig_to_noise)
    print "# Saturated:\t{0:d}".format(saturated)
    print "# Remaining:\t{0:d}".format(len(indices))
    return np.array( [y[i,:] for i in indices] )
        
if __name__ == "__main__":
    import argparse
    import matplotlib.pyplot as plt
    from FileReader import read_h5
    parser = argparse.ArgumentParser("Script to run diagnostrics on pulse shapes")
    parser.add_argument('infile', type=str,
                        help="File(s) to be read in")
    parser.add_argument('outfile', type=str,
                        help="Name of output (root) file to be generated")
    parser.add_argument('-n', '--no_events', type=int, default=10000,
                        help="Number of events to display")
    args = parser.parse_args()

    # Run root in batch mode
    ROOT.gROOT.SetBatch(True)

    # Make rootfile to save data to
    outfile = ROOT.TFile(args.outfile, "RECREATE")
    
    # Read in data and loop over each save channel
    x,y_dict = read_h5(args.infile, nevents=args.no_events)

    # Make some Canvases to hold results
    if len(y_dict.keys()) > 1:
        rise_can = ROOT.TCanvas("Rise_c", "rise")
        fall_can = ROOT.TCanvas("Fall_c", "fall")
        width_can = ROOT.TCanvas("Width_c", "width")
        integral_can = ROOT.TCanvas("Integral_c", "integral")
        peak_can = ROOT.TCanvas("Peak_c","peak")

        rise_leg = ROOT.TLegend(0.65, 0.55, 0.87, 0.77)
        fall_leg = ROOT.TLegend(0.65, 0.55, 0.87, 0.77)
        width_leg = ROOT.TLegend(0.65, 0.55, 0.87, 0.77)
        integral_leg = ROOT.TLegend(0.13, 0.65, 0.35, 0.87)
        peak_leg =  ROOT.TLegend(0.13, 0.65, 0.35, 0.87)

    average_transients = {}
    histos = [] # cos root sucks a bag of dicks
    for i, key in enumerate(y_dict.keys()):
        print "#############################"
        print "# Ch{0} Summary ".format(key)
        print "#############################"
        y = dataCleaning(y_dict[key])
        
        rise_mean, rise_rms, rise = calcRise(x,y)
        fall_mean, fall_rms, fall = calcFall(x,y)
        width_mean, width_rms, width = calcWidth(x,y)
        integral_mean, integral_rms, integral = calcArea(x,y)
        peak_mean, peak_rms, peak = calcPeak(x,y)

        print "########################"
        print "# Pulse measurements"
        print "########################"
        print "#"
        print "# Rise:\t\t {:.1f} +/- {:.1f} [ns]".format(rise_mean, rise_rms)
        print "# Fall:\t\t {:.1f} +/- {:.1f} [ns]".format(fall_mean, fall_rms)
        print "# Width:\t {:.1f} +/- {:.1f} [ns]".format(width_mean, width_rms)
        print "# Integral:\t {:.1f} +/- {:.1f} [V.ns]".format(integral_mean, integral_rms)
        print "# Pulse Height:\t {:.1f} +/- {:.1f} [V]".format(peak_mean, peak_rms)

        summed_pulse = np.zeros( len(y[1,:]) )
        for j, area in enumerate(integral):
            summed_pulse = summed_pulse + (y[j,:])
        summed_pulse = np.array([summed_pulse / j])
        average_transients["Channel {:d}".format(key)] = summed_pulse
        
        rise_avg_mean, rise_avg_rms, _ = calcRise(x, summed_pulse)
        fall_avg_mean, fall_avg_rms, _ = calcFall(x, summed_pulse)
        width_avg_mean, width_avg_rms, _ = calcWidth(x, summed_pulse)
        integral_avg_mean, integral_avg_rms, _ = calcArea(x,y)
        peak_avg_mean, peak_avg_rms, _ = calcPeak(x,y)

        print "#"
        print "# Averaged pulse:"
        print "#"
        print "# Rise:\t\t {:.1f} [ns]".format(rise_avg_mean)
        print "# Fall:\t\t {:.1f} [ns]".format(fall_avg_mean)
        print "# Width:\t {:.1f} [ns]".format(width_mean, width_rms)
        print "# Integral:\t {:.1f} [V.ns]".format(integral_mean, integral_rms)
        print "# Pulse Height:\t {:.1f} [V]".format(peak_mean, peak_rms)
        print ""
        #print "# Making histos to go in {0}...".format(args.outfile)

        rise_h = ROOT.TH1D("Ch{0}_RiseTime".format(key), "",
                           200,
                           0,
                           50)
        fall_h = ROOT.TH1D("Ch{0}_FallTime".format(key), "",
                           200,
                           0,
                           100)
        width_h = ROOT.TH1D("Ch{0}_Width".format(key), "",
                            200,
                            0,
                            100)
        integral_h = ROOT.TH1D("Ch{0}_Integral".format(key), "",
                               200,
                               -200,
                               0)
        peak_h = ROOT.TH1D("Ch{0}_Peak".format(key), "",
                           100,
                           -10,
                           0)
        rise_h.GetXaxis().SetTitle("Rise time [ns]")
        fall_h.GetXaxis().SetTitle("Fall time [ns]")
        width_h.GetXaxis().SetTitle("Pulse width [ns]")
        integral_h.GetXaxis().SetTitle("Pulse integral [V.ns]")
        peak_h.GetXaxis().SetTitle("Pulse height [V]")

        rise_h.SetLineColor( colors[i] )
        fall_h.SetLineColor( colors[i] )
        width_h.SetLineColor( colors[i] )
        integral_h.SetLineColor( colors[i] )
        peak_h.SetLineColor( colors[i] )

        for j, r in enumerate(rise):
            rise_h.Fill( rise[j] )
            fall_h.Fill( fall[j] )
            width_h.Fill( width[j] )
            integral_h.Fill( integral[j] )
            peak_h.Fill( peak[j] )

        # Save histos to disk
        rise_h.Write()
        fall_h.Write()
        width_h.Write()
        integral_h.Write()
        peak_h.Write()

        # Keep histos in memory
        histos.append(rise_h)
        histos.append(fall_h)
        histos.append(width_h)
        histos.append(integral_h)
        histos.append(peak_h)
        
        # Do we also want to draw on canvases?
        opt = ""
        leg_entry = "Channel {0:d}".format(key)
        if len(y_dict.keys()) > 1:
            if i > 0:
                opt = "SAME"
            rise_can.cd()
            rise_h.GetYaxis().SetRangeUser(0, 1.2*rise_h.GetMaximum())
            rise_h.Draw(opt)
            rise_leg.AddEntry(rise_h, leg_entry)
            
            fall_can.cd()
            fall_h.GetYaxis().SetRangeUser(0, 1.2*fall_h.GetMaximum())            
            fall_h.Draw(opt)
            fall_leg.AddEntry(fall_h, leg_entry)
            
            width_can.cd()
            width_h.GetYaxis().SetRangeUser(0, 1.2*width_h.GetMaximum())            
            width_h.Draw(opt)
            width_leg.AddEntry(width_h, leg_entry)
            
            integral_can.cd()
            #integral_h.GetYaxis().SetRangeUser(0, 1.2*integral_h.GetMaximum())
            integral_h.Draw(opt)
            integral_leg.AddEntry(integral_h, leg_entry)

            peak_can.cd()
            #peak_h.GetYaxis().SetRangeUser(0, 1.2*peak_h.GetMaximum())
            peak_h.Draw(opt)
            peak_leg.AddEntry(peak_h, leg_entry)
            
    # Finish making comparitive canvanses
    if len(y_dict.keys()) > 1:
        rise_can.cd()
        rise_leg.Draw()
        rise_can.Update()
        rise_can.Write()

        fall_can.cd()
        fall_leg.Draw()
        fall_can.Update()
        fall_can.Write()

        width_can.cd()
        width_leg.Draw()
        width_can.Update()
        width_can.Write()

        integral_can.cd()
        integral_can.SetLogy()
        integral_leg.Draw()
        integral_can.Update()
        integral_can.Write()

        peak_can.cd()
        peak_can.SetLogy()
        peak_leg.Draw()
        peak_can.Update()
        peak_can.Write()

    plt.figure()
    for key in average_transients:
        plt.plot(x, average_transients[key][0,:], '.', label=key)
    plt.title("Average pulses".format(key))
    plt.xlabel("Time [ns]")
    plt.ylabel("Arbitrary Units")
    plt.legend()
    plt.show()
