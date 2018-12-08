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
        #print "More than one item needed to calculate RMS, returning 0"
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

def rootify_xy(x, y, scaling=1, name="", title=""):
    '''Turn xy array into a TH1D'''
    hist = ROOT.TH1D(name, title, len(x), x)
    for i, entry in enumerate(y):
        hist.SetBinContent(i, entry*scaling)
    return hist

def makeSummedPulseHisto(x, y, nEarlyBins=100):
    """Average all pulses"""
    tmp_hist = ROOT.TH1D("", "", len(x), x[0], x[(len(x)-1)])
    summed_hist = ROOT.TH1D("", "", len(x), x[0], x[(len(x)-1)])
    summed_hist.GetXaxis().SetTitle("Time [ns]")
    summed_hist.GetYaxis().SetTitle("Summed Voltage / {0:.2f} ns".format(x[1]-x[0]))

    for i in range(len(y[:,0])):
        tmp_hist.SetContent(y[i,:-1])
        summed_hist.Add(tmp_hist)
        tmp_hist.Clear()

    # Apply offset correction using first 100 bins
    tmp_sum = 0
    for i in range(nEarlyBins):
        tmp_sum = tmp_sum + summed_hist.GetBinContent(i)
    for i in range(summed_hist.GetNbinsX()):
        summed_hist.SetBinContent(i, summed_hist.GetBinContent(i)-(tmp_sum/nEarlyBins))

    return summed_hist

def normaliseSummedPulseHisto(histo):
    """Perform a sensible normalisation of the summed pulses
    """
    tmp = ROOT.TH1D(histo)
    tmp.Scale(1./ tmp.Integral("width"))
    lastBin = tmp.FindLastBinAbove(0.001)
    tmp.Delete()

    dx = histo.GetBinCenter(1) - histo.GetBinCenter(0)
    normalised = ROOT.TH1D("","",lastBin, 0, histo.GetBinLowEdge(lastBin))
    normalised.GetYaxis().SetTitle("Probabilitiy / {:.2} ns".format(dx))
    for i in range(lastBin):
        normalised.SetBinContent(i, histo.GetBinContent(i))

    normalised.Scale( 1. / normalised.Integral("width") )
    return normalised

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
    import FileReader as fr
    parser = argparse.ArgumentParser("Script to run diagnostrics on pulse shapes")
    parser.add_argument('infile', type=str,
                        help="Basepath to file(s) to be read in")
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
    myFileReader = fr.FileReader('')
    extension = args.infile.split("/")[-1].split(".")[-1]
    if extension == "h5":
        myFileReader = fr.Hdf5FileReader(args.infile)
    else:
        myFileReader = fr.TraceFileReader(args.infile)
    x, y_dict = myFileReader.get_xy_data(nevents=args.no_events)
    
    # Make some Canvases to hold results
    if len(y_dict.keys()) > 1:
        rise_can = ROOT.TCanvas("Rise_c", "rise")
        fall_can = ROOT.TCanvas("Fall_c", "fall")
        width_can = ROOT.TCanvas("Width_c", "width")
        integral_can = ROOT.TCanvas("Integral_c", "integral")
        peak_can = ROOT.TCanvas("Peak_c","peak")
        ratio_can = ROOT.TCanvas("ratio_c", "ratio")
        average_pulses_can = ROOT.TCanvas("avg_pulse_c", "avg_pulse")
        normalised_pulses_can = ROOT.TCanvas("norm_pulse_c", "norm_pulse")
        
        rise_leg = ROOT.TLegend(0.65, 0.55, 0.87, 0.77)
        fall_leg = ROOT.TLegend(0.65, 0.55, 0.87, 0.77)
        width_leg = ROOT.TLegend(0.65, 0.55, 0.87, 0.77)
        integral_leg = ROOT.TLegend(0.13, 0.65, 0.35, 0.87)
        peak_leg =  ROOT.TLegend(0.13, 0.65, 0.35, 0.87)
        ratio_leg = ROOT.TLegend(0.65, 0.55, 0.87, 0.77)
        average_pulses_leg = ROOT.TLegend(0.65, 0.55, 0.87, 0.7)
        normalised_pulses_leg = ROOT.TLegend(0.65, 0.55, 0.87, 0.7)
        
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
        summed_pulses_h = makeSummedPulseHisto(x,y)

        average_pulses_h = ROOT.TH1D(summed_pulses_h)
        average_pulses_h.Scale( 1. / len(rise) )
        average_pulses_h.GetYaxis().SetTitle("Average Voltage / {0:.2f} ns".format(x[1]-x[0]))
        
        normalised_pulses_h = normaliseSummedPulseHisto(summed_pulses_h)
        
        average_pulses_h.SetName("Ch{0}_Average_pulse".format(key))
        normalised_pulses_h.SetName("Ch{0}_Normalised_pulse".format(key))

        
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

        rise_h = ROOT.TH1D("Ch{0}_RiseTime".format(key),
                           "",
                           200,
                           0,
                           50)
        fall_h = ROOT.TH1D("Ch{0}_FallTime".format(key),
                           "",
                           200,
                           0,
                           100)
        width_h = ROOT.TH1D("Ch{0}_Width".format(key),
                            "",
                            200,
                            0,
                            100)
        integral_h = ROOT.TH1D("Ch{0}_Integral".format(key),
                               "",
                               200,
                               integral_mean+(10*integral_rms),
                               0)
        peak_h = ROOT.TH1D("Ch{0}_Peak".format(key),
                           "",
                           100,
                           peak_mean+(10*peak_rms),
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
        average_pulses_h.SetLineColor( colors[i] )
        normalised_pulses_h.SetLineColor( colors[i] )
        
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
        average_pulses_h.Write()
        normalised_pulses_h.Write()
        
        # Keep histos in memory
        histos.append(rise_h)
        histos.append(fall_h)
        histos.append(width_h)
        histos.append(integral_h)
        histos.append(peak_h)
        histos.append(average_pulses_h)
        histos.append(normalised_pulses_h)
        
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

            average_pulses_can.cd()
            average_pulses_h.Draw(opt)
            average_pulses_leg.AddEntry(average_pulses_h, leg_entry)

            normalised_pulses_can.cd()
            normalised_pulses_h.Draw(opt)
            normalised_pulses_leg.AddEntry(normalised_pulses_h, leg_entry)
            
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

        average_pulses_can.cd()
        average_pulses_leg.Draw()
        average_pulses_can.Update()
        average_pulses_can.Write()

        normalised_pulses_can.cd()
        normalised_pulses_leg.Draw()
        normalised_pulses_can.Update()
        normalised_pulses_can.Write()
        
        # Make ratio plots
        keys = y_dict.keys()
        ratio_mean, ratio_rms, ratio = calcPulseRatios(x,
                                                       y_dict[keys[0]],
                                                       y_dict[keys[1]])
        print "########################"
        print "# Pulse charge ratio"
        print "# Chan {} / Chan {}".format(keys[0], keys[1])
        print "########################"
        print "#"
        print "# Charge ratio:\t {:.2f} +/- {:.2f}".format(ratio_mean, ratio_rms)

        ratio_h = ROOT.TH1D("charge_ratio", "", 50, 0, 3)
        ratio_h.GetXaxis().SetTitle("Charge ratio")
        for r in ratio:
            ratio_h.Fill( r )
        ratio_h.Write()

    plt.figure()
    for key in average_transients:
        plt.plot(x, average_transients[key][0,:], '.', label=key)
    plt.title("Average pulses".format(key))
    plt.xlabel("Time [ns]")
    plt.ylabel("Arbitrary Units")
    plt.legend()
    plt.show()
