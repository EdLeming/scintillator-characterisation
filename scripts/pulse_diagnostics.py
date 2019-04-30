#!/usr/bin/env python
import ROOT
import numpy as np
import utils.transient_calculations as calc

colors = [ROOT.kBlack, ROOT.kRed-4, 
          ROOT.kGreen+2, ROOT.kCyan+2,
          ROOT.kBlue, ROOT.kMagenta+1,
          ROOT.kGray, ROOT.kOrange+7]

def makeSummedPulseHisto(x, y, nEarlyBins=100):
    """Average all pulses"""
    tmp_hist = ROOT.TH1D("", "", len(x), x[0], x[(len(x)-1)])
    summed_hist = ROOT.TH1D("", "", len(x), x[0], x[(len(x)-1)])
    summed_hist.GetXaxis().SetTitle("Time [ns]")
    summed_hist.GetYaxis().SetTitle("Summed Voltage / {0:.2f} ns".format(x[1]-x[0]))

    for i in range(len(y[:,0])):
        tmp_hist.SetContent( np.array( y[i,:-1], dtype=np.float64) )
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
    f = calc.positiveCheck(y)
    for i, event in enumerate(y):
        if f == True:
            sn = np.abs( np.min(event) / calc.rms(event[:20]) )
            if sn <= 4:
                sig_to_noise = sig_to_noise + 1
                continue
            if np.where( event == max(event))[0].size > 10:
                saturated = saturated + 1
                continue
        else:
            sn = np.abs( np.min(event) / calc.rms(event[:20]) )
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
    import utils.file_reader as file_reader
    parser = argparse.ArgumentParser("Script to run diagnostrics on pulse shapes")
    parser.add_argument('infile', type=str,
                        help="Basepath to file(s) to be read in")
    parser.add_argument('outfile', type=str,
                        help="Name of output (root) file to be generated")
    parser.add_argument('-n', '--no_events', type=int, default=10000,
                        help="Number of events to display")
    parser.add_argument('-r', '--termination_resistance', type=float, default=50.,
                        help="Termination resistance at the oscilloscope [50 ohms]")
    args = parser.parse_args()

    # Run root in batch mode
    ROOT.gROOT.SetBatch(True)

    # Make rootfile to save data to
    outfile = ROOT.TFile(args.outfile, "RECREATE")
    
    # Read in data and loop over each save channel
    myFileReader = file_reader.FileReader('')
    extension = args.infile.split("/")[-1].split(".")[-1]
    if extension == "h5":
        myFileReader = file_reader.Hdf5FileReader(args.infile)
    else:
        myFileReader = file_reader.TraceFileReader(args.infile)
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

        rise_mean, rise_rms, rise = calc.calcRise(x,y)
        fall_mean, fall_rms, fall = calc.calcFall(x,y)
        width_mean, width_rms, width = calc.calcWidth(x,y)
        peak_mean, peak_rms, peak = calc.calcPeak(x,y)
        integral_mean, integral_rms, integral = calc.calcArea(x*1e-9 ,y*(1e12/args.termination_resistance))
        integral_mean, intrgral_rms, integral = np.abs(integral_mean), np.abs(integral_rms), np.abs(integral)
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
        print "# Charge:\t {:.1f} +/- {:.1f} [pC]".format(np.abs(integral_mean), integral_rms)
        print "# Pulse Height:\t {:.1f} +/- {:.1f} [V]".format(peak_mean, peak_rms)

        summed_pulse = np.zeros( len(y[1,:]) )
        for j, area in enumerate(integral):
            summed_pulse = summed_pulse + (y[j,:])
        summed_pulse = np.array([summed_pulse / j])
        average_transients["Channel {:d}".format(key)] = summed_pulse
        
        rise_avg_mean, rise_avg_rms, _ = calc.calcRise(x, summed_pulse)
        fall_avg_mean, fall_avg_rms, _ = calc.calcFall(x, summed_pulse)
        width_avg_mean, width_avg_rms, _ = calc.calcWidth(x, summed_pulse)
        integral_avg_mean, integral_avg_rms, _ = calc.calcArea(x*1e-9,
                                                               summed_pulse*(1e12/args.termination_resistance))
        integral_avg_mean = np.abs(integral_avg_mean)
        peak_avg_mean, peak_avg_rms, _ = calc.calcPeak(x,summed_pulse)

        print "#"
        print "# Averaged pulse:"
        print "#"
        print "# Rise:\t\t {:.1f} [ns]".format(rise_avg_mean)
        print "# Fall:\t\t {:.1f} [ns]".format(fall_avg_mean)
        print "# Width:\t {:.1f} [ns]".format(width_avg_mean)
        print "# Charge:\t {:.1f} [pC]".format(np.abs(integral_avg_mean))
        print "# Pulse Height:\t {:.1f} [V]".format(peak_avg_mean)
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
        integral_h = ROOT.TH1D("Ch{0}_Charge".format(key),
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
        integral_h.GetXaxis().SetTitle("Pulse integral [pC]")
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
        ratio_mean, ratio_rms, ratio = calc.calcPulseRatios(x,
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
    #plt.show()
    
    print "Results written to: {0}".format(args.outfile)
