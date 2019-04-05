import ROOT
import numpy as np
import utils.transient_calculations as calc
import utils.digital_filters as digital_filters

def calcCharge(x, y, termination=50.):
    '''Calculate the 
    '''
    yc = (y - np.mean(y[:100])) / termination
    return np.abs(np.trapz(yc,x*1e-9))

if __name__ == "__main__":
    import sys
    import argparse
    import matplotlib.pyplot as plt
    import utils.file_reader as file_reader
    parser = argparse.ArgumentParser("Script to run diagnostrics on pulse shapes")
    parser.add_argument('infile', type=str,
                        help="File(s) to be read in")
    parser.add_argument('outfile', type=str,
                        help="Name of output (root) file to be generated")
    parser.add_argument('-n', '--no_events', type=int, default=10000,
                        help="Number of events to display")
    parser.add_argument('-c', '--charge_channel', type=int, default=1,
                        help="Scope channel with large NEMO tube")
    parser.add_argument('-s', '--signal_channel', type=int, default=2,
                        help="Scope channel with the single photon trace")
    parser.add_argument('-t', '--trigger_channel', type=int, default=3,
                        help="Scope channel with trigger channel from scintillating fibre")
    parser.add_argument('-m', '--max_threshold', type=float, default=-0.05,
                        help="Highest threshold to be used on the trigger signal")
    args = parser.parse_args()

    # Run root in batch mode
    ROOT.gROOT.SetBatch(True)
    
    ##########################
    # Read in data and loop over each save channel
    myFileReader = file_reader.FileReader('')
    extension = args.infile.split("/")[-1].split(".")[-1]
    if extension == "h5":
        myFileReader = file_reader.Hdf5FileReader(args.infile)
    else:
        myFileReader = file_reader.TraceFileReader(args.infile)
    x, y_dict = myFileReader.get_xy_data(nevents=args.no_events)

    ################
    no_events = len(y_dict[y_dict.keys()[0]][:,0]) - 1
    dx = x[1] - x[0]
    fs = 1./ (dx*1e-9)
    
    ##########################
    # Check we can use the data in this file
    keys = sorted(y_dict.keys())
    if len(keys) < 3:
        print "Only found {0:d} active channels in: {1}".format(len(y_dict.keys()), args.infile)
        sys.exit(0)
    try:
        charge_traces = y_dict[args.charge_channel]
    except KeyError as e:
        print "Couldn't find channel {0:d} in the passed data set".format(args.charge_channel)
        sys.exit(0)
    try:
        signal_traces = y_dict[args.signal_channel]
    except KeyError as e:
        print "Couldn't find channel {0:d} in the passed data set".format(args.signal_channel)
        sys.exit(0)
    try:
        trigger_traces = y_dict[args.trigger_channel]
    except KeyError as e:
        print "Couldn't find channel {0:d} in the passed data set".format(args.trigger_channel)
        sys.exit(0)

    ##########################
    # Get an initial charge measurement for setting histogram limits
    check_events = 10000
    if no_events < check_events:
        check_events = no_events
    NEMO_charge, trigger_charge = np.zeros(check_events), np.zeros(check_events)
    for i in range(check_events):
        NEMO_charge[i] = calcCharge(x, charge_traces[i,:])
        trigger_charge[i] = calcCharge(x, trigger_traces[i,:])        
    NEMO_charge_end = round( (np.median(sorted(NEMO_charge)[-100:])*1e10), 1)*1e2       #pC
    trigger_charge_end = round( (np.median(sorted(trigger_charge)[-100:])*1e10), 1)*1e2 #pC
    charge_cuts = np.arange(20, NEMO_charge_end+10, NEMO_charge_end*0.3)*1e-12
    
    ##########################
    # Define a range of constant fraction thresholds to apply on the trigger channel
    thresholds = np.arange(0.1, 0.4, 0.1)
    
    #######################
    # Book some histograms
    mean_NEMO_h = ROOT.TH1D("AverageNEMOPulse", "", len(x), x[0], x[(len(x)-1)])
    mean_NEMO_h.GetXaxis().SetTitle("Time [ns]")
    mean_NEMO_h.GetYaxis().SetTitle("Average Voltage / {0:.2f} ns".format(dx))

    mean_trigger_h = ROOT.TH1D("AverageTriggerPulse", "", len(x), x[0], x[(len(x)-1)])
    mean_trigger_h.GetXaxis().SetTitle("Time [ns]")
    mean_trigger_h.GetYaxis().SetTitle("Average Voltage / {0:.2f} ns".format(dx))

    
    NEMO_charge_h = ROOT.TH1D("NemoCharge", "", 200, 0, NEMO_charge_end*1.3)
    NEMO_charge_h.GetXaxis().SetTitle("Pulse integral [pC]")

    trigger_charge_h = ROOT.TH1D("TriggerCharge", "", 200, 0, trigger_charge_end*1.3)
    trigger_charge_h.GetXaxis().SetTitle("Pulse integral [pC]")

    trigger_vs_NEMO_charge_h = ROOT.TH2D("TriggerVsNEMOCharge", "",
                                         50, 0, trigger_charge_end*1.3,
                                         50, 0, NEMO_charge_end*1.3)
    trigger_vs_NEMO_charge_h.GetXaxis().SetTitle("Trigger Tube [pC]")
    trigger_vs_NEMO_charge_h.GetYaxis().SetTitle("NEMO [pC]")

    
    histos, basename = {}, ""
    for i, charge_cut in enumerate(charge_cuts):
        thresh_histos = []
        if i < len(charge_cuts)-1:
            basename = "Charge{:.0f}-{:.0f}".format(charge_cut*1e12, charge_cuts[i+1]*1e12)
        else:
            basename = "Charge{:.0f}".format(charge_cut*1e12)
        for j, thresh in enumerate(thresholds):
            thresh_histos.append(ROOT.TH1D("{0}_t{1:d}".format(basename, int(thresh*100)),
                                           "Charge Cut {:.1f} pC, trigger signal threshold {:d} %".format(charge_cut*1e12, int(thresh*100)),
                                           int((600./dx)), -100, 500))
            thresh_histos[-1].GetXaxis().SetTitle("Pulse separation [ns]")
        histos[charge_cut] = thresh_histos
        
    ##########################
    # Loop over each pulse set
    # and apply some thresholds
    counter = 0
    trigger_h = ROOT.TH1D("", "", len(x), x[0], x[(len(x)-1)])
    NEMO_h = ROOT.TH1D("", "", len(x), x[0], x[(len(x)-1)])
    trigger_times = np.zeros(len(thresholds))
    for i in range(no_events):
        if i % 5000 == 0 and i > 0:
            print "Evaluated {:d} pulse pairs, {:d} had signals".format(i, counter)
            
        # Clean single photon signal and find peaks
        signal_clean = digital_filters.butter_lowpass_filter(signal_traces[i,:], fs, cutoff=500e6)
        try:
            signal_peaks = calc.peakFinder(x, signal_clean, thresh=-0.075, min_deltaT=8.)
        except IndexError as e:
            print "Event {0:d}: Index error {1}".format(i, e)
            continue
        except ValueError as e:
            print "Event {0:d}: Index error {1}".format(i, e)
            continue
        if not signal_peaks: # If the fast 'signal' pmt didn't see anything, continue
            continue
        
            # Clean trigger signal and find peaks
        trigger_clean = digital_filters.butter_lowpass_filter(trigger_traces[i,:], fs, cutoff=350e6)
        try:
            trigger_peaks = calc.peakFinder(x, trigger_clean, thresh=-0.075, min_deltaT=8.)
        except IndexError as e:
            print "Event {0:d}: Index error {1}".format(i, e)
            continue
        except ValueError as e:
            print "Event {0:d}: Index error {1}".format(i, e)
            continue
        if not len(trigger_peaks)==1: # If the trigger pmt didn't see anything, continue
            continue
        
        # Calc timestamps for signal events
        idx_error = False
        signal_times = []
        for peak in signal_peaks:
            thresh = signal_clean[peak]*0.4
            try:
                signal_times.append(calc.calcLeadingEdgeTimestamp(x, signal_clean, peak, thresh))
            except IndexError:
                idx_error = True
                break
        if idx_error:
            continue

        # Calc timestamps for trigger events
        for j, thresh in enumerate(thresholds):
            constant_fraction = trigger_clean[peak]*thresh
            try:
                trigger_times[j] = calc.calcLeadingEdgeTimestamp(x,
                                                                 trigger_clean,
                                                                 trigger_peaks[0],
                                                                 constant_fraction)
            except IndexError:
                idx_error = True
                break
        if idx_error:
            continue

        # Clean NEMO pulse and find charge for binning results
        charge_clean = digital_filters.butter_lowpass_filter(charge_traces[i,:], fs, cutoff=250e6)
        Q = calcCharge(x, charge_clean)
        trigger_Q = calcCharge(x, trigger_clean)
        
        # Add entries into histograms
        for j, cut in enumerate(charge_cuts):
            accept = False
            if j < len(charge_cuts)-1:
                if Q >= cut and Q < charge_cuts[j+1]:
                    accept = True
            else:
                if Q >= cut:
                    accept = True
            if accept:
                for k, trig in enumerate(trigger_times):
                    for sig in signal_times:
                        t = sig - trig
                        histos[cut][k].Fill(t)

        # Add event trigger characteristics into histograms
        NEMO_charge_h.Fill(Q*1e12)
        trigger_charge_h.Fill(trigger_Q*1e12)
        trigger_vs_NEMO_charge_h.Fill(trigger_Q*1e12, Q*1e12)

        NEMO_h.SetContent( np.array( charge_clean, dtype=np.float64) )
        trigger_h.SetContent( np.array( trigger_clean, dtype=np.float64) )

        mean_NEMO_h.Add(NEMO_h)
        mean_trigger_h.Add(trigger_h)

        counter = counter + 1 

    print "Evaluated {:d} pulse pairs, {:d} had signals".format(i, counter)

    # Make rootfile to save data to
    outfile = ROOT.TFile(args.outfile, "RECREATE")
    
    # Write histograms to file
    NEMO_charge_h.Write()
    trigger_charge_h.Write()    
    trigger_vs_NEMO_charge_h.Write()

    mean_NEMO_h.Scale( 1. / counter )
    mean_trigger_h.Scale( 1. / counter )
    mean_NEMO_h.Write()
    mean_trigger_h.Write()

    for cut in sorted(histos.keys()):
        for hist in histos[cut]:
            hist.Write()

    print "Results written to: {0}".format(args.outfile)
