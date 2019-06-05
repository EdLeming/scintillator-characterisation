#!/usr/bin/env python
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
import numpy as np
import matplotlib.pyplot as plt

import config
import utils.transient_calculations as calc
import utils.digital_filters as digital_filters


if __name__ == "__main__":
    import argparse
    import time
    import sys
    import utils.file_reader as file_reader
    parser = argparse.ArgumentParser("")
    parser.add_argument('infile', type=str,
                        help="File(s) to be read in")
    parser.add_argument('outfile', type=str,
                        help="Name of output (root) file to be generated")
    parser.add_argument('-n', '--no_events', type=int, default=10000,
                        help="Number of events to display")
    parser.add_argument('-s', '--signal_channel', type=int, default=2,
                        help="Scope channel with the single photon trace")
    parser.add_argument('-t', '--trigger_channel', type=int, default=3,
                        help="Scope channel with high charge trigger signal")
    parser.add_argument("--plot", action='store_true', 
                        help="Plot multi-peak events")
    args = parser.parse_args()    

    # Run root in batch mode
    ROOT.gROOT.SetBatch(True)
    
    # Read in data
    myFileReader = file_reader.FileReader('')
    extension = args.infile.split("/")[-1].split(".")[-1]
    if extension == "h5":
        myFileReader = file_reader.Hdf5FileReader(args.infile)
    else:
        myFileReader = file_reader.TraceFileReader(args.infile)
    x, y_dict = myFileReader.get_xy_data(nevents=args.no_events)
    
    ################
    no_events = len(y_dict[y_dict.keys()[0]][:,0]) -1
    dx = x[1] - x[0]
    fs = 1./ (dx*1e-9)
    
    ##########################
    # Check we can use the data in this file
    keys = sorted(y_dict.keys())
    if len(keys) < 2:
        print "Only found {0:d} active channels in: {1}".format(len(y_dict.keys()), args.infile)
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
    print "Calculating appropriate charge range for histos..."
    check_events = 10000
    if no_events < check_events:
        check_events = no_events
    trigger_charge, signal_charge = np.zeros(check_events), np.zeros(check_events)
    for i in range(check_events):
        trigger_clean = digital_filters.butter_lowpass_filter(trigger_traces[i,:], fs,
                                                              cutoff=config.trigger_BW)
        trigger_charge[i] = calc.calcWindowedCharge(x, trigger_clean,
                                                    thresh=config.trigger_thresh,
                                                    window=config.trigger_window)

        signal_clean = digital_filters.butter_lowpass_filter(signal_traces[i,:], fs,
                                                             cutoff=config.signal_BW)
        signal_charge[i] = calc.calcWindowedCharge(x, signal_clean,
                                                   thresh=config.signal_thresh,
                                                   window=config.signal_window)
        
    trigger_charge_end = round( (np.median(sorted(trigger_charge)[-150:])*1e10), 1)*1e2 #pC
    signal_charge_end = round( (np.median(sorted(signal_charge)[-150:])*1e10), 1)*1e2   #pC

    ##########################
    # Define a range of thresholds to apply on the trigger channel
    cf_thresholds = np.arange(0.05, 0.45, 0.05)
    abs_thresholds = np.arange(-0.01, -0.11, -0.01)
        
    #####################
    # ROOT stuff - book some histos
    mean_signal_h = ROOT.TH1D("AverageSignalPulse", "", len(x), x[0], x[(len(x)-1)])
    mean_signal_h.GetXaxis().SetTitle("Time [ns]")
    mean_signal_h.GetYaxis().SetTitle("Average Voltage / {0:.2f} ns".format(dx))

    mean_trigger_h = ROOT.TH1D("AverageTriggerPulse", "", len(x), x[0], x[(len(x)-1)])
    mean_trigger_h.GetXaxis().SetTitle("Time [ns]")
    mean_trigger_h.GetYaxis().SetTitle("Average Voltage / {0:.2f} ns".format(dx))

    trigger_charge_h = ROOT.TH1D("TriggerCharge", "", 200, 0, trigger_charge_end*1.5)
    trigger_charge_h.GetXaxis().SetTitle("Pulse integral [pC]")

    signal_charge_h = ROOT.TH1D("SignalCharge", "", 200, 0, signal_charge_end*1.5)
    signal_charge_h.GetXaxis().SetTitle("Pulse integral [pC]")
    
    cf_histos, abs_histos = [], []
    cf_ntuple_strings, abs_ntuple_strings = [], []
    for j, fraction in enumerate(cf_thresholds):
        cf_histos.append(ROOT.TH1D("dt_thresh{0:d}".format(int(fraction*100)),
                                   "Coincidence timing resolution: trigger signal threshold {:d} %".format(int(fraction*100)),
                                   int((600./dx)), -100, 500))
        cf_histos[-1].GetXaxis().SetTitle("Pulse separation [ns]")
        cf_ntuple_strings.append("cf_{:.2f}".format(fraction).replace(".","p"))
    # Make fixed threshold histograms
    abs_histos = []
    for j, thresh in enumerate(abs_thresholds):
        abs_histos.append(ROOT.TH1D("dt_thresh{0:.2}V".format(thresh),
                                   "Coincidence timing resolution: trigger signal threshold {:.2f} V".format(thresh),
                                   int((600./dx)), -100, 500))
        abs_histos[-1].GetXaxis().SetTitle("Pulse separation [ns]")
        abs_ntuple_strings.append("led_{:.2f}V".format(np.abs(thresh)).replace(".","p"))

    ##########################
    # Make an NTuple for storing data
    base_ntuple_string = "eventID:TriggerQ:SignalQ"
    cf_string = ":".join(cf_ntuple_strings)
    abs_string = ":".join(abs_ntuple_strings)
    ntuple_string ="{0}:{1}:{2}".format(base_ntuple_string,cf_string,abs_string)
    ntuple = ROOT.TNtuple( 'ntuple', 'ntuple', ntuple_string)
    #################################
    # Some variables for use when looping
    signal_h = ROOT.TH1D("", "", len(x), x[0], x[(len(x)-1)])    
    trigger_h = ROOT.TH1D("", "", len(x), x[0], x[(len(x)-1)])
    cf_trigger_times = np.zeros(len(cf_thresholds))
    abs_trigger_times = np.zeros(len(abs_thresholds))
    # For the ntuple
    ntuple_head = np.zeros(3, dtype=np.float32)
    ntuple_cf_dt = np.zeros(len(cf_thresholds), dtype=np.float32)
    ntuple_abs_dt = np.zeros(len(abs_thresholds), dtype=np.float32)
    ntuple_array = np.zeros( (len(ntuple_head)
                              + len(ntuple_cf_dt)
                              + len(ntuple_abs_dt)), dtype=np.float32)
    counter, start = 0, time.time()
    #################################
    # Loop over events
    for i in range(no_events):
        signal_clean = digital_filters.butter_lowpass_filter(signal_traces[i,:], fs, cutoff=config.signal_BW)
        try:
            peaks = calc.peakFinder(x, signal_clean, thresh=config.signal_thresh, min_deltaT=8., plot=False)
        except IndexError as e:
            print "Event {0:d}: Signal peaks index error {1}".format(i, e)
            continue
        except ValueError as e:
            print "Event {0:d}: Signal peaks value error {1}".format(i, e)
            continue
        if not peaks.any(): # If the fast 'signal' pmt didn't see anything, continue
            continue

        # Calc timestamps for fast tube
        idx_error = False
        fast_times = []
        for peak in peaks:
            thresh = signal_clean[peak]*0.4
            try:
                fast_times.append(calc.calcLeadingEdgeTimestamp(x, signal_clean, peak, thresh))
            except IndexError:
                idx_error = True
                break
        if idx_error:
            continue
            
        # Clean trigger tube signal and find timestamps
        trigger_clean = digital_filters.butter_lowpass_filter(trigger_traces[i,:], fs, cutoff=config.trigger_BW)
        try:
            trigger_peaks = calc.peakFinder(x, trigger_clean, thresh=config.trigger_thresh, min_deltaT=20., plot=False)
        except IndexError as e:
            print "Event {0:d}: Trigger peaks index error {1}".format(i, e)
            continue
        except ValueError as e:
            print "Event {0:d}: Trigger peaks value error {1}".format(i, e)
            continue
        # Only continue is there is a single peak
        #if not len(trigger_peaks) == 1:
        if not trigger_peaks.any():
            continue

        # Calc timestamps for trigger events
        for j, fraction in enumerate(cf_thresholds):
            thresh = trigger_clean[trigger_peaks[0]]*fraction
            try:
                cf_trigger_times[j] = calc.calcLeadingEdgeTimestamp(x,
                                                                    trigger_clean,
                                                                    trigger_peaks[0],
                                                                    thresh,
                                                                    plot=False)
            except IndexError as e:
                idx_error = True
                print "Event {0:d}: Constant fraction evaluation error {1}".format(i, e)
                break
        if idx_error:
            continue

        # Fixed threshold
        for j, thresh in enumerate(abs_thresholds):
            try:
                abs_trigger_times[j] = calc.calcLeadingEdgeTimestamp(x,
                                                                     trigger_clean,
                                                                     trigger_peaks[0],
                                                                     thresh)
            except IndexError as e:
                idx_error = True
                print "Event {0:d}: Fixed threshold evaluation error {1}".format(i, e)
                break
        if idx_error:
            continue
        # Calc charges for this event - these have to be in arrays for the ntuple to get filled correctly
        trigger_Q = calc.calcWindowedCharge(x, trigger_clean,
                                            thresh=config.trigger_thresh,
                                            window=config.trigger_window)
        signal_Q = calc.calcWindowedCharge(x, signal_clean,
                                           thresh=config.signal_thresh,
                                           window=config.signal_window)
        # Fill histograms
        signal_h.SetContent( np.array( signal_clean, dtype=np.float64) )
        trigger_h.SetContent( np.array( trigger_clean, dtype=np.float64) )
        mean_signal_h.Add(signal_h)
        mean_trigger_h.Add(trigger_h)
        signal_charge_h.Fill(signal_Q*1e12)
        trigger_charge_h.Fill(trigger_Q*1e12)
        for j, trig in enumerate(cf_trigger_times):
            for sig in fast_times:
                dt = sig - trig
                cf_histos[j].Fill(dt)
                ntuple_cf_dt[j] = dt
        for j, trig in enumerate(abs_trigger_times):
            for sig in fast_times:
                dt = sig - trig
                abs_histos[j].Fill(dt)
                ntuple_abs_dt[j] = dt
        # Fill the ntuple
        ntuple_head[0] = i
        ntuple_head[1] = trigger_Q*1e12
        ntuple_head[2] = signal_Q*1e12
        np.concatenate((ntuple_head, ntuple_cf_dt, ntuple_abs_dt), out=ntuple_array)
        ntuple.Fill(ntuple_array)
        # Increment counter and print
        counter = counter + 1
        if i % 10000 == 0 and i > 0:
            print "Evaluated {:d} pulse pairs, {:d} had signals".format(i, counter)

    ##############################
    # Save histos and exit
    end = time.time()
    print "{0} events took {1:.2f}s to process [{2:.4f} s/evt]".format(i+1,
                                                                       (end-start),
                                                                       (end-start)/(i+1))
    outfile = ROOT.TFile(args.outfile, "RECREATE")
    mean_signal_h.Scale( 1. / counter )
    mean_trigger_h.Scale( 1. / counter )
    
    mean_signal_h.Write()
    mean_trigger_h.Write()
    trigger_charge_h.Write()
    signal_charge_h.Write()
    for histo in cf_histos:
        histo.Write()
    for histo in abs_histos:
        histo.Write()
    # Write ntuple
    ntuple.Write()
    print "Results written to: {0}".format(args.outfile)
