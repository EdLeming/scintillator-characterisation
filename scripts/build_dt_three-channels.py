#!/usr/bin/env python
import time
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
import numpy as np
import utils.transient_calculations as calc
import utils.digital_filters as digital_filters

def calcCharge(x, y, termination=50.):
    '''Calculate the 
    '''
    yc = (y - np.mean(y[:200])) / termination
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
    no_events = len(y_dict[y_dict.keys()[0]][:,0])-1
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
    check_events = 20000
    if no_events < check_events:
        check_events = no_events
    NEMO_charge, trigger_charge, signal_charge = np.zeros(check_events), np.zeros(check_events), np.zeros(check_events)
    for i in range(check_events):
        NEMO_charge[i] = calcCharge(x, charge_traces[i,:])
        trigger_charge[i] = calcCharge(x, trigger_traces[i,:])
        signal_charge[i] = calcCharge(x, signal_traces[i,:])
    NEMO_charge_end = round( (np.median(sorted(NEMO_charge)[-150:])*1e10), 1)*1e2       #pC
    trigger_charge_end = round( (np.median(sorted(trigger_charge)[-150:])*1e10), 1)*1e2 #pC
    signal_charge_end = round( (np.median(sorted(signal_charge)[-150:])*1e10), 1)*1e2   #pC

    # Set an array of charge cuts.
    NEMO_charge_cuts = np.array([50])*1e-12 # Just one at 50 pC - others can be done in post processing

    ##########################
    # Define a range of thresholds to apply on the trigger channel
    cf_thresholds = np.arange(0.05, 0.45, 0.10)
    abs_thresholds = np.arange(-0.01, -0.10, -0.02)

    #######################
    # Book histograms    
    NEMO_charge_h = ROOT.TH1D("NemoCharge", "", 200, 0, NEMO_charge_end*1.3)
    NEMO_charge_h.GetXaxis().SetTitle("Pulse integral [pC]")

    trigger_charge_h = ROOT.TH1D("TriggerCharge", "", 200, 0, trigger_charge_end*1.3)
    trigger_charge_h.GetXaxis().SetTitle("Pulse integral [pC]")

    signal_charge_h = ROOT.TH1D("SignalCharge", "", 200, 0, signal_charge_end*1.3)
    trigger_charge_h.GetXaxis().SetTitle("Pulse integral [pC]")

    trigger_vs_NEMO_charge_h = ROOT.TH2D("TriggerVsNEMOCharge", "",
                                         50, 0, trigger_charge_end*1.3,
                                         50, 0, NEMO_charge_end*1.3)
    trigger_vs_NEMO_charge_h.GetXaxis().SetTitle("Trigger Tube [pC]")
    trigger_vs_NEMO_charge_h.GetYaxis().SetTitle("NEMO [pC]")

    mean_NEMO_h = ROOT.TH1D("AverageNEMOPulse", "", len(x), x[0], x[(len(x)-1)])
    mean_NEMO_h.GetXaxis().SetTitle("Time [ns]")
    mean_NEMO_h.GetYaxis().SetTitle("Average Voltage / {0:.2f} ns".format(dx))

    mean_trigger_h = ROOT.TH1D("AverageTriggerPulse", "", len(x), x[0], x[(len(x)-1)])
    mean_trigger_h.GetXaxis().SetTitle("Time [ns]")
    mean_trigger_h.GetYaxis().SetTitle("Average Voltage / {0:.2f} ns".format(dx))

    mean_signal_h = ROOT.TH1D("AverageSignalPulse", "", len(x), x[0], x[(len(x)-1)])
    mean_signal_h.GetXaxis().SetTitle("Time [ns]")
    mean_signal_h.GetYaxis().SetTitle("Average Voltage / {0:.2f} ns".format(dx))
    
    # Loop over charge cuts and differnt thresholds to pick off
    cf_histos, abs_histos, basename = {}, {}, ""
    cf_ntuple_strings, abs_ntuple_strings = [], []
    for i, charge_cut in enumerate(NEMO_charge_cuts):
        
        if i < len(NEMO_charge_cuts)-1:
            basename = "dt_cut_{:.0f}-{:.0f}pC".format(charge_cut*1e12, NEMO_charge_cuts[i+1]*1e12)
        else:
            basename = "dt_cut_{:.0f}pC".format(charge_cut*1e12)

        # Make constant fraction histograms
        histo_arr = []
        for j, fraction in enumerate(cf_thresholds):
            histo_arr.append(ROOT.TH1D("{0}_t{1:d}".format(basename, int(fraction*100)),
                                        "Charge Cut {:.1f} pC, trigger signal threshold {:d} %".format(charge_cut*1e12, int(fraction*100)),
                                        int((600./dx)), -100, 500))
            histo_arr[-1].GetXaxis().SetTitle("Pulse separation [ns]")
            cf_ntuple_strings.append("cf_{:.2f}".format(fraction).replace(".","p"))
        cf_histos[charge_cut] = histo_arr

        # Make fixed threshold histograms
        histo_arr = []
        for j, thresh in enumerate(abs_thresholds):
            histo_arr.append(ROOT.TH1D("{0}_t{1:.2}V".format(basename, thresh),
                                        "Charge Cut {:.1f} pC, trigger signal threshold {:.2f} V".format(charge_cut*1e12, thresh),
                                        int((600./dx)), -100, 500))
            histo_arr[-1].GetXaxis().SetTitle("Pulse separation [ns]")
            abs_ntuple_strings.append("led_{:.2f}V".format(np.abs(thresh)).replace(".","p"))
        abs_histos[charge_cut] = histo_arr

    # Make charge vs dt histos
    dt_vs_NEMO_charge_cf_h = ROOT.TH2D("NEMOChargeVsdt-cf", "",
                                    int(500), 0, 500,
                                    50, 0, NEMO_charge_end*1.3)
    dt_vs_NEMO_charge_cf_h.GetXaxis().SetTitle("dt [ns]")
    dt_vs_NEMO_charge_cf_h.GetYaxis().SetTitle("NEMO Tube [pC]")

    dt_vs_trigger_charge_cf_h = ROOT.TH2D("TriggerChargeVsdt-cf", "",
                                       int(500), 0, 500,
                                       50, 0, trigger_charge_end*1.3)
    dt_vs_trigger_charge_cf_h.GetXaxis().SetTitle("dt [ns]")
    dt_vs_trigger_charge_cf_h.GetYaxis().SetTitle("Trigger Tube [pC]")

    dt_vs_signal_charge_cf_h = ROOT.TH2D("SignalChargeVsdt-cf", "",
                                       int(500), 0, 500,
                                       50, 0, signal_charge_end*1.3)
    dt_vs_signal_charge_cf_h.GetXaxis().SetTitle("dt [ns]")
    dt_vs_signal_charge_cf_h.GetYaxis().SetTitle("Signal Tube [pC]")

    ##########################
    # Make an NTuple for storing data
    base_ntuple_string = "eventID:NemoQ:TriggerQ:SignalQ"
    cf_string = ":".join(cf_ntuple_strings)
    abs_string = ":".join(abs_ntuple_strings)
    ntuple_string ="{0}:{1}:{2}".format(base_ntuple_string,cf_string,abs_string)
    ntuple = ROOT.TNtuple( 'ntuple', 'ntuple', ntuple_string)
    #################################
    # Define some variables to use when looping
    counter = 0
    start = time.time()
    NEMO_h = ROOT.TH1D("", "", len(x), x[0], x[(len(x)-1)])
    trigger_h = ROOT.TH1D("", "", len(x), x[0], x[(len(x)-1)])
    signal_h = ROOT.TH1D("", "", len(x), x[0], x[(len(x)-1)])
    cf_trigger_times = np.zeros(len(cf_thresholds))
    abs_trigger_times = np.zeros(len(abs_thresholds))
    # For the ntuple
    ntuple_head = np.zeros(4, dtype=np.float32)
    ntuple_cf_dt = np.zeros(len(cf_thresholds), dtype=np.float32)
    ntuple_abs_dt = np.zeros(len(abs_thresholds), dtype=np.float32)
    ntuple_array = np.zeros( (len(ntuple_head)
                              + len(ntuple_cf_dt)
                              + len(ntuple_abs_dt)), dtype=np.float32)
    #################################
    # Loop over events
    for i in range(no_events):
        if i % 5000 == 0 and i > 0:
            print "Evaluated {:d} pulse pairs, {:d} had signals".format(i, counter)
            
        # Clean single photon signal and find peaks
        signal_clean = digital_filters.butter_lowpass_filter(signal_traces[i,:], fs, cutoff=500e6)
        try:
            signal_peaks = calc.peakFinder(x, signal_clean, thresh=-0.075, min_deltaT=8.)
        except IndexError as e:
            print "Event {0:d}: Signal peak finder index error {1}".format(i, e) 
            continue
        except ValueError as e:
            print "Event {0:d}: Signal peak finder index error {1}".format(i, e)
            continue
        if not signal_peaks: # If the fast 'signal' pmt didn't see anything, continue
            continue

        # Clean trigger signal and find peaks
        trigger_clean = digital_filters.butter_lowpass_filter(trigger_traces[i,:], fs, cutoff=500e6)
        try:
            trigger_peaks = calc.peakFinder(x,
                                            trigger_clean,
                                            thresh=-0.5,
                                            min_deltaT=15.)
        except IndexError as e:
            print "Event {0:d}: Trigger peak finder index error {1}".format(i, e)
            continue
        except ValueError as e:
            print "Event {0:d}: Trigger peak findex error {1}".format(i, e)
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

        #####################################
        # Calc timestamps for trigger events
        # Constant fraction
        for j, fraction in enumerate(cf_thresholds):
            thresh = trigger_clean[trigger_peaks[0]]*fraction
            try:
                cf_trigger_times[j] = calc.calcLeadingEdgeTimestamp(x,
                                                                    trigger_clean,
                                                                    trigger_peaks[0],
                                                                    thresh,
                                                                    plot=False)
            except IndexError:
                idx_error = True
                break
        if idx_error:
            continue
        # Fixed threshold
        for j, thresh in enumerate(abs_thresholds):
            try:
                abs_trigger_times[j] = calc.calcLeadingEdgeTimestamp(x,
                                                                     trigger_clean,
                                                                     trigger_peaks[0],
                                                                     thresh,
                                                                     plot=False)
            except IndexError:
                idx_error = True
                break
        if idx_error:
            continue
        
        # Clean NEMO pulse and find charge for binning results
        charge_clean = digital_filters.butter_lowpass_filter(charge_traces[i,:], fs, cutoff=250e6)
        # Calc charges for this event - these have to be in arrays for the ntuple to get filled correctly
        NEMO_Q = calcCharge(x, charge_clean)
        trigger_Q = calcCharge(x, trigger_clean)
        signal_Q = calcCharge(x, signal_clean)
        # Add entries into histograms
        for j, cut in enumerate(NEMO_charge_cuts):
            accept = False
            if j < len(NEMO_charge_cuts)-1:
                if NEMO_Q >= cut and NEMO_Q < NEMO_charge_cuts[j+1]:
                    accept = True
            else:
                if NEMO_Q >= cut:
                    accept = True
            if accept:
                for k, trig in enumerate(cf_trigger_times):
                    for sig in signal_times:
                        dt = sig - trig
                        ntuple_cf_dt[k] = dt
                        cf_histos[cut][k].Fill(dt)
                        
                for k, trig in enumerate(abs_trigger_times):
                    for sig in signal_times:
                        dt = sig - trig
                        ntuple_abs_dt[k] = dt
                        abs_histos[cut][k].Fill(dt)
        # Add event characteristics into histograms
        NEMO_charge_h.Fill(NEMO_Q*1e12)
        trigger_charge_h.Fill(trigger_Q*1e12)
        signal_charge_h.Fill(signal_Q*1e12)
        trigger_vs_NEMO_charge_h.Fill(trigger_Q*1e12, NEMO_Q*1e12)
        NEMO_h.SetContent( np.array( charge_clean, dtype=np.float64) )
        trigger_h.SetContent( np.array( trigger_clean, dtype=np.float64) )
        signal_h.SetContent( np.array( signal_clean, dtype=np.float64) )
        mean_NEMO_h.Add(NEMO_h)
        mean_trigger_h.Add(trigger_h)
        mean_signal_h.Add(signal_h)
        # Fill 2D histos
        for k, trig in enumerate(cf_trigger_times):
            for sig in signal_times:
                dt = sig - trig
                dt_vs_NEMO_charge_cf_h.Fill(dt, NEMO_Q*1e12)
                dt_vs_trigger_charge_cf_h.Fill(dt, trigger_Q*1e12)
                dt_vs_signal_charge_cf_h.Fill(dt, signal_Q*1e12)
        # Fill the ntuple
        ntuple_head[0] = i
        ntuple_head[1] = NEMO_Q*1e12
        ntuple_head[2] = trigger_Q*1e12
        ntuple_head[3] = signal_Q*1e12
        np.concatenate((ntuple_head, ntuple_cf_dt, ntuple_abs_dt), out=ntuple_array)
        ntuple.Fill(ntuple_array)
        # Increment counter
        counter = counter + 1 
    ##############################
    # Save histos and exit
    end = time.time()
    print "Processing {0} events took {1:.2f} mins to process [{2:.4f} s/evt]".format(i+2,
                                                                                      (end-start)/60.,
                                                                                      (end-start)/(i+2))
    
    # Make rootfile to save data to
    outfile = ROOT.TFile(args.outfile, "RECREATE")

    # Normalise mean pulse histos
    mean_NEMO_h.Scale( 1. / counter )
    mean_trigger_h.Scale( 1. / counter )
    mean_signal_h.Scale( 1. / counter )
    # Write histograms to file
    mean_NEMO_h.Write()
    mean_trigger_h.Write()
    mean_signal_h.Write()
    NEMO_charge_h.Write()
    trigger_charge_h.Write()
    signal_charge_h.Write()    
    for cut in sorted(cf_histos.keys()):
        for hist in cf_histos[cut]:
            hist.Write()
    for cut in sorted(abs_histos.keys()):
        for hist in abs_histos[cut]:
            hist.Write()
    trigger_vs_NEMO_charge_h.Write()
    dt_vs_NEMO_charge_cf_h.Write()
    dt_vs_trigger_charge_cf_h.Write()
    dt_vs_signal_charge_cf_h.Write()
    # Write ntuple
    ntuple.Write()
    print "Results written to: {0}".format(args.outfile)
