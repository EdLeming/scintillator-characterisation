import ROOT
import numpy as np
import matplotlib.pyplot as plt

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

    #####################
    # ROOT stuff - book some histos
    ROOT.gROOT.SetBatch(True)
    mean_signal_h = ROOT.TH1D("AverageSignalPulse", "", len(x), x[0], x[(len(x)-1)])
    mean_signal_h.GetXaxis().SetTitle("Time [ns]")
    mean_signal_h.GetYaxis().SetTitle("Average Voltage / {0:.2f} ns".format(dx))

    mean_trigger_h = ROOT.TH1D("AverageTriggerPulse", "", len(x), x[0], x[(len(x)-1)])
    mean_trigger_h.GetXaxis().SetTitle("Time [ns]")
    mean_trigger_h.GetYaxis().SetTitle("Average Voltage / {0:.2f} ns".format(dx))

    dt_h = ROOT.TH1D("Convolved-response",
                     "Time resolution between scintillating fibre trigger and Cerenkov light",
                     int(1000*(1./dx)), -500, 500)
    dt_h.GetXaxis().SetTitle("Pulse separation [ns]")
    
    # Loop over events and calculate timestamps for observed events
    signal_h = ROOT.TH1D("", "", len(x), x[0], x[(len(x)-1)])    
    trigger_h = ROOT.TH1D("", "", len(x), x[0], x[(len(x)-1)])
    start = time.time()
    for i in range(args.no_events):            
        signal_clean = digital_filters.butter_lowpass_filter(signal_traces[i,:], fs, cutoff=500e6)
        try:
            peaks = calc.peakFinder(x, signal_clean, thresh=-0.07, min_deltaT=8.)
        except IndexError as e:
            print "Event {0:d}: Signal peaks index error {1}".format(i, e)
            continue
        except ValueError as e:
            print "Event {0:d}: Signal peaks value error {1}".format(i, e)
            continue
        if not peaks: # If the fast 'signal' pmt didn't see anything, continue
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
        trigger_clean = digital_filters.butter_lowpass_filter(trigger_traces[i,:], fs, cutoff=350e6)
        try:
            trigger_peaks = calc.peakFinder(x, trigger_clean, thresh=-0.1, min_deltaT=8.)
        except IndexError as e:
            print "Event {0:d}: Trigger peaks index error {1}".format(i, e)
            continue
        except ValueError as e:
            print "Event {0:d}: Trigger peaks value error {1}".format(i, e)
            continue
        # Only continue is there is a single peak
        if not len(trigger_peaks) == 1:
            continue
        
        # Calc timestamps for trigger tube
        trigger_time = 0
        thresh = trigger_clean[trigger_peaks[0]]*0.1
        try:
            trigger_time = calc.calcLeadingEdgeTimestamp(x, trigger_clean, trigger_peaks[0], thresh)
        except IndexError as e:
            print "Event {0:d}: Trigger thresh index error {1}".format(i, e)
            continue

        # Fill histograms
        signal_h.SetContent( np.array( signal_clean, dtype=np.float64) )
        trigger_h.SetContent( np.array( trigger_clean, dtype=np.float64) )

        mean_signal_h.Add(signal_h)
        mean_trigger_h.Add(trigger_h)
        for sig in fast_times:
            dt_h.Fill(sig - trigger_time)
            
        if i % 10000 == 0 and i is not 0:
            print "{0} events processed".format(i) 
    end = time.time()
    print "{0} events took {1:.2f}s to process [{2:.4f} s/evt]".format(i+1,
                                                                       (end-start),
                                                                       (end-start)/(i+1))
    outfile = ROOT.TFile(args.outfile, "RECREATE")
    mean_signal_h.Write()
    mean_trigger_h.Write()
    dt_h.Write()
    
    print "Results written to: {0}".format(args.outfile)
