import ROOT
import matplotlib.pyplot as plt

import utils.transient_calculations as calc
import utils.digital_filters as digital_filters


if __name__ == "__main__":
    import argparse
    import time
    import sys
    import utils.file_reader as file_reader
    parser = argparse.ArgumentParser("Calculate the time response of a PMT from LED data")
    parser.add_argument('infile', type=str,
                        help="File(s) to be read in")
    parser.add_argument('outfile', type=str,
                        help="Name of output (root) file to be generated")
    parser.add_argument('-n', '--no_events', type=int, default=10000,
                        help="Number of events to display")
    parser.add_argument('-ne', '--nemo_channel', type=int, default=1,
                        help="Scope channel with high charge trigger signal")
    parser.add_argument('-s', '--signal_channel', type=int, default=2,
                        help="Scope channel with the single photon trace")
    parser.add_argument('-t', '--trigger_channel', type=int, default=3,
                        help="Scope channel with high charge trigger signal")
    parser.add_argument("--plot", action='store_true', 
                        help="Plot multi-peak events")
    args = parser.parse_args()    
    
    # Read in data and loop over each save channel
    myFileReader = file_reader.FileReader('')
    extension = args.infile.split("/")[-1].split(".")[-1]
    if extension == "h5":
        myFileReader = file_reader.Hdf5FileReader(args.infile)
    else:
        myFileReader = file_reader.TraceFileReader(args.infile)
    x, y_dict = myFileReader.get_xy_data(nevents=args.no_events)

    ################
    no_events = len(y_dict[y_dict.keys()[0]][:,0])
    dx = x[1] - x[0]
    fs = 1./ (dx*1e-9)
    
    ##########################
    # Check we can use the data in this file
    keys = sorted(y_dict.keys())
    if len(keys) < 3:
        print "Only found {0:d} active channels in: {1}".format(len(y_dict.keys()), args.infile)
        sys.exit(0)
    try:
        fast = y_dict[args.signal_channel]
    except KeyError as e:
        print "Couldn't find channel {0:d} in the passed data set".format(args.signal_channel)
        sys.exit(0)
    try:
        nemo = y_dict[args.nemo_channel]
    except KeyError as e:
        print "Couldn't find channel {0:d} in the passed data set".format(args.trigger_channel)
        sys.exit(0)
    try:
        trigger = y_dict[args.trigger_channel]
    except KeyError as e:
        print "Couldn't find channel {0:d} in the passed data set".format(args.trigger_channel)
        sys.exit(0)

    #####################
    # ROOT stuff
    ROOT.gROOT.SetBatch(True)
    outfile = ROOT.TFile(args.outfile, "RECREATE")
    dt_h = ROOT.TH1D("Convolved-response", "", int(1000*(1./dx)), -500, 500)
    fast_h = ROOT.TH1D("fast-response", "", int(1000*(1./dx)), -500, 500)
    nemo_h = ROOT.TH1D("nemo-response", "", int(1000*(1./dx)), -500, 500)
    dt_h.GetXaxis().SetTitle("Pulse separation [ns]")
    fast_h.GetXaxis().SetTitle("Pulse separation [ns]")
    nemo_h.GetXaxis().SetTitle("Pulse separation [ns]")

    # Loop over events and calculate timestamps for observed events
    start = time.time()
    for i in range(args.no_events):            
        signal_clean = digital_filters.butter_lowpass_filter(fast[i,:], fs, cutoff=500e6)
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
            
        # Clean nemo tube signal and find timestamps
        nemo_clean = digital_filters.butter_lowpass_filter(nemo[i,:], fs, cutoff=500e6)
        try:
            nemo_peaks = calc.peakFinder(x, nemo_clean, thresh=-0.06, min_deltaT=8.)
        except IndexError as e:
            print "Event {0:d}: Nemo peaks index error {1}".format(i, e)
            continue
        except ValueError as e:
            print "Event {0:d}: Nemo peaks value error {1}".format(i, e)
            continue
        if not nemo_peaks or len(nemo_peaks) > 1: # Only if there's more than one peak in the trigger, ignore
            continue
        
        # Calc timestamps for nemo tube
        nemo_time = 0
        thresh = nemo_clean[nemo_peaks[0]]*0.4
        try:
            nemo_time = calc.calcLeadingEdgeTimestamp(x, nemo_clean, nemo_peaks[0], thresh)
        except IndexError as e:
            print "Event {0:d}: Trigger thresh index error {1}".format(i, e)
            continue

        # Calc timestamp from trigger signal
        trigger_time = calc.interpolateThreshold(x, trigger[i,:], 0.5, rise=False)

        for sig in fast_times:
            dt_h.Fill(sig - nemo_time)
            fast_h.Fill(sig - trigger_time)
        nemo_h.Fill(nemo_time - trigger_time)
        if i % 10000 == 0 and i is not 0:
            print "{0} events processed".format(i) 
    end = time.time()
    print "{0} events took {1:.2f}s to process [{2:.4f} s/evt]".format(i+1,
                                                                       (end-start),
                                                                       (end-start)/(i+1))
    fast_h.Write()
    nemo_h.Write()
    dt_h.Write()
    
    print "Results written to: {0}".format(args.outfile)
