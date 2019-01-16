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
    parser.add_argument('-s', '--signal_channel', type=int, default=2,
                        help="Scope channel with the single photon trace")
    parser.add_argument('-t', '--trigger_channel', type=int, default=1,
                        help="Scope channel with high charge trigger signal")
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
    no_events = len(y_dict[y_dict.keys()[0]][:,0])
    dx = x[1] - x[0]
    fs = 1./ (dx*1e-9)
    
    ##########################
    # Check we can use the data in this file
    keys = sorted(y_dict.keys())
    if len(keys) < 2:
        print "Only found {0:d} active channels in: {1}".format(len(y_dict.keys()), args.infile)
        sys.exit(0)
    try:
        signal = y_dict[args.signal_channel]
    except KeyError as e:
        print "Couldn't find channel {0:d} in the passed data set".format(args.signal_channel)
        sys.exit(0)
    try:
        trigger = y_dict[args.trigger_channel]
    except KeyError as e:
        print "Couldn't find channel {0:d} in the passed data set".format(args.trigger_channel)
        sys.exit(0)

    ##########################
    # Get a measurement of the charge for thresholding
    check_events = 10000
    if no_events < check_events:
        check_events = no_events
    charge = sorted([calcCharge(x, trigger[i,:]) for i in range(check_events)])
    charge_end = round( (np.median(charge[-100:])*1e10), 1)*1e-10
    charge_cuts = np.arange(0, 700, 100)*1e-12
    
    min_threshold = round(calc.rms( trigger[i,:200] ),2)*-2.
    thresholds = np.linspace(min_threshold, args.max_threshold, 4)

    #######################
    # Book some histograms
    summed_trigger_h = ROOT.TH1D("AverageTriggerPulse", "", len(x), x[0], x[(len(x)-1)])
    summed_trigger_h.GetXaxis().SetTitle("Time [ns]")
    summed_trigger_h.GetYaxis().SetTitle("Average Voltage / {0:.2f} ns".format(dx))

    charge_h = ROOT.TH1D("TriggerCharge", "", 200, 0, charge_end*1.3*1e12)
    charge_h.GetXaxis().SetTitle("Pulse integral [pC]")

    histos, basename = {}, ""
    for i, charge_cut in enumerate(charge_cuts):
        thresh_histos = []
        if i < len(charge_cuts)-1:
            basename = "Charge{:.0f}-{:.0f}".format(charge_cut*1e12, charge_cuts[i+1]*1e12)
        else:
            basename = "Charge{:.0f}".format(charge_cut*1e12)
        for j, thresh in enumerate(thresholds):
            thresh_histos.append(ROOT.TH1D("{}_t{:.3f}".format(basename, thresh),
                                           "Charge Cut {:.1f} pC, trigger signal threshold {:.3f}V".format(charge_cut*1e12, thresh),
                                           int((500./dx)), -400, 100))
            thresh_histos[-1].GetXaxis().SetTitle("Pulse separation [ns]")
        histos[charge_cut] = thresh_histos
        
    ##########################
    # Loop over each pulse set
    # and apply some thresholds
    counter = 0
    trigger_h = ROOT.TH1D("", "", len(x), x[0], x[(len(x)-1)])
    trigger_times = np.zeros(len(thresholds))
    for i in range(no_events):
        if i % 5000 == 0 and i > 0:
            print "Evaluated {:d} pulse pairs, {:d} had signals".format(i, counter)
            
        signal_clean = digital_filters.butter_lowpass_filter(signal[i,:], fs, cutoff=500e6)
        try:
            peaks = calc.peakFinder(x, signal_clean, thresh=-0.07, min_deltaT=8.)
        except IndexError as e:
            print "Event {0:d}: Index error {1}".format(i, e)
            continue
        except ValueError as e:
            print "Event {0:d}: Index error {1}".format(i, e)
            continue
        if not peaks: # If the fast 'signal' pmt didn't see anything, continue
            continue

        # Calc timestamps for signal events
        idx_error = False
        signal_times = []
        for peak in peaks:
            thresh = signal_clean[peak]*0.4
            try:
                signal_times.append(calc.calcLeadingEdgeTimestamp(x, signal_clean, peak, thresh))
            except IndexError:
                idx_error = True
                break
        if idx_error:
            continue
            
        # Clean trigger signal and find timestamps
        trigger_clean = digital_filters.butter_lowpass_filter(trigger[i,:], fs, cutoff=250e6)
        for j, thresh in enumerate(thresholds):
            try:
                trigger_times[j] = calc.interpolateThreshold(x, trigger_clean, thresh, rise=False)
            except IndexError:
                idx_error = True
                break
        if idx_error:
            continue
        
        # Add entries into histograms
        Q = calcCharge(x, trigger_clean)
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

        # For debugging
        if t < -250:
            plt.plot(x,signal_clean)
            plt.plot(x,trigger[i,:])
            plt.plot(x,trigger_clean)
            #plt.show()

        # Add event trigger characteristics into histograms
        charge_h.Fill(Q*1e12)
        trigger_h.SetContent( np.array( trigger_clean, dtype=np.float64) )
        summed_trigger_h.Add(trigger_h)
        counter = counter + 1 

    print "Evaluated {:d} pulse pairs, {:d} had signals".format(i, counter)

    # Make rootfile to save data to
    outfile = ROOT.TFile(args.outfile, "RECREATE")
    
    # Write histograms to file
    charge_h.Write()
    summed_trigger_h.Scale( 1. / counter )
    summed_trigger_h.Write()
    for cut in sorted(histos.keys()):
        for hist in histos[cut]:
            hist.Write()

    print "Results written to: {0}".format(args.outfile)
