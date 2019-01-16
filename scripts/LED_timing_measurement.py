import ROOT

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
    parser.add_argument("-c", "--channel", type=int, default=2,
                        help="Which scope channel should be evaluated?")
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

    y=y_dict[args.channel]
    fs = 1./ ((x[1]-x[0])*1e-9)

    # ROOT stuff
    ROOT.gROOT.SetBatch(True)
    outfile = ROOT.TFile(args.outfile, "RECREATE")
    peaks_h = ROOT.TH1D("EventTimes", "", len(x)-1, x)
    peaks_h.GetXaxis().SetTitle("Time [ns]")

    # Loop over events and calculate timestamps for observed events
    start = time.time()
    for i in range(args.no_events):
        clean = digital_filters.butter_lowpass_filter(y[i,:], fs, cutoff=500e6)
        try:
            peaks = calc.peakFinder(x, clean, thresh=-0.07, min_deltaT=8., plot=args.plot)
        except IndexError as e:
            print "Event {0:d}: Index error {1}".format(i, e)
            continue
        except ValueError as e:
            print "Event {0:d}: Index error {1}".format(i, e)
            continue
        
        for peak in peaks:
            thresh = clean[peak]*0.4
            timestamp = calc.calcLeadingEdgeTimestamp(x, clean, peak, thresh)
            peaks_h.Fill(timestamp)
        if i % 10000 == 0 and i is not 0:
            print "{0} events processed".format(i) 
    end = time.time()
    print "{0} events took {1:.2f}s to process [{2:.4f} s/evt]".format(i+1,
                                                                       (end-start),
                                                                       (end-start)/(i+1))
    peaks_h.Write()
    print "Results written to: {0}".format(args.outfile)
