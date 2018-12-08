import ROOT
import itertools
import numpy as np
import scipy.signal as signal
import matplotlib.pyplot as plt

import pulse_diagnostics as pd
import digital_filters as filt

def peakFinder(x, y, thresh=-0.07, positive=False, min_deltaT=8., plot=False):
    '''
    A peak finding algorithm
    '''
    # Find all points above / below threashold
    if positive:
        above_thresh = np.where( y > thresh, 1., 0. )
    else:
        above_thresh = np.where( y < thresh, 1., 0. )
    # Differentiate to find first point of each pulse to cross thresh
    diff = np.diff(above_thresh)
    start_indicies = np.where( diff > 0.5 )[0]
    if not start_indicies.any() or (min(y) < -1):
        return []
    stop_indicies = np.where( diff < -0.5)[0]
    peak_indicies = np.zeros( len(start_indicies) )
    # Find the peak at each threshold crossing
    for i in range( len(start_indicies) ):
        if positive:
            peak = max(y[start_indicies[i]:stop_indicies[i]])
        else:
            peak = min(y[start_indicies[i]:stop_indicies[i]])
        peak_indicies[i] = (np.where(y[start_indicies[i]:stop_indicies[i]] == peak )[0][0]
                            + start_indicies[i])
    # Loop through indicies and remove any within min_deltaT of each other
    dx = x[1] - x[0]
    peak_check = True
    while peak_check:
        too_close = np.where( np.diff(peak_indicies)*dx < min_deltaT )[0]
        if too_close.size > 0:
            peak_indicies = np.delete(peak_indicies, too_close[0]+1)
        else:
            peak_check = False

    # Some plotting stuff - probably delete after debugging
    if plot and peak_indicies.size > 1:
        x_select = [x[int(i)] for i in peak_indicies]
        y_select = [y[int(i)] for i in peak_indicies]
        plt.plot(x, y)
        plt.plot(x_select, y_select, 'x')
        plt.axhline(thresh, x[0], x[-1], linestyle="--", linewidth=0.5, alpha=0.5)
        plt.show()
    return peak_indicies

if __name__ == "__main__":
    import argparse
    import time
    import sys
    import FileReader as fr
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
    
    # Read in data and loopg over each save channel
    myFileReader = fr.FileReader('')
    extension = args.infile.split("/")[-1].split(".")[-1]
    if extension == "h5":
        myFileReader = fr.Hdf5FileReader(args.infile)
    else:
        myFileReader = fr.TraceFileReader(args.infile)
    x, y_dict = myFileReader.get_xy_data(nevents=args.no_events)

    y=y_dict[args.channel]
    fs = 1./ ((x[1]-x[0])*1e-9)

    # ROOT stuff
    ROOT.gROOT.SetBatch(True)
    outfile = ROOT.TFile(args.outfile, "RECREATE")
    peaks_h = ROOT.TH1D("EventTimes", "", len(x)-1, x)
    peaks_h.GetXaxis().SetTitle("Time [ns]")

    # Do looping
    start = time.time()
    for i in range(args.no_events):
        clean = filt.butter_lowpass_filter(y[i,:], fs, cutoff=0.5e9)
        try:
            #peaks = peakFinder(x, y[i,:])
            peaks = peakFinder(x, clean, plot=args.plot)
        except Exception as e:
            print e
            continue
        for peak in peaks:
            peaks_h.Fill(x[int(peak)])
        if i % 10000 == 0 and i is not 0:
            print "{0} events processed".format(i) 
    end = time.time()
    print "{0} events took {1:.2f}s to process [{2:.4f} s/evt]".format(i+1,
                                                                       (end-start),
                                                                       (end-start)/(i+1))
    peaks_h.Write()
