import ROOT
import numpy as np
import matplotlib.pyplot as plt
import utils.transient_calculations as calc
import pulse_diagnostics as pd

if __name__ == "__main__":
    import argparse
    import utils.file_reader as file_reader
    
    parser = argparse.ArgumentParser("Script to calculate the ratio of early and late light")
    parser.add_argument('infile', type=str,
                        help="File(s) to be read in")
    parser.add_argument('outfile', type=str,
                        help="Name of output (root) file to be generated")
    parser.add_argument('-n', '--no_events', type=int, default=10000,
                        help="Number of events to display")
    parser.add_argument('-w', '--window', type=int, default=20,
                        help="Time window for early light")
    parser.add_argument('-a', '--areaCut', type=float, default=1,
                        help="Area cut threshold, given as a fraction pf the event you want to keep")
    args = parser.parse_args()
    fname = "{0}_{1}ns_{2}pc.root".format(args.outfile,args.window,int(100*(args.areaCut)))
    outfile = ROOT.TFile("rootFiles/{0}".format(fname), "RECREATE")
    extension = args.infile.split("/")[-1].split(".")[-1]
    if extension == "h5":
        myFileReader = file_reader.Hdf5FileReader(args.infile)
    else:
        myFileReader = file_reader.TraceFileReader(args.infile)

    print "Reading files"
    x, y_dict = myFileReader.get_xy_data(nevents=args.no_events)
       
    for i, key in enumerate(y_dict.keys()):

        ratio_h = ROOT.TH1D("ratio_channel{0}_h".format(key),"Ratio of early to late light in channel {0}, with a {1}ns window".format(key,args.window),100,0,1)
        ratio_h.GetXaxis().SetTitle("Early light / Late Light")
        ratio_h.GetYaxis().SetTitle("Frequency")
        y = y_dict[key]
        print "Calculating threshold"
        threshold = calc.threshold_finder(x,y)
        print "Calculating early light for channel {0}".format(key)
        early_light = calc.calcPartialArea(x,y,no_events=args.no_events,window=args.window,threshold=threshold,areaCut=args.areaCut,early=True)
        print "Calculating late light for channel {0}".format(key)
        late_light = calc.calcPartialArea(x,y,no_events=args.no_events,window=args.window,threshold=threshold,areaCut=args.areaCut,early=False)
        print "Calculating ratios"
        for j in range(len(early_light)):
            ratio = early_light[j]/late_light[j]
            ratio_h.Fill(ratio)
            
        ratio_h.Write()
            
    print "Results written to {0}".format(fname)
    
