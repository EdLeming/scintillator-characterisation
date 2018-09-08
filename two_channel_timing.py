import ROOT
import pulse_diagnostics as pd
import numpy as np

def calculate_leading_edge_timestamp(x, y, thresh, rise=False):
    '''
    Step back from peak to get the most noise-free timestamp
    '''
    if rise:
        m = max(y)
    else:
        m = min(y)
    m_index = np.where(y == m)[0][0]
    y_reverse = y[m_index+1:0:-1]
    try:
        timestamp = pd.interpolate_threshold(x[:m_index+1],
                                             y_reverse,
                                             thresh,
                                             rise=rise)
    except:
        timestamp = -999
    return timestamp

if __name__ == "__main__":
    import sys
    import argparse
    import matplotlib.pyplot as plt
    from FileReader import read_h5
    parser = argparse.ArgumentParser("Script to run diagnostrics on pulse shapes")
    parser.add_argument('infile', type=str,
                        help="File(s) to be read in")
    parser.add_argument('outfile', type=str,
                        help="Name of output (root) file to be generated")
    parser.add_argument('-n', '--no_events', type=int, default=10000,
                        help="Number of events to display")
    args = parser.parse_args()

    # Run root in batch mode
    ROOT.gROOT.SetBatch(True)
    
    # Make rootfile to save data to
    outfile = ROOT.TFile(args.outfile, "RECREATE")

    ##########################
    # Read in data and loop over each save channel
    x,y_dict = read_h5(args.infile, nevents=args.no_events)

    ##########################
    # Check we can use the data in this file
    keys = sorted(y_dict.keys())
    if len(keys) < 2:
        print "Only found {0:d} active channels in: {1}".format(len(y_dict.keys()), args.infile)
        sys.exit(0)
    nPulses_1 = len(y_dict[keys[0]][:,0])
    nPulses_2 = len(y_dict[keys[1]][:,0])
    if nPulses_1 != nPulses_2:
        print "Active channels have different size datasets"
        sys.exit(0)
    # Check the pulse polarity and set rise / fall appropriately
    f1 = pd.positive_check(y_dict[keys[0]])
    f2 = pd.positive_check(y_dict[keys[1]])
    rise_1 = False
    if f1:
        rise_1 = True
    rise_2 = False
    if f2:
        rise_2 = True

    #######################
    # Book some histograms
    cfd_histos = []
    cfd_thresholds = np.arange(0.05, 0.8, 0.05)
    for i in cfd_thresholds:
        cfd_histos.append(ROOT.TH1D("cfd_{:d}".format(int(i*100)),
                                    "Constant fraction discrimiator: {:.0f} %".format(i*100),
                                    100, -15, 15))
        cfd_histos[-1].GetXaxis().SetTitle("Pulse separation [ns]")
        
    led_histos = []
    led_thresholds = np.arange(10e-3, 260e-3, 10e-3)
    for i in led_thresholds:
        led_histos.append(ROOT.TH1D("led_{:d}mV".format(int(i*1e3)),
                                    "Leading edge discrimiator: {:.0f}mV".format(i*1e3),
                                    100, -15, 15))
        led_histos[-1].GetXaxis().SetTitle("Pulse separation [ns]")

    ##########################
    # Loop over each pulse set
    # and apply some thresholds
    #
    for i in range(nPulses_1):
        # Select y data for this set
        y1 = y_dict[keys[0]][i,:]
        y2 = y_dict[keys[1]][i,:]

        for j, thresh in enumerate(cfd_thresholds):
            t1 = calculate_leading_edge_timestamp(x, y1, thresh*y1.max(), rise=rise_1)
            t2 = calculate_leading_edge_timestamp(x, y2, thresh*y2.max(), rise=rise_2)
            cfd_histos[j].Fill(t1-t2)

        for j, thresh in enumerate(led_thresholds):
            try:
                t1 = calculate_leading_edge_timestamp(x, y1, thresh, rise=rise_1)
                t2 = calculate_leading_edge_timestamp(x, y2, thresh, rise=rise_2)
            except:
                continue
            led_histos[j].Fill(t1-t2)

        if i % 1000 == 0 and i > 0:
            print "\rEvaluated {:d} pulse pairs".format(i)
            
    ###########################
    # Make some plots and save
    # to file
    gausFit = ROOT.TF1("myGaus", "[0]*TMath::Gaus(x, [1], [2])", -2, 2)
    expoGausFit = ROOT.TF1("expoGaus",
                           ("[0]*( ((1-[1])/(TMath::Sqrt(2*TMath::Pi())*[3]))"
                            "*TMath::Gaus(x, [2], [3])"
                            "+ ([1]/2*[4])*exp(-TMath::Abs(x - [2])/[4]))"),
                           -15, 15)

    x_cfd, x_err_cfd = np.zeros(len(cfd_histos)), np.zeros(len(cfd_histos))
    mean_cfd, sigma_cfd = np.zeros(len(cfd_histos)), np.zeros(len(cfd_histos))
    mean_err_cfd, sigma_err_cfd = np.zeros(len(cfd_histos)), np.zeros(len(cfd_histos))
    for i, hist in enumerate(cfd_histos):
        expoGausFit.SetParameters(hist.GetEntries(), 0.1, hist.GetMean(), 1, 5)
        hist.Fit(expoGausFit, "QEMR")
        x_cfd[i] = float(hist.GetName().split("_")[-1])
        mean_cfd[i] = expoGausFit.GetParameter(2)
        mean_err_cfd[i] = expoGausFit.GetParError(2)
        sigma_cfd[i] = expoGausFit.GetParameter(3)
        sigma_err_cfd[i] = expoGausFit.GetParError(3)
        hist.Write()
        print "Chi2 / NDF: {:.2f} / {:.2f}".format(expoGausFit.GetChisquare(), expoGausFit.GetNDF())
        
    x_led, x_err_led = np.zeros(len(led_histos)), np.zeros(len(led_histos))
    mean_led, sigma_led = np.zeros(len(led_histos)), np.zeros(len(led_histos))
    mean_err_led, sigma_err_led = np.zeros(len(led_histos)), np.zeros(len(led_histos))
    for i, hist in enumerate(led_histos):
        expoGausFit.SetParameters(hist.GetEntries(), 0.1, hist.GetMean(), 1, 5)
        hist.Fit(expoGausFit, "QEMR")
        x_led[i] = int(hist.GetName().split("_")[-1][:-2])
        mean_led[i] = expoGausFit.GetParameter(2)
        mean_err_led[i] = expoGausFit.GetParError(2)
        sigma_led[i] = expoGausFit.GetParameter(3)
        sigma_err_led[i] = expoGausFit.GetParError(3)
        hist.Write()
        
    #############################    
    # Actually make those fuckers
    cfd_mean_graph = ROOT.TGraphErrors(len(x_cfd), x_cfd, mean_cfd, x_err_cfd, mean_err_cfd)
    cfd_mean_graph.GetXaxis().SetTitle("Threshold [%]")
    cfd_mean_graph.GetYaxis().SetTitle("Pulse separation [ns]")
    cfd_mean_graph.SetTitle("")
    cfd_mean_graph.SetName("cfd_means")
    cfd_mean_graph.Write()

    cfd_sigma_graph = ROOT.TGraphErrors(len(x_cfd), x_cfd, sigma_cfd, x_err_cfd, sigma_err_cfd)
    cfd_sigma_graph.GetXaxis().SetTitle("Threshold [%]")
    cfd_sigma_graph.GetYaxis().SetTitle("Pulse separation [ns]")
    cfd_sigma_graph.SetTitle("")
    cfd_sigma_graph.SetName("cfd_sigmas")
    cfd_sigma_graph.Write()

    led_mean_graph = ROOT.TGraphErrors(len(x_led), x_led, mean_led, x_err_led, mean_err_led)
    led_mean_graph.GetXaxis().SetTitle("Threshold [mV]")
    led_mean_graph.GetYaxis().SetTitle("Pulse separation [ns]")
    led_mean_graph.SetTitle("")
    led_mean_graph.SetName("led_means")
    led_mean_graph.Write()

    led_sigma_graph = ROOT.TGraphErrors(len(x_led), x_led, sigma_led, x_err_led, sigma_err_led)
    led_sigma_graph.GetXaxis().SetTitle("Threshold [mV]")
    led_sigma_graph.GetYaxis().SetTitle("Pulse separation [ns]")
    led_sigma_graph.SetTitle("")
    led_sigma_graph.SetName("led_sigmas")
    led_sigma_graph.Write()
