import ROOT
import numpy as np

def bufferToArray(buff, N):
    buff.SetSize(N)
    arr = np.array(buff,copy=True)
    return arr
            
def getHistoLowBins(histo):
    '''Get low edges of all bins in histo object
    '''
    n = histo.GetNbinsX()
    width = histo.GetBinWidth(1)
    edges, content, error = np.zeros(n), np.zeros(n), np.zeros(n)
    for i in range(n):
        edges[i] = histo.GetBinLowEdge(i+1) + width
        content[i] = histo.GetBinContent(i+1)
        error[i] = histo.GetBinError(i+1)
    return edges, content, error
    
def makeDiffHisto(hist_1, hist_2):
    '''Make the MC - Data / Data histo
    '''
    bins, contents_1, contents_1_err = getHistoLowBins(hist_1)
    b2, contents_2, contents_2_err = getHistoLowBins(hist_2)
    if not len(contents_1) == len(contents_2):
        raise ValueError

    histo = ROOT.TH1D("difference", "", len(contents_1)-1, bins)
    for i in range(0,len(bins)):
        if contents_1[i] == 0 or contents_2[i] == 0:
            new_contents = 0
            new_err = 0
        else:
            new_contents = (contents_1[i] - contents_2[i])
            new_err = ( np.sqrt( (contents_1_err[i]*contents_1_err[i])
                                 + (contents_2_err[i]*contents_2_err[i]) ))
        histo.SetBinContent(i, new_contents)
        histo.SetBinError(i, new_err)
    return histo

def makeChi2Histo(data_h, fit_h):
    '''Make the Chi2 per bin histogram
    '''
    data_bins, data, data_err = getHistoLowBins(data_h)
    fit_bins, fit, fit_err = getHistoLowBins(fit_h)
    if not len(data) == len(fit):
        raise ValueError

    histo = ROOT.TH1D("Chi2_per_bin", "Chi2", len(data)-1, data_bins)
    for i in range(0,len(data_bins)):
        if data[i] == 0 or fit[i] == 0:
            new_contents = 0
        else:
            diff = data[i] - fit[i]
            new_contents = diff*diff / fit[i]
        new_err = 0
        histo.SetBinContent(i, new_contents)
        histo.SetBinError(i, new_err)
    return histo

def makeResidualComparison(hist_1, hist_2, diff_hist, name="c1", line=None):
    '''
    Plot two histos on top of eachother plus provide an
    inlay of the residuals.
    '''
    can = ROOT.TCanvas(name, name, 800, 700)
    can.Draw()
    p1 = ROOT.TPad("p1", "A pad 75% of the height",0.0,0.37,1.0,1.0)
    p2 = ROOT.TPad("p2", "a pad 35% of the height",0.0,0.0,1.0,0.35)
            
    p1.SetTopMargin(0.02)
    p1.SetBottomMargin(0.02)
    p2.SetTopMargin(0.02)
    p2.SetBottomMargin(0.3)

    p1.Draw()
    p2.Draw()
    
    p1.cd()
    hist_1.Draw("E1")
    hist_2.Draw("HIST SAME")
    hist_1.SetTitleOffset(5, 'x')
    hist_1.SetLabelColor(ROOT.kWhite, 'x')
    hist_1.SetLabelOffset(0.015, 'x')
    hist_1.SetLabelSize(0.01, 'x')
    hist_1.SetTickLength(0.015, "xy")
    p1.SetTickx(1)
    p1.SetTicky(1)


    p2.cd()
    diff_hist.Draw("AXIS")
    if line:
        line.Draw()
    diff_hist.DrawCopy("P SAME")
    diff_hist.SetTitleSize(.07, "xy")
    diff_hist.SetLabelSize(.06, "xy")
    diff_hist.SetTitleOffset(1.2, "x")
    diff_hist.SetTitleOffset(0.45, "y")
    diff_hist.SetLabelOffset(0.02, "x")
    diff_hist.SetTickLength(0.015, "xy")
    p2.SetTickx(1)
    p2.SetTicky(1)
    
    can.Update()

    return can, p1, p2


if __name__ == "__main__":
    import sys
    import argparse
    parser = argparse.ArgumentParser("%prog file [...]")
    parser.add_argument('fname', type=str,
                        help="Path to file containing ntuple")
    parser.add_argument('-t', '--threshold', type=str, default="cf_0p05",
                        help="Threshold required for dt calculation")
    parser.add_argument('-q', '--charge_cut', type=float, default=50,
                        help="Charge cut to be applied [50 pC]")
    args = parser.parse_args()

    # Read in file and check ntuple exists
    try:
        infile = ROOT.TFile(args.fname, "READ")
        ntuple = infile.Get('ntuple')
        nEvents = ntuple.GetEntriesFast()
        print "Read ntuple with {:d} entries".format(nEvents)
    except Exception as e:
        print "Problem reading ntuple from: {0}".format(args.fname)
        sys.exit(0)        
    # See if threshold branch exists
    if not ntuple.GetBranchStatus(args.threshold):
        print "Can't find branch {0} in ntuple".format(args.threshold)
        sys.exit(0)

    # Make containers to hold event by event results and set branch addresses
    # - These have to be arrays because of root-python interface
    threshold = np.zeros( 1, dtype=np.float32)
    charge = np.zeros( 1, dtype=np.float32)
    ntuple.SetBranchAddress(args.threshold, threshold)
    #ntuple.SetBranchAddress("NemoQ", charge)
    ntuple.SetBranchAddress("TriggerQ", charge)

    # Put the requested data into a histogram and plot
    can_dt = ROOT.TCanvas("dt","Custom plot")
    dx = np.arange(-50, 400, 0.2)
    hist_dt = ROOT.TH1D("dt_custom",
                        "Thresh: {0}, NemoQ > {1:.1f}".format(args.threshold,args.charge_cut),
                        len(dx)-1, dx)
    hist_dt.GetXaxis().SetTitle("Pulse separation [ns]")
    for i, entry in enumerate(ntuple):
        if charge[0] > args.charge_cut:
            hist_dt.Fill(threshold[0])
    hist_dt.Draw()
    can_dt.Update()

    '''
    # Draw using standard draw utility
    can = ROOT.TCanvas("c1","c1")
    ntuple.Draw(args.threshold, "NemoQ > {0:.2f}".format(args.charge_cut))
    can.Update()
    '''
    raw_input("Hit enter...")
