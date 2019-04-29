import sys
import ROOT
import numpy as np

def bufferToArray(buff, N):
    buff.SetSize(N)
    arr = np.array(buff,copy=True)
    return arr
            
def getPythonPars(mini):
    par_buff = mini.X()
    try:
        par_buff.SetSize(mini.NDim())
        pars = np.array(par_buff,copy=True)
    except:
        print "Couldn't find pars"
        pars = np.zeros(mini.NDim())

    err_buff = mini.Errors()
    try:
        err_buff.SetSize(mini.NDim())
        err = np.array(err_buff,copy=True)
    except:
        print "Couldn't find pars"
        err = np.zeros(mini.NDim())

    return pars, err

def getXY(histo):
    '''
    Turn a histogram into x,y arrays
    '''
    x, y = [],[]
    for i in range(histo.GetNbinsX()):
        x.append(histo.GetBinLowEdge(i+1) + histo.GetBinWidth(i+1))
        y.append(histo.GetBinContent(i+1))
    return x,y

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


def plotFromNtuple(fname,
                   x=np.arange(-50,400,0.2),
                   threshold="cf_0p05",
                   trigger_cut=False,
                   nemo_cut=False):
    '''Read ntuple from file and make pulse separation plot
    '''
    # Read in file and check ntuple exists
    try:
        infile = ROOT.TFile(fname, "READ")
        ntuple = infile.Get('ntuple')
        nEvents = ntuple.GetEntriesFast()
        print "Read ntuple with {:d} entries".format(nEvents)
    except Exception as e:
        print "Problem reading ntuple from: {0}\n{1}".format(fname, e)
        sys.exit(0)        
    # See if threshold branch exists
    if not ntuple.GetBranchStatus(threshold):
        print "Can't find branch {0} in ntuple".format(threshold)
        sys.exit(0)
    if trigger_cut and not ntuple.GetBranchStatus("TriggerQ"):
        print "Can't find branch TriggerQ in ntuple"
        sys.exit(0)
    if nemo_cut and not ntuple.GetBranchStatus("NemoQ"):
        print "Can't find branch NemoQ in ntuple"
        sys.exit(0)

    # Make containers to hold event by event results and set branch addresses
    # - These have to be arrays because of root-python interface
    thresh = np.zeros( 1, dtype=np.float32)
    ntuple.SetBranchAddress(threshold, thresh)
    if trigger_cut:
        triggerQ = np.zeros( 1, dtype=np.float32)
        ntuple.SetBranchAddress("TriggerQ", triggerQ)
    if nemo_cut:
        nemoQ = np.zeros( 1, dtype=np.float32)
        ntuple.SetBranchAddress("NemoQ", nemoQ)

    # Put the requested data into a histogram and plot
    hist_dt = ROOT.TH1D("dt", "", len(x)-1, x)
    title = "Thresh: {0}".format(threshold)
    if trigger_cut:
        title = "{0}, TriggerQ > {1:.1f}".format(title, trigger_cut)
    if nemo_cut:
        title = "{0}, NemoQ > {1:.1f}".format(title, nemo_cut)
    hist_dt.SetTitle(title)
    hist_dt.GetXaxis().SetTitle("Pulse separation [ns]")
    for i, entry in enumerate(ntuple):
        if trigger_cut:
            if triggerQ[0] < trigger_cut:
                continue
        if nemo_cut:
            if nemoQ[0] < nemo_cut:
                continue
        hist_dt.Fill(thresh[0])

    hist_dt.DirectoryAutoAdd(0)
    infile.Close()
    return hist_dt

def plotCorrelationMatrix(minimizer):

    pars, errors = getPythonPars(minimizer)
    nPars = len(pars)
    matrix_h = ROOT.TH2D("CorrelationMatrix",
                         "CorrelationMatrix",
                         len(pars), 1, nPars,
                         len(pars), 1, nPars)
    for i in range(nPars):
        matrix_h.GetXaxis().SetBinLabel(i+1, minimizer.VariableName(i))
        matrix_h.GetYaxis().SetBinLabel(i+1, minimizer.VariableName(i))
        for j in range(nPars):
            matrix_h.SetBinContent(i+1, j+1, minimizer.Correlation(i,j))
    return matrix_h


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser("%prog file [...]")
    parser.add_argument('fname', type=str,
                        help="Path to file containing ntuple")
    parser.add_argument('-t', '--threshold', type=str, default="cf_0p05",
                        help="Threshold required for dt calculation")
    parser.add_argument('-q', '--nemo_cut', type=float, default=False,
                        help="Charge cut to be applied in pC [None]")
    parser.add_argument('-c', '--trigger_cut', type=float, default=False,
                        help="Charge cut to be applied in pC [None]")

    args = parser.parse_args()

    can_dt = ROOT.TCanvas("c1", "c1")
    hist_dt =  plotFromNtuple(args.fname,
                             x=np.arange(-50,400,0.2),
                             threshold=args.threshold,
                             trigger_cut=args.trigger_cut,
                             nemo_cut=args.nemo_cut)

    hist_dt.Draw()
    can_dt.Update()

    '''
    # Draw using standard draw utility
    can = ROOT.TCanvas("c1","c1")
    ntuple.Draw(args.threshold, "NemoQ > {0:.2f}".format(args.charge_cut))
    can.Update()
    '''
    raw_input("Hit enter...")
