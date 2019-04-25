import sys
import ROOT
import matplotlib.pyplot as plt
import numpy as np

import utils.root_plotting as rootplot

class MinuitMinimize( ROOT.TPyMultiGenFunction ):
    ''' A class to hold a minuit minimiser
    '''
    def __init__(self, pdf, time_response, lead_time=20):
        '''
        '''
        # Some required definitions
        self._nDim = 5
        ROOT.TPyMultiGenFunction.__init__(self, self)
        # Hold the current parameters
        self._pars = []
        # Get data histogram's contents and store as private variables for fitting
        x,y = get_xy(pdf)
        errors, contents = [], []
        for i, ent in enumerate(y):
            errors.append(np.sqrt(ent))
        self._x = np.array(x, dtype=np.float128)
        self._dx = round( (x[1] - x[0])*100)*0.01 # round to 10ps resolution
        self._bin_contents = np.array(y, dtype=np.float128)
        self._bin_contents_errors = np.array(errors, dtype=np.float128)
        # How much 'lead time' should we include before the fitted t0?
        self._lead_time = lead_time #ns 
        self._lead_time_offset = int( self._lead_time / self._dx )
        # Some vairable for the chi2 calculations
        self._nActive_bins=0
        self._chi2 = 0
        # Get and store system timing resolution histogram a private array
        det_res_x, det_res_y = get_coincidence_resolution(time_response)
        if np.abs(self._dx - (det_res_x[1] - det_res_x[0])) > 1e-2:
            print "PDF and resolution histos have different binning: {0:e}, {1:e}".format(self._dx, det_res_x[1] - det_res_x[0])
            sys.exit()
        self._resolution = det_res_y
        
    def NDim(self):
        '''REQUIRED FOR MINIMIZER TO WORK
        How many dimensions in this fit
        '''
        return self._nDim

    def DoEval(self, pars):
        '''REQUIRED FOR MINIMIZER TO WORK
        Function called by ROOT minimizer
        '''
        _, fit = self.FitFunc(pars)
        self._chi2 = self.Chi2(fit)
        self._chi2_NDF = self._chi2 / (self._nActive_bins - self._nDim)
        return self._chi2

    def FitFunc(self, pars):
        '''Make functional form given the passed parameters
        '''
        # Save current parameters for use later
        self._pars = pars
        # Define new timebase
        dx = self._dx
        shift_x = self._x - pars[0]
        zero_index = np.where(shift_x > 0)[0]-1
        raw_x = shift_x[zero_index]
        # Define optical model
        scint = (1-pars[4])*self.scintillator_response(raw_x, pars[2], pars[3])
        ceren = (pars[4])*self.delta_response(raw_x)
        total = scint + ceren
        # Define variables to select appropriate part of convolved array
        diff = len(self._resolution) - len(raw_x)
        half_index = int((len(raw_x)+diff)/2.)
        # Do convolution
        total_response = np.convolve(total, self._resolution)[half_index-self._lead_time_offset:-half_index]
        # Normalise to one and scale response to fit NEntries in data histogram
        total_response = pars[1]*(total_response / np.trapz(total_response, dx=1))
        total_x = np.arange(-self._lead_time, (len(total_response)*dx)-self._lead_time, dx)
        return total_x, total_response

    def Chi2(self, fit):
        '''Calculate the Pearsons Chi2
        '''
        chi2, self._chi2, self._nActive_bins  = 0, 0, 0
        offset = len(self._x) - len(fit) + self._lead_time_offset    
        data = self._bin_contents[offset:]
        for i, entry in enumerate(data):
            if entry < 1:
                continue
            diff = (entry - fit[i])
            chi2 = chi2 + (diff*diff / fit[i])
            self._nActive_bins = self._nActive_bins + 1
        return chi2

    def scintillator_response(self, x, rise, fall, normalise=False):
        f = (np.exp(-x/fall) - np.exp(-x/rise)) / (fall - rise)
        if normalise:
            f = f*self._dx
        return f

    def delta_response(self, x, normalise=False):
        f = np.zeros(len(x))
        f[ 0 ] = 1 * (1/self._dx)
        if normalise:
            f = f*self._dx
        return f

    def gaus(self, x, mu=0, sigma=10, normalise=False):
        f =  (np.exp( -np.power(x-mu, 2) / (2*sigma*sigma) ) / np.sqrt(2*np.pi*sigma*sigma))
        if normalise:
            f = f / np.trapz(f,x)
        return f

    
def get_xy(histo):
    '''
    Turn a histogram into x,y arrays
    '''
    x, y = [],[]
    for i in range(histo.GetNbinsX()):
        x.append(histo.GetBinLowEdge(i+1) + histo.GetBinWidth(i+1))
        y.append(histo.GetBinContent(i+1))
    return x,y

def xy_to_hist(x_arr,y_arr):
    '''
    Turn a histogram into x,y arrays
    '''
    return histo

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

def get_coincidence_resolution(histo):
    '''Get histogram of coincidence timing resolution
    '''
    # Get xy values
    x,y = get_xy(histo)
    dx = x[1] - x[0]
    peak_bin = response_h.GetBinLowEdge(response_h.GetMaximumBin())
    y_trimmed = np.trim_zeros(np.array(y))
    # Find the max bin
    y_max_bin = np.argmax(y_trimmed)
    # Make the y array symetric around this bin
    padding = (0,0)
    left = y_max_bin
    right = len(y_trimmed) - left
    if left > right:
        padding = (0, left-right)
    else:
        padding = (right-left, 0)
    # Make symetric arrays
    y_symmetric = np.pad(y_trimmed, padding,'constant',constant_values=(0,0))
    y_symmetric = (y_symmetric / np.trapz(y_symmetric, dx=dx))
    n = len(y_symmetric)
    x_symmetric = np.arange(-(dx*n)/2., (dx*n)/2., dx)
    return x_symmetric, y_symmetric    


def plot_correlation_matrix(minimizer):

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
    import time
    import sys
    import array
    usage = "usage: %prog <filename/prefix>"
    parser = argparse.ArgumentParser(usage)
    parser.add_argument('infile', type=str,
                        help="File(s) to be read in")
    parser.add_argument('-n', '--hist_name', type=str, default="Charge50_t5",
                        help="Name of the histogram to grab from infile [Charge200-300_t-0.020]")
    parser.add_argument('-r', '--resolution_file', type=str,
                        default="./results/system_time_response/tests/cerenkov_thresh-NEMO-reflec.root",
                        help="Path to rootfile containing")
    args = parser.parse_args()

    # ROOT stuff
    ROOT.gROOT.SetBatch(False)
    ROOT.gStyle.SetOptStat(0)
    
    # Get data file
    infile = ROOT.TFile(args.infile)
    histos = []
    try:
        pdf_h = infile.Get("{0}".format(args.hist_name))
        pdf_h.SetDirectory(0)
        pdf_h.Integral()
    except Exception as e:
        raise e
    histos.append(pdf_h)
    
    x,y = get_xy(pdf_h)
    dx = x[1] - x[0]

    # Get system response
    thresh = args.hist_name.split("t")[-1]
    infile = ROOT.TFile(args.resolution_file)
    try:
        response_h = infile.Get("dt_thresh{0}".format(thresh))
        response_h.SetDirectory(0)
        response_h.Integral()
    except Exception as e:
        raise e
    histos.append(response_h)
    
    # Set up minimiser object with custom class object where we define our function
    mini = ROOT.Math.Factory.CreateMinimizer("Minuit2", "Migrad")
    mini.SetMaxFunctionCalls(10000000)
    mini.SetMaxIterations(10000)
    mini.SetTolerance(1)
    mini.SetPrintLevel(2)
    mini.SetStrategy(2)
    mini.SetValidError(True)
    myMinimizer = MinuitMinimize(pdf_h, response_h, lead_time=20)
    mini.SetFunction(myMinimizer)

    # Get some useful values from the hist for seeding fit
    #N = pdf_h.Integral('width')
    N = pdf_h.Integral()
    peak_centre = pdf_h.GetXaxis().GetBinLowEdge( pdf_h.GetMaximumBin() )
    trace_start = x[pdf_h.FindFirstBinAbove( 5 )]
    trace_end = pdf_h.GetXaxis().GetBinLowEdge( pdf_h.GetXaxis().GetNbins() )
    fit_start = trace_start - myMinimizer._lead_time

    pars = [fit_start+1., # Account for offset in fit function
            N,
            2.,
            40.,
            0.01]

    # Set initial values and constraints for fit
    mini.SetVariable(0,"t_0", pars[0], dx/10.)
    mini.SetVariable(1,"N", pars[1], 1)
    mini.SetVariable(2,"Rise", pars[2], 0.001)
    mini.SetVariable(3,"Fall", pars[3], 0.001)
    mini.SetVariable(4,"R_cs", pars[4], 0.00001)
    mini.SetVariableLimits(0, fit_start+0.2, fit_start+1.2)
    mini.SetVariableLimits(1, N*0.95, N*1.05)
    mini.SetVariableLimits(2, 1.25, 2.5)
    mini.SetVariableLimits(3, 35., 50.)
    mini.SetVariableLimits(4, 0.1, 0.001)

    start = time.time()
    mini.Minimize()
    end = time.time()
    mini.Hesse()
    
    # Get final parameters and their errors
    pars, errors = getPythonPars(mini)

    # Print final results to screen
    print "\nResults:"
    print "Chi2/NDF: \t{0:.0f} / {1:.0f}".format(myMinimizer._chi2, (myMinimizer._nActive_bins -
                                                                    myMinimizer._nDim))
    for i, par in enumerate(pars):
        print "{0}: \t{1:.2E} +/- {2:.2E}".format(mini.VariableName(i), par, errors[i])
    print ""    
    print "t0 = initial guess + {0:.3f} ns".format(pars[0] - fit_start)
    print "Fitting took {0:.1f}s".format(end-start)
    print ""

    # Get fit and data - offset appropriately for plotting
    fit_x, fit = myMinimizer.FitFunc( pars )
    offset = len(myMinimizer._x) - len(fit_x) + myMinimizer._lead_time_offset
    data = myMinimizer._bin_contents[offset:]
    plot_x = fit_x[:len(data)]

    # Make data and fit histograms for plotting    
    data_h= ROOT.TH1D("Data_h","",len(plot_x)-1, np.array(plot_x, dtype=np.float64))
    data_h.SetContent( np.array(data, dtype=np.float64) )    
    data_h.GetXaxis().SetTitle("Time residuals [ns]")
    data_h.GetYaxis().SetTitle("Counts / {:.2f} ns".format(dx))
    data_h.SetLineColor( ROOT.kBlack )
    data_h.SetMarkerColor( ROOT.kBlack )

    fit_h= ROOT.TH1D("Fit_h","",len(plot_x)-1, np.array(plot_x, dtype=np.float64))
    for i in range(len(plot_x)-1):
        fit_h.SetBinError(i, 0)
    fit_h.SetContent( np.array(fit, dtype=np.float64) )
    fit_h.SetLineColor( ROOT.kRed )
    fit_h.SetMarkerColor( ROOT.kRed )

    # Make diff and Chi2 histos for inlay plot
    diff_h = rootplot.makeDiffHisto(data_h, fit_h)
    diff_h.SetLineColor( ROOT.kBlack )
    diff_h.SetMarkerColor( ROOT.kBlack )
    diff_h.GetXaxis().SetTitle("Time residuals [ns]")
    diff_h.GetYaxis().SetTitle("Data - fit")
    diff_h.GetXaxis().SetTitleOffset(2.)
    chi2_per_bin_h = rootplot.makeChi2Histo(data_h, fit_h)
    chi2_per_bin_h.SetLineColor( ROOT.kBlack )
    chi2_per_bin_h.SetMarkerStyle( 7 )
    chi2_per_bin_h.SetMarkerSize( 0.5 )
    chi2_per_bin_h.SetMarkerColor( ROOT.kBlack )
    chi2_per_bin_h.GetXaxis().SetTitle("Time residuals [ns]")
    chi2_per_bin_h.GetYaxis().SetTitle("Chi2 per bin")
    chi2_per_bin_h.GetXaxis().SetTitleOffset(2.)
    
    # Only show the interesting range of the histograms we've just built
    last_bin = data_h.FindLastBinAbove( 0 )
    data_h.GetXaxis().SetRange(0, last_bin)
    fit_h.GetXaxis().SetRange(0, last_bin)
    diff_h.GetXaxis().SetRange(0, last_bin)
    chi2_per_bin_h.GetXaxis().SetRange(0, last_bin)
    
    # Make comparison plot
    #can, p1, p2 = rootplot.makeResidualComparison(data_h, fit_h, diff_h)
    can, p1, p2 = rootplot.makeResidualComparison(data_h, fit_h, chi2_per_bin_h)
    
    # Draw a box summarising the results in the main plot
    p1.cd()
    tPave = ROOT.TPaveText(0.63, 0.6, 0.88, 0.96, "NDC")
    tPave.SetFillColor(ROOT.kWhite)
    tPave.SetBorderSize(1)
    tPave.SetTextAlign(12)
    tPave.AddText("Chi2/NDF = {0:.0f} / {1:.0f}".format(myMinimizer._chi2, (myMinimizer._nActive_bins -
                                                                            myMinimizer._nDim)))
    tPave.AddLine(.0,.85,1.,.85);
    tPave.AddText("t_0       = {0:.1f} +/- {1:.2f} ns".format(pars[0], errors[0]))
    tPave.AddText("N_e      = {0:d} +/- {1:d}".format(int(pars[1]), int(errors[1])))
    tPave.AddText("#tau_r        = {0:.2f} +/- {1:.2f} ns".format(pars[2], errors[2]))
    tPave.AddText("#tau_f        = {0:.1f} +/- {1:.1f} ns".format(pars[3], errors[3]))
    tPave.AddText("Ceren/Scint    = {0:.4f} +/- {1:.4f}".format(pars[4], errors[4]))
    tPave.Draw()

    # Draw a line on the inlay plot to guide the eye
    p2.cd()
    last_value = data_h.GetBinLowEdge( last_bin )
    line = ROOT.TLine(-10,0,last_value,0)
    line.SetLineColorAlpha(ROOT.kRed, 0.9)
    line.SetLineStyle(2)
    line.Draw()
    
    # Some extra plots
    can_chi2 = ROOT.TCanvas("Chi2", "Chi2")
    chi2_h = ROOT.TH1D("Chi2","Chi2 per bin distribution: Best fit",100, 0, 10)
    for i, entry in enumerate(data):
        diff = entry - fit[i]
        chi2_h.Fill((diff*diff)/fit[i])
    chi2_h.GetXaxis().SetTitle("Chi2")
    chi2_h.Draw("")

    can_corr = ROOT.TCanvas("Correlations","Correlations")
    corr_matrix_h = plot_correlation_matrix(mini)
    corr_matrix_h.Draw("TEXT")

    can_chi2.Update()
    can_corr.Update()
    can.Update()
    
    raw_input("Hit enter...")
