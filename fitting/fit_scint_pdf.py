import ROOT
import matplotlib.pyplot as plt
import numpy as np

class MinuitMinimize( ROOT.TPyMultiGenFunction ):
    ''' A class to hold a minuit minimiser
    '''
    def __init__(self, hist):
        '''
        '''
        x,y = get_xy(hist)
        errors, contents = [], []
        for i, ent in enumerate(y):
            errors.append(np.sqrt(ent))
        self._x = np.array(x, dtype=np.float128)
        self._dx = x[1] - x[0]
        self._bin_contents = np.array(y, dtype=np.float128)
        self._bin_contents_errors = np.array(errors, dtype=np.float128)
        self._pars = []

        # How much 'lead time' should we include before the fitted t0?
        self._lead_time = 20 #ns 
        self._lead_time_offset = int( self._lead_time / self._dx )
        
        self._nActive_bins=0
        self._chi2 = 0
        self._nDim = 6
        
        ROOT.TPyMultiGenFunction.__init__(self, self)

    def NDim(self):
        return self._nDim

    def DoEval(self, pars):
        ''' Function to be minimized
        '''
        x, func = self.FitFunc(pars)
        self._chi2 = self.Chi2(x, func)
        self._chi2_NDF = self._chi2 / (self._nActive_bins - self._nDim)
        return self._chi2

    def FitFunc(self, pars):
        '''
        '''
        # Save parameters for use later
        self._pars = pars

        # Define timebases
        dx = self._dx

        shift_x = self._x - pars[0]
        t = shift_x[np.where(shift_x > 0)[0]]

        len_x = self._dx*len(t)
        x = np.arange(-len_x/2., len_x/2., dx)

        # Make characteristic shapes
        gaussian = self.gaus(x,0,pars[1])
        scint = self.scintillator_response(t, pars[3], pars[4])
        ceren = self.ceren_response(x, order=0.01)

        # For selecting appropriate regions of the convolution
        half_index = int(len(t)/2.)
        quarter_index = int(len(t)/4.)
        
        # Convolve optical response with detector response
        scint_response = np.convolve(scint, gaussian)[half_index-self._lead_time_offset:-half_index]
        ceren_response = np.convolve(ceren, gaussian)[len(t) - self._lead_time_offset:-quarter_index]
        # Normalize
        scint_response = scint_response / np.trapz(scint_response, dx=dx)
        ceren_response = ceren_response / np.trapz(ceren_response, dx=dx)
        # Scale by estimated ceren-scint ratio
        scint_response = (1-pars[5])*scint_response
        ceren_response = pars[5]*ceren_response
        
        # Zero pad the signals
        zero_padding = np.zeros( len(scint_response) - len(ceren_response) )
        ceren_response = np.append(ceren_response, zero_padding)
        pad_x = np.arange(0, dx*len(ceren_response), dx) - self._lead_time
        
        # Sum the two signal arrays and scale by predicted intensity to form the total response
        total_response = np.array(scint_response)
        for i, ent in enumerate(ceren_response):
            total_response[i] = (total_response[i] + ent)*pars[2]

        #plt.plot(pad_x, scint_response)
        #plt.plot(pad_x, ceren_response)
        #plt.plot(pad_x, total_response,'--',label='total')
        #plt.legend()
        #plt.show()
        return pad_x, total_response

    def scintillator_response(self, x, rise, fall, normalise=False):
        f = (np.exp(-x/fall) - np.exp(-x/rise)) / (fall - rise)
        if normalise:
            f = f / np.trapz(f,x)
        return f
    
    def ceren_response(self, x, order=0.1, normalise=False):
        f = np.exp(- ((x)*(x) )/ order*order ) / (order*np.sqrt(np.pi))
        if normalise:
            f = f / np.trapz(f,x)
        return f

    def gaus(self, x, mu=0, sigma=10, normalise=False):
        f =  (np.exp( -np.power(x-mu, 2) / (2*sigma*sigma) ) / np.sqrt(2*np.pi*sigma*sigma))
        if normalise:
            f = f / np.trapz(f,x)
        return f

    def Chi2(self, x, array):
        '''Calculate the chi2
        '''
        chi2 = 0
        self._nActive_bins = 0
        general_offset = len(self._x) - len(x) + self._lead_time_offset
        no_steps = len(x) - self._lead_time_offset
        diff_a = []
        for i in range(no_steps):
            raw_index = i + general_offset
            if self._bin_contents[raw_index] < 1.0:
                diff_a.append(0)
                continue
            diff = (self._bin_contents[raw_index] - array[i]) / self._bin_contents_errors[raw_index]
            diff_a.append(diff)
            chi2 = chi2 + diff*diff
            self._nActive_bins = self._nActive_bins + 1

        #plt.plot(x, array)
        #plt.plot(x[:-self._lead_time_offset]+self._lead_time, self._bin_contents[general_offset:])
        #plt.plot(x[:-self._lead_time_offset], diff_a)
        #plt.plot(self._bin_contents[general_offset:])

        #plt.plot(array[:-self._lead_time_offset])
        #plt.plot(self._bin_contents[general_offset:])

        #plt.plot(array)
        #plt.plot(self._bin_contents[general_offset:])
        #plt.plot(diff_a)
        #plt.show()
        return chi2
    
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

if __name__ == "__main__":
    import argparse
    import time
    import sys
    import array
    usage = "usage: %prog <filename/prefix>"
    parser = argparse.ArgumentParser(usage)
    parser.add_argument('infile', type=str,
                        help="File(s) to be read in")
    parser.add_argument('-n', '--hist_name', type=str, default="Charge200-300_t-0.020",
                        help="Name of the histogram to grab from infile [Charge200-300_t-0.020]")
    args = parser.parse_args()

    # ROOT stuff
    ROOT.gROOT.SetBatch(False)
    ROOT.gStyle.SetOptStat(0)
    infile = ROOT.TFile(args.infile)

    try:
        time_h = infile.Get("{0}".format(args.hist_name))
    except Exception as e:
        raise e

    x,y = get_xy(time_h)
    dx = x[1] - x[0]

    # Get some useful values from the hist we've just read
    N = time_h.Integral('width')
    peak_centre = time_h.GetXaxis().GetBinLowEdge( time_h.GetMaximumBin() )
    trace_start = x[time_h.FindFirstBinAbove( 3 )]
    trace_end = time_h.GetXaxis().GetBinLowEdge( time_h.GetXaxis().GetNbins() )

    pars = [trace_start-10,
            1.,
            N,
            20.,
            45.,
            0.001]
    
    mini = ROOT.Math.Factory.CreateMinimizer("Minuit2", "Migrad")
    #mini = ROOT.Math.Factory.CreateMinimizer("Minuit2", "Fumili")
    #mini = ROOT.Math.Factory.CreateMinimizer("GSLSimAn", "")
    #mini = ROOT.Math.Factory.CreateMinimizer("GSLMultiMin", "BFGS")
    mini.SetMaxFunctionCalls(1000000)
    mini.SetMaxIterations(100000)
    mini.SetTolerance(0.00000001)
    mini.SetPrintLevel(1)
    mini.SetStrategy(2)
    minuitMini = MinuitMinimize(time_h)
    mini.SetFunction(minuitMini)
    mini.SetVariable(0,"trigger_offset", pars[0], dx)
    mini.SetVariable(1,"det_resolution", pars[1], 0.001)
    mini.SetVariable(2,"N", pars[2], 1)
    mini.SetVariable(3,"Rise", pars[3], 0.00001)
    mini.SetVariable(4,"Fall", pars[4], 0.00001)
    mini.SetVariable(5,"R_cs", pars[5], 0.0000000001)
    mini.SetVariableLimits(0, trace_start-20, trace_start-5)
    mini.SetVariableLimits(1, 0.5, 3.5)
    mini.SetVariableLimits(2, N*0.8, N*1.2)
    mini.SetVariableLimits(3, 12., 25.)
    mini.SetVariableLimits(4, 35., 65.)
    mini.SetVariableLimits(5, 0.06, 0.001)
    start = time.time()
    mini.Minimize()
    end = time.time()
    
    pars, errors = getPythonPars(mini)

    for i, par in enumerate(pars):
        print "{0}: \t{1:.2E} +/- {2:.2E}".format(mini.VariableName(i), par, errors[i])
    print trace_start
    fit_x, fitted = minuitMini.FitFunc( pars )
    general_offset = len(minuitMini._x) - len(fit_x) + minuitMini._lead_time_offset

    data = minuitMini._bin_contents[general_offset:]
    fitted = fitted[:-minuitMini._lead_time_offset]
    plot_x = np.arange(-minuitMini._lead_time, (len(fitted)-minuitMini._lead_time_offset-1)*dx, dx)
    
    fitted_h= ROOT.TH1D("Fit_h","",len(plot_x)-1, np.array(plot_x, dtype=np.float64))
    data_h= ROOT.TH1D("Data_h","",len(plot_x)-1, np.array(plot_x, dtype=np.float64))
    fitted_h.SetContent( np.array(fitted, dtype=np.float64) )
    data_h.SetContent( np.array(data, dtype=np.float64) )    

    data_h.GetXaxis().SetTitle("Time residuals [ns]")
    data_h.GetYaxis().SetTitle("Counts / {:.2f} ns".format(dx))
    
    can = ROOT.TCanvas("c1", "c1", 1200, 800)
    data_h.Draw("E")
    fitted_h.SetLineColor(ROOT.kRed)
    fitted_h.Draw("SAME")
    can.Update()

    tPave = ROOT.TPaveText(0.63, 0.6, 0.88, 0.885, "NDC")
    tPave.SetFillColor(ROOT.kWhite)
    tPave.SetBorderSize(1)
    tPave.SetTextAlign(12)
    tPave.AddText("Chi2/NDF = {0:.0f} / {1:.0f}".format(minuitMini._chi2, (minuitMini._nActive_bins -
                                                                          minuitMini._nDim)))
    tPave.AddLine(.0,.85,1.,.85);
    tPave.AddText("t_0     = {0:.1f} +/- {1:.2f} ns".format(pars[0], errors[0]))
    tPave.AddText("Det resolution = {0:.2f} +/- {1:.2f} ns".format(pars[1], errors[1]))
    tPave.AddText("N_e     = {0:.1f} +/- {1:.1f}".format(pars[2], errors[2]))
    tPave.AddText("#tau_r        = {0:.2f} +/- {1:.2f} ns".format(pars[3], errors[3]))
    tPave.AddText("#tau_f        = {0:.1f} +/- {1:.1f} ns".format(pars[4], errors[4]))
    tPave.AddText("Ceren/Scint    = {0:.4f} +/- {1:.4f}".format(pars[5], errors[5]))
    tPave.Draw()
                                                
    
    raw_input("Hit enter...")
