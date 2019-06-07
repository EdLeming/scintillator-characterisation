import sys
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
import numpy as np
import matplotlib.pyplot as plt

import utils.root_tools as root_tools

class ResponseModel( ROOT.TPyMultiGenFunction ):
    ''' A base class compatible with the miunuit minimier
        It generates a model from a set of passed parameters and returns a teststatisic taken
        against the passed data set.
    '''
    def __init__(self, data_h, time_response_h, lead_time=5, nDim=5):
        ''' Initialise a response model datafile
        '''
        # Some required definitions
        self._nDim = nDim
        ROOT.TPyMultiGenFunction.__init__(self, self)
        # Hold the current parameters
        self._pars = []
        # Get data histogram's contents and store as private variables
        x,y = root_tools.getXY(data_h)
        errors, contents = [], []
        for i, ent in enumerate(y):
            errors.append(np.sqrt(ent))
        self._x = np.array(x, dtype=np.float128)
        self._dx = round( (x[1] - x[0])*100)*0.01 # round to 10ps resolution
        self._bin_contents = np.array(y, dtype=np.float128)
        self._bin_contents_errors = np.array(errors, dtype=np.float128)
        # Get and store system timing resolution histogram a private array
        det_res_x, det_res_y = self.GetResolution(time_response_h)
        if np.abs(self._dx - (det_res_x[1] - det_res_x[0])) > 1e-2:
            print "PDF and resolution histos have different binning: {0:e}, {1:e}".format(self._dx,
                                                                                          det_res_x[1] - det_res_x[0])
            sys.exit()
        self._resolution = det_res_y
        # How much 'lead time' should we include before the fitted t0?
        self._lead_time = lead_time #ns 
        self._lead_time_offset = int( self._lead_time / self._dx )
        # Some vairable for the chi2 calculations
        self._nActive_bins=0
        self._chi2 = 0
        
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
        ''' A place holder method - to be superceeded by methods in daughter classes
        '''
        pass

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
    
    def ScintillatorModel(self, x, rise, fall, normalise=False):
        f = (np.exp(-x/fall) - np.exp(-x/rise)) / (fall - rise)
        if normalise:
            f = f*self._dx
        return f

    def CerenkovModel(self, x, normalise=False):
        f = np.zeros(len(x))
        f[ 0 ] = 1 * (1/self._dx)
        if normalise:
            f = f*self._dx
        return f

    def Gaus(self, x, mu=0, sigma=10, normalise=False):
        f =  (np.exp( -np.power(x-mu, 2) / (2*sigma*sigma) ) / np.sqrt(2*np.pi*sigma*sigma))
        if normalise:
            f = f / np.trapz(f,x)
        return f

    def GetResolution(self, time_resolution_h):
        '''Get histogram of coincidence timing resolution and make it symmetric for
           convoling with optical model
        '''
        # Get xy values
        x,y = root_tools.getXY(time_resolution_h)
        dx = x[1] - x[0]
        peak_bin = time_resolution_h.GetBinLowEdge(time_resolution_h.GetMaximumBin())
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


class  SingleExponential(ResponseModel):
    ''' A class to fit a sigle exponential optical model
    '''
    def __init__(self, data, time_response, lead_time=5):
        '''Initialise a single exponential model
        '''
        super(SingleExponential, self).__init__(data,
                                                time_response,
                                                nDim=5, # This model has five dimensions
                                                lead_time=lead_time)
        
    def FitFunc(self, pars):
        '''Make functional form considering a single exponential optical model considering
           the passed parameters
        '''
        # Save current parameters for use later
        self._pars = pars
        # Define new timebase
        dx = self._dx
        shift =  round(pars[1]*100)*0.01 # round to 50ps
        shift_x = self._x - shift
        zero_index = np.where(shift_x > 0)[0]-1
        raw_x = shift_x[zero_index]
        #print pars[1], shift, raw_x
        # Define optical model
        scint = (1-pars[4])*self.ScintillatorModel(raw_x, pars[2], pars[3])
        ceren = (pars[4])*self.CerenkovModel(raw_x)
        total = scint + ceren
        # Define variables to select appropriate part of convolved array
        diff = len(self._resolution) - len(raw_x)
        half_index = int((len(raw_x)+diff)/2.)
        # Do convolution
        total_response = np.convolve(total, self._resolution)[half_index-self._lead_time_offset:-half_index]
        # Normalise to one and scale response to fit NEntries in data histogram
        total_response = pars[0]*(total_response / np.trapz(total_response, dx=1))
        total_x = np.arange(-self._lead_time, (len(total_response)*dx)-self._lead_time, dx)
        return total_x, total_response
