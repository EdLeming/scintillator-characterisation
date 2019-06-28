import sys
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
import matplotlib.pyplot as plt
import numpy as np

import utils.fit_models as fit_models
import utils.root_tools as root_tools

if __name__ == "__main__":
    import argparse
    import time
    import sys
    import array
    usage = "usage: %prog <filename/prefix>"
    parser = argparse.ArgumentParser(usage)
    parser.add_argument('dataFile', type=str,
                        help="Path to a data file to be fit")
    parser.add_argument('responseFile', type=str,
                        help="Path to data file containing a measurement of the system timing resolution")
    parser.add_argument('-t', '--threshold', type=str, default='cf_0p05',
                        help='Leading edge threshold used to calculate tirgger and signal pulse separation')
    parser.add_argument('-q', '--nemo_cut', type=float, default=20,
                        help='Select only events in which the NEMO tube measured a charge greater than this cut [50pC]')
    parser.add_argument('-c', '--trigger_cut', type=float, default=400,
                        help='Select only events in the trigger tube measured a charge greater than this cut [100pC]')
    parser.add_argument('-s', '--signal_cut', type=float, default=40,
                        help="Charge cut to be applied in pC [None]")
    
    args = parser.parse_args()

    # ROOT stuff
    ROOT.gROOT.SetBatch(False)
    ROOT.gStyle.SetOptStat(0)

    # Get data file
    x_array = np.arange(0, 300, 0.12)
    pdf_h =  root_tools.plotFromNtuple(args.dataFile,
                                       x=x_array,
                                       threshold=args.threshold,
                                       trigger_cut=args.trigger_cut,
                                       nemo_cut=args.nemo_cut,
                                       signal_cut=args.signal_cut)
    x,y = root_tools.getXY(pdf_h)
    dx = x[1] - x[0]

    # Get system response
    response_h =  root_tools.plotFromNtuple(args.responseFile,
                                            x=x_array,
                                            threshold=args.threshold,
                                            trigger_cut=args.trigger_cut,
                                            signal_cut=args.signal_cut)
    
    #can=ROOT.TCanvas("c1","c1")
    #pdf_h.Draw("")
    #can.Update()
    # Set up minimiser object with custom class object where we define our function
    mini = ROOT.Math.Factory.CreateMinimizer("Minuit2", "Migrad")
    mini.SetMaxFunctionCalls(10000000)
    mini.SetMaxIterations(10000)
    mini.SetTolerance(1)
    mini.SetPrintLevel(2)
    mini.SetStrategy(2)
    mini.SetValidError(True)
    myMinimizer = fit_models.SingleExponential(pdf_h, response_h, lead_time=5)
    mini.SetFunction(myMinimizer)

    # Get some useful values from the hist for seeding fit
    N = pdf_h.Integral()
    peak_centre = pdf_h.GetXaxis().GetBinLowEdge( pdf_h.GetMaximumBin() )
    trace_start = x[pdf_h.FindFirstBinAbove( 10 )]
    trace_end = pdf_h.GetXaxis().GetBinLowEdge( pdf_h.GetXaxis().GetNbins() )
    fit_start = trace_start - myMinimizer._lead_time

    pars = [N,
            fit_start+.7, # Account for offset in fit function
            1.7,
            46.,
            0.01]

    # Set initial values and constraints for fit
    #mini.SetFixedVariable(0,"N", pars[0])
    mini.SetVariable(0,"N", pars[0], 1.)
    mini.SetVariable(1,"t_0", pars[1], dx) # ns
    mini.SetVariable(2,"Rise", pars[2], 0.001)
    mini.SetVariable(3,"Fall", pars[3], 0.01)
    mini.SetVariable(4,"R_cs", pars[4], 0.0001)
    mini.SetVariableLimits(0, N*0.999, N)
    mini.SetVariableLimits(1, fit_start+0.3, fit_start+2.)
    mini.SetVariableLimits(2, 1., 3.5)
    mini.SetVariableLimits(3, 30., 50.)
    mini.SetVariableLimits(4, 0.001, 0.2)

    start = time.time()
    mini.Minimize()
    end = time.time()
    mini.Hesse()
    
    # Get final parameters and their errors
    pars, errors = root_tools.getPythonPars(mini)

    # Print final results to screen
    print "\nResults:"
    print "Chi2/NDF: \t{0:.0f} / {1:.0f}".format(myMinimizer._chi2, (myMinimizer._nActive_bins -
                                                                    myMinimizer._nDim))
    for i, par in enumerate(pars):
        print "{0}: \t{1:.2E} +/- {2:.2E}".format(mini.VariableName(i), par, errors[i])
    print ""    
    print "t0 = initial guess + {0:.3f} ns".format(pars[1] - fit_start)
    print "Fitting took {0:.1f}s".format(end-start)
    print ""

    # Get fit and data - offset appropriately for plotting
    fit_x, fit = myMinimizer.FitFunc( pars )
    offset = len(myMinimizer._x) - len(fit_x)
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
    diff_h = root_tools.makeDiffHisto(data_h, fit_h)
    diff_h.SetLineColor( ROOT.kBlack )
    diff_h.SetMarkerColor( ROOT.kBlack )
    diff_h.GetXaxis().SetTitle("Time residuals [ns]")
    diff_h.GetYaxis().SetTitle("Data - fit")
    diff_h.GetXaxis().SetTitleOffset(2.)
    chi2_per_bin_h = root_tools.makeChi2Histo(data_h, fit_h)
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
    #can, p1, p2 = root_tools.makeResidualComparison(data_h, fit_h, diff_h)
    can, p1, p2 = root_tools.makeResidualComparison(data_h, fit_h, chi2_per_bin_h)
    
    # Draw a box summarising the results in the main plot
    p1.cd()
    tPave = ROOT.TPaveText(0.63, 0.6, 0.88, 0.96, "NDC")
    tPave.SetFillColor(ROOT.kWhite)
    tPave.SetBorderSize(1)
    tPave.SetTextAlign(12)
    tPave.AddText("Chi2/NDF = {0:.0f} / {1:.0f}".format(myMinimizer._chi2, (myMinimizer._nActive_bins -
                                                                            myMinimizer._nDim)))
    tPave.AddLine(.0,.85,1.,.85);
    tPave.AddText("N_e      = {0:d} +/- {1:d}".format(int(pars[0]), int(errors[0])))
    tPave.AddText("t_0       = {0:.1f} +/- {1:.2f} ns".format(pars[1], errors[1]))
    tPave.AddText("#tau_r        = {0:.2f} +/- {1:.2f} ns".format(pars[2], errors[2]))
    tPave.AddText("#tau_f        = {0:.1f} +/- {1:.1f} ns".format(pars[3], errors[3]))
    tPave.AddText("Ceren/Scint    = {0:.4f} +/- {1:.4f}".format(pars[4], errors[4]))
    tPave.Draw()

    # Draw a line on the inlay plot to guide the eye
    p2.cd()
    last_value = data_h.GetBinLowEdge( last_bin )
    line = ROOT.TLine(-myMinimizer._lead_time,0,last_value,0)
    line.SetLineColorAlpha(ROOT.kRed, 0.9)
    line.SetLineStyle(2)
    #line.Draw()
    
    # Some extra plots
    can_chi2 = ROOT.TCanvas("Chi2", "Chi2")
    chi2_h = ROOT.TH1D("Chi2","Chi2 per bin distribution: Best fit",100, 0, 10)
    for i, entry in enumerate(data):
            diff = entry - fit[i]
            chi2_h.Fill((diff*diff)/fit[i])
    chi2_h.GetXaxis().SetTitle("Chi2")
    chi2_h.Draw("")

    can_corr = ROOT.TCanvas("Correlations","Correlations")
    corr_matrix_h = root_tools.plotCorrelationMatrix(mini)
    corr_matrix_h.Draw("TEXT")

    can_chi2.Update()
    can_corr.Update()
    can.Update()
    
    raw_input("Hit enter...")
