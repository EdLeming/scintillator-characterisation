import sys
import ROOT
import argparse


if __name__ == "__main__":
    parser = argparse.ArgumentParser("Compare an identically named plot from two root files")
    parser.add_argument("file1", type=str, help="First root file")
    parser.add_argument("file2", type=str, help="Second root file")
    parser.add_argument("plotName",
                        nargs="*",
                        type=str,
                        help="Name(s) of the plot to be grabbed from each file")
    parser.add_argument("--normalise", dest="normalise",
                        action="store_true",
                        help="Do you want to normalise histogram integrals?")
    parser.add_argument("--subtract", dest="subtract",
                        action="store_true",
                        help="Do you want to subtract the results?")
    args = parser.parse_args()

    ROOT.gStyle.SetOptStat(0)
    
    first_file = ROOT.TFile(args.file1)
    second_file = ROOT.TFile(args.file2)

    ############################
    # Grab the multigraphs
    try:
        if len(args.plotName) == 1:
            first_plot = first_file.Get("{0}".format(args.plotName[0]))
            second_plot = second_file.Get("{0}".format(args.plotName[0]))
        elif len(args.plotName) == 2:
            first_plot = first_file.Get("{0}".format(args.plotName[0]))
            second_plot = second_file.Get("{0}".format(args.plotName[1]))
        else:
            print "Don't pass more than two plot names!"
            sys.exit()
            
        first_plot.Integral()
        second_plot.Integral()
    except Exception as e:
        print "Looks like there's no histogram matching the passed name(s)"
        sys.exit()

    if args.normalise:
        first_plot.Sumw2()
        second_plot.Sumw2()
        
        first_plot.Scale( 1. / first_plot.Integral("width") )
        second_plot.Scale( 1. / second_plot.Integral("width") )

        
    first_plot.SetLineStyle(1)
    first_plot.SetLineWidth(1)
    first_plot.SetMarkerColor( ROOT.kBlue )
    first_plot.SetLineColor( ROOT.kBlue )
    
    second_plot.SetMarkerColor( ROOT.kRed )
    second_plot.SetLineColor( ROOT.kRed )

    second_plot.SetLineStyle(1)
    second_plot.SetLineWidth(1)

    first_max = first_plot.GetMaximum()
    second_max= second_plot.GetMaximum()
    m = first_max
    if m < second_max:
        m = second_max
    
    can = ROOT.TCanvas("c1", "c1") 
    
    first_plot.Draw("HIST E")
    second_plot.Draw("HIST E SAME")
    first_plot.GetYaxis().SetRangeUser(0, m*1.1)
    can.Update()
    
    if args.subtract:
        subtraction_h = ROOT.TH1D(first_plot)
        can_sub = ROOT.TCanvas("c2", "c2")
        subtraction_h.Add(second_plot, -1)
        subtraction_h.Draw("")
        can_sub.Update()
    
    raw_input("Hit enter...")
