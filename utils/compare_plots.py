import ROOT
import argparse


if __name__ == "__main__":
    parser = argparse.ArgumentParser("Compare an identically named plot from two root files")
    parser.add_argument("file1", type=str, help="First root file")
    parser.add_argument("file2", type=str, help="Second root file")
    parser.add_argument("plotName", type=str, help="Name of the plot to be grabbed from each file")
    parser.add_argument("--normalise", dest="normalise",
                        action="store_true",
                        help="Do you want to normalise histogram integrals?")
    args = parser.parse_args()

    ROOT.gStyle.SetOptStat(0)
    
    first_file = ROOT.TFile(args.file1)
    second_file = ROOT.TFile(args.file2)
    
    ############################
    # Grab the multigraphs
    try:
        first_plot = first_file.Get("{0}".format(args.plotName))
        second_plot = second_file.Get("{0}".format(args.plotName))

        first_plot.Integral()
        second_plot.Integral()
    except Exception as e:
        raise e

    if args.normalise:
        first_plot.Sumw2()
        second_plot.Sumw2()
        
        first_plot.Scale( 1. / first_plot.Integral() )
        second_plot.Scale( 1. / second_plot.Integral() )
        

    second_plot.SetLineColor( ROOT.kBlue )

    can = ROOT.TCanvas("c1", "c1") 
    
    first_plot.Draw("")
    second_plot.Draw("HIST SAME")
    
    raw_input("Hit enter...")
