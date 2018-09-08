import re
import sys
import h5py
import glob
import time
import numpy as np
import matplotlib.pyplot as plt

attributes = ["c{0}_horiz_offset",
              "c{0}_horiz_scale",
              "c{0}_num_samples",
              "c{0}_samples",
              "c{0}_vert_offset",
              "c{0}_vert_scale"]

def read_h5(filename, nevents=1):
    """ Read in a h5py file and return the xy arrays contained within
    """
    try:
        f = h5py.File(filename, 'r')
    except:
        print "File {0} does not seem to exist".format(filename)
    
    # Cast find all active channels in this dataset
    channel_prefix = list( set( map(int, re.findall(r"\d+", " ".join([str(name) for name in f]))) ) )

    y = {}
    for i, chan in enumerate(channel_prefix):
        
        y_off       = f[attributes[4].format(chan)]
        x_off       = f[attributes[0].format(chan)]
        
        dy          = f[attributes[5].format(chan)][0]
        dx          = f[attributes[1].format(chan)][0]

        nsamples    = int(f[attributes[2].format(chan)][0])
        data        = f[attributes[3].format(chan)]

        #x = np.tile(np.linspace(x_off, x_off+dx*nsamples, num=nsamples)*1e9, (nevents, 1))
        x = np.linspace(0, dx*nsamples, num=nsamples)*1e9
        y[chan] = data[:nevents,:nsamples]*dy - y_off[0]

    return x, y
        

def pulse_height_spectra(y, dx=0.1e-9):
    '''
    Use the passed x, y scope traces to form a pulse height spectrum
    '''
    return [np.trapz(event, dx=dx) for event in y]
        

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser("%prog file [...]")
    parser.add_argument('files', type=str,
                        help="File(s) to be read in",
                        nargs="+")
    parser.add_argument('-n', '--no_events', type=int, default=10000,
                        help="Number of events to display")
    args = parser.parse_args()

    bins = np.linspace(-6e-8, -0.001e-8, num=50)
    plt.figure()
    for f in args.files:
        start = time.time()
        x, y = read_h5(f, nevents=args.no_events)
        print "Took {0:.1f}s to read {1}".format(time.time()- start, (len(y.keys()), y[y.keys()[0]].shape))

        plt.hist(pulse_height_spectra(y[y.keys()[0]], dx=(x[2]-x[1])*1e-9), bins=bins, label=f.split("/")[-1], histtype="step" )
    plt.xlabel("Pulse integral [Vs]")
    plt.gca().set_yscale("log", nonposy="clip")
    #plt.gca().set_xscale("log", nonposy="clip")
    plt.legend(loc='upper left')
    plt.show()
