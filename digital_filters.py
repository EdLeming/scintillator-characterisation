import numpy as np
import scipy.signal as signal
import matplotlib.pyplot as plt

def butter_lowpass(cutoff, fs, order=5):
    nyq = 0.5 * fs
    normal_cutoff = cutoff / nyq
    b, a = signal.butter(order, normal_cutoff, btype='low', analog=False)
    return b, a

def butter_lowpass_filter(data, fs, cutoff=750e6, order=5):
    b, a = butter_lowpass(cutoff, fs, order=order)
    y = signal.filtfilt(b, a, data)
    return y


if __name__ == "__main__":
    import argparse
    import time
    from FileReader import read_h5
    parser = argparse.ArgumentParser("Calculate the time response of a PMT from LED data")
    parser.add_argument('infile', type=str,
                        help="Data file to be read in")
    parser.add_argument('-i', '--event_index', type=int, default=0,
                        help="Index of event to evaluate")
    parser.add_argument("-c", "--channel", type=int, default=2,
                        help="Which scope channel should be evaluated?")
    parser.add_argument("-o", "--order", type=int, default=5,
                        help="The order of the filter to be made [5]")
    parser.add_argument("-t", "--cut_off", type=float, default=500e6,
                        help="Cut off frequency of filter [500 MHz]")
    args = parser.parse_args()

    # Read in data and loop over each save channel
    x,y_dict = read_h5(args.infile, nevents=args.event_index+1)
    
    # Filter requirements.
    order = args.order
    fs = 1. / ((x[1]-x[0])*1e-9) # sample rate, Hz
    cutoff = args.cut_off        # desired cutoff frequency of the filter, Hz

    # Get the filter coefficients so we can check its frequency response.
    b, a = butter_lowpass(cutoff, fs, order)
    
    # Plot the frequency response.
    w, h = signal.freqz(b, a, worN=8000)
    plt.subplot(2, 1, 1)
    plt.plot(0.5*fs*w/np.pi, np.abs(h), 'b')
    plt.plot(cutoff, 0.5*np.sqrt(2), 'ko')
    plt.axvline(cutoff, color='k')
    plt.xlim(0, 0.5*fs)
    plt.title("Lowpass Filter Frequency Response")
    plt.xlabel('Frequency [Hz]')
    plt.grid()
                                
    # Filter the data, and plot both the original and filtered signals.
    data = y_dict[args.channel]
    data = data[args.event_index,:]
    start = time.time()
    y = butter_lowpass_filter(data, cutoff, fs, order=order)
    print "Took {0:.3}s to run filter".format(time.time() - start)
    
    plt.subplot(2, 1, 2)
    plt.plot(x, data, 'b-', label='data')
    plt.plot(x, y, 'g-', linewidth=2, label='filtered data')
    plt.xlabel('Time [sec]')
    plt.grid()
    plt.legend()

    plt.subplots_adjust(hspace=0.35)
    plt.show()
