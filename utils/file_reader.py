import re
import os
import sys
import h5py
import glob
import time
import struct
import itertools
import time
import numpy as np
import matplotlib.pyplot as plt

attributes = ["c{0}_horiz_offset",
              "c{0}_horiz_scale",
              "c{0}_num_samples",
              "c{0}_samples",
              "c{0}_vert_offset",
              "c{0}_vert_scale"]

class FileReader(object):
    '''A generic file class
    '''
    def __init__(self, fname):
        '''Initialise a file reader object
        '''
        self._x = np.array([])
        self._channel_data = {}
        self._files = {}
        self._extension = ''
        self.set_file_path(fname)
        
    def set_file_path(self, fname):
        '''Set the file name
        '''
        pass
        
    def get_files(self):
        '''Return the pass file name
        '''
        return [self._files[chan] for chan in self._files]

    def get_xy_data(self, nevents=1):
        '''Return the xy data for active channels
        '''
        pass
    

class Hdf5FileReader(FileReader):
    '''A class to read .h5 files produced by the fetch.py script
    '''
    def __init__(self, fname):
        '''Initialise an hdf5 file reader'''
        super(Hdf5FileReader, self).__init__(fname)
    
    def set_file_path(self, fname):
        '''Set the file path to be read from
        '''
        self._extension = fname.split("/")[-1].split(".")[-1]
        self._files["all"] = fname

    def get_xy_data(self, nevents=1):
        '''Read in all available .h5 files and return the xy arrays contained within
        '''
        start = time.time()
        for channel in self._files:
            self.read_hdf5_file(self._files[channel], nevents=nevents)
        print "Took {0:.3f}s to read {1:d} events from:".format(time.time()-start, nevents)
        for f in self.get_files():
            print f
        return self._x, self._channel_data

    def read_hdf5_file(self, fname, nevents=1):
        """ Read in a h5py file and return the xy arrays contained within
        """
        if self._extension != "h5":
            raise ValueError("The extension of file {0} is not .h5, will not be able to read")

        f = h5py.File(fname, 'r')
            
        # Cast find all active channels in this dataset
        channel_prefix = list( set( map(int, re.findall(r"\d+", " ".join([str(name) for name in f])))))

        y = {}
        for i, channel in enumerate(channel_prefix):
            
            y_off       = f[attributes[4].format(channel)]
            x_off       = f[attributes[0].format(channel)]
            
            dy          = f[attributes[5].format(channel)][0]
            dx          = f[attributes[1].format(channel)][0]
            
            nsamples    = int(f[attributes[2].format(channel)][0])
            data        = f[attributes[3].format(channel)]
            
            self._channel_data[channel] = data[:nevents,:nsamples]*dy - y_off[0]
        self._x = np.linspace(0, dx*nsamples, num=nsamples)*1e9

        f.close()


class TraceFileReader(FileReader):
    '''A class to read in a directory containing .trace files
    '''
    def __init__(self, fname):
        '''Initialise a Trace file reader
        '''
        super(TraceFileReader, self).__init__(fname)
        self._memory_maps = []

    def __del__(self):
        '''If we've made any tempory (memory map) files, delete them
        '''
        for tmp_file in self._memory_maps:
            print "Deleting tempory storage file {0}".format(tmp_file)
            os.remove(tmp_file)
        
    def set_file_path(self, fname):
        '''Set the file path
        '''
        files = []
        if os.path.isfile(fname):
            self._extension = fname.split("/")[-1].split(".")[-1]
            if self._extension != "traces":
                print self._extension
                raise ValueError("{0} is not a .traces file".format(fname))
            files.append(fname)
        else:
            files = glob.glob("{0}*.traces".format(fname))
            if files == []:
                raise ValueError("Can't find any .traces files matching {0}*.traces".format(fname))
            base = files[0].split("/")[-1].split(".")[-3]
            for f in files:
                if f.split("/")[-1].split(".")[-3] != base:
                    raise ValueError("There are multiple datasets in {0}".format(fname))
                
        for f in files:
            channel = int(f.split(".")[-2][-1])
            self._files[channel] = f
        
    def get_xy_data(self, nevents=1):
        '''Read in all available .traces files and return the xy arrays contained within
        '''
        start = time.time()
        for channel in self._files:
            self.read_channel_trace(channel, self._files[channel], nevents=nevents)
        nevents = len(self._channel_data[self._files.keys()[0]])
        print "Took {0:.3f}s to read {1:d} events from:".format(time.time()-start, nevents)
        for f in self.get_files():
            print f
        return self._x, self._channel_data
        
    def read_channel_trace(self, channel, fname, nevents=1):
        """ Read in a .trace file and return the xy arrays contained within
        """
        f = open(fname, 'rb') # Read in the file as binary
            
        header_pattern = '=IBdddd' # (num_samples, sample_bytes, v_off, v_scale, h_off, h_scale, [samples])
        header_size = struct.calcsize(header_pattern)
        binary_header = f.read(header_size)
        header = struct.unpack(header_pattern, binary_header)
        
        # Assign header values to local variables
        nsamples = header[0]
        sample_bytes = header[1]
        
        y_off = header[2]
        x_off = header[4]
        
        dx = header[5]
        dy = header[3]
        
        # Use header values to define global file pattern
        sample_pattern = '={0:d}c'.format(nsamples)
        file_pattern = '{0}{1:d}c'.format(header_pattern, nsamples)
        sample_size = struct.calcsize(sample_pattern)
        fp_size = struct.calcsize(file_pattern)

        # How many traces in this file?
        f.seek(0,2) # move the cursor to the end of the file
        file_size = f.tell()
        nTraces = file_size / fp_size
        if nevents > nTraces:
            nevents = nTraces
            print "Returning {0:d} traces - the full dataset in this file". format(nevents)
            
        # Make containers for the data sets
        # If we're dealing with a large dataset use a memory map
        if nevents*nsamples*len(self._files) < 5e8:
            y = np.zeros( (nevents, nsamples), dtype=np.float32 )
        else:
            time_stamp = time.localtime()
            tmp_file_name = 'memmapped_ch{0}_{1}m{2}s.dat'.format(channel,time_stamp.tm_min,time_stamp.tm_sec)
            y = np.memmap(tmp_file_name,
                          dtype=np.float32,
                          mode='w+',
                          shape=(nevents, nsamples))
            self._memory_maps.append(tmp_file_name)                          
        
        # Make file reading a loop so binary_data is never too large
        counter = 0
        max_read = int(5e4)
        event_chunks = np.full((nevents / max_read), max_read)
        if nevents % max_read != 0:
            event_chunks = np.append(event_chunks, nevents % max_read)

        # Set file pointer back to the start
        f.seek(0,0)
        # Loop over event chunks
        for chunk in event_chunks:
            # Read the data from file
            binary_data = f.read(fp_size*chunk)
            # Unpack the binary data string
            # Note: This bit is SLOW, limited by python looping itself - replacing the numpy stuff
            # with pass yields a < 10% speedup.
            for section in chunked(binary_data, fp_size):
                try:
                    y[counter, :] = np.fromstring(section[header_size:], count=nsamples, dtype=np.int8)*dy - y_off
                except Exception as e:
                    print "Problem reading trace {0}: {1}".format(counter, e)
                counter = counter + 1

        self._channel_data[channel] = y
        self._x = np.linspace(0, dx*nsamples, num=nsamples)*1e9 
        f.close()
        
def chunked(iterable, n, fillvalue=''):
    args = [iter(iterable)] * n
    return itertools.imap(''.join, itertools.izip_longest(*args, fillvalue=fillvalue))    


def roundrobin(*iterables):
    pending = len(iterables)
    nexts = itertools.cycle(iter(it).next for it in iterables)
    
    while pending:
        try:
            for next in nexts:
                yield next()
        except StopIteration:
            pending -= 1
            nexts = itertools.cycle(itertools.islice(nexts, pending))

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser("%prog file [...]")
    parser.add_argument('fname', type=str,
                        help="File(s) to be read in")
    parser.add_argument('-n', '--no_events', type=int, default=10,
                        help="Number of events to display")
    args = parser.parse_args()
    
    fileReader = FileReader('')
    extension = args.fname.split("/")[-1].split(".")[-1]
    if extension == "h5":
        fileReader = Hdf5FileReader(args.fname)
    else:
        fileReader = TraceFileReader(args.fname)
    x, y = fileReader.get_xy_data(nevents=args.no_events)

    events = 100 if args.no_events > 100 else args.no_events

    figures = []
    fig = plt.figure(figsize=(10,6))
    for i, chan in enumerate(y):
        plt.subplot(len(y.keys()), 1, i+1)
        for j in range(events-1):
            plt.plot(x,y[chan][j,:])
        plt.xlabel('Time (ns)')
        plt.ylabel('Voltage (V)')
        plt.title("Channel {0}".format(chan))
        figures.append(fig)
    plt.grid()
    plt.show()
