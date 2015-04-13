import Errors
import numpy as np
def Abs2ChainRef(data,channels,chains,t0i=None,t1i=None):
    # This function converts absolute referenced EEG data to chain referenced EEG data
    # data: np.array of multi-channel data
    # channels: names of each channel
    # chains: channel relationships (('Fp1","F7"),("F7","T3"),...)
    from Utility import DataSegment

    if not isinstance(data,np.ndarray):
        data=np.array(data)

    Ny, Nx = data.shape
    tidx = DataSegment(np.arange(0,Nx),t0i,t1i)
    chain_data = np.empty((len(chains),len(tidx)))
    chain_names = []

    for chidx, chpair in enumerate(chains):
        if (chpair[0] not in channels) or (chpair[1] not in channels):
            raise Errors.ParamError('chains','content',
                    'Undefined channel name {}'.format(chpair))
        idx0 = channels.index(chpair[0])
        idx1 = channels.index(chpair[1])
        chain_data[chidx,:] = data[idx0,tidx] - data[idx1,tidx]
        chain_names.append('-'.join(chpair))

    return chain_data, chain_names

def normalize(data,amplitude=1.0):
    dmax=np.max(data)
    data=data/dmax*float(amplitude)
    return data

def filter(x,Wn,order=3,btype=None):
    import scipy.signal as sig
    if btype is None:
        if not isinstance(Wn,list):
            btype='lowpass'
        else:
            btyp='bandpass'
    B,A = sig.iirfilter(order,Wn,btype=btype)

    from Utility import DataSegment

    return sig.lfilter(B,A,x)

def SignalStat(data,r=0.9,Print=True):
    d = np.sort(np.ravel(data))
    idx = int(float(len(d))*(1.0-r)/2.0)
    dmin = d[idx]
    dmax = d[-idx]
    duse = d[np.logical_and(d>=dmin, d<=dmax)]
    dmean = np.mean(duse)
    drms = np.std(duse)
    if Print:
        print "Mean: {}".format(dmean)
        print "RMS:  {}".format(drms)
        print "Max:  {}".format(dmax)
        print "Min:  {}".format(dmin)
    return dmean,drms


def AutoGain(t,data,win_size=100.0,showProg=True,showResult=None):
    def showAll(t,d0,d1,var,si=0):
        import matplotlib.pyplot as plt
        plt.ion()
        plt.subplot(3,1,1)
        plt.plot(t,d0[si,:])
        plt.subplot(3,1,2)
        plt.plot(t,var[si,:])
        plt.subplot(3,1,3)
        plt.plot(t,d1[si,:])

    Ts = (t[-1]-t[0])/(len(t)-1)
    n=int(float(win_size)/Ts)
    hn=n/2
    N=data.shape[-1]
    i_list = np.append(np.arange(0,N,n/10,dtype=int),N-1)
    data_var = np.zeros(data.shape[:-1]+(len(i_list),))
    if showProg:
        proc0 = 0.0
        print "Calculating gain:{:4.1f}%".format(0.0),
    for i,ti in enumerate(i_list):
        if showProg:
            proc1 = float(ti+1)/float(N)*100.0
            if proc1-proc0>0.1:
                print "\rCalculating gain:{:4.1f}%".format(proc1),
                proc0 = proc1
        t0i = ti-hn
        t1i = ti+hn
        if t0i<0:
            t0i=0
        if t1i>=N:
            t1i=N
        data_w = np.take(data,range(t0i,t1i),data.ndim-1)
        var_w = np.sqrt(np.var(data_w,-1))
        data_var[:,i]=var_w
    if showProg:
        print

    from scipy.interpolate import interp1d
    f_var = interp1d(t[i_list],data_var,bounds_error=False)
    data_var = f_var(t)
    data_new = data/data_var
    if showResult is not None:
        showAll(t,data,data_new,data_var,showResult)
    return data_new

if __name__=='__main__':
    import FileIO
    epirec = FileIO.getEPIRecord('CHIMIC')
    EEG = FileIO.LoadEPI(epirec)
    data = EEG['data']
    t = EEG['t']
    win_size = 100.0
    data=AutoGain(t,data,win_size,showResult=5)
