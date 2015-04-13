import os,sys
import Errors
EEG_chain_bipolar = [
                ("Fp1", "F7"), ("F7", "T3"), ("T3", "T5"), ("T5", "O1"),
                ("Fp1", "F3"), ("F3", "C3"), ("C3", "P3"), ("P3", "O1"),
                ("Fp2", "F8"), ("F8", "T4"), ("T4", "T6"), ("T6", "O2"),
                ("Fp2", "F4"), ("F4", "C4"), ("C4", "P4"), ("P4", "O2"),
                ("Fz" , "Cz"), ("Cz", "Pz")
                ]
def LoadTxtEEG(path,ign_col=[],ampsat=10.0):
    import re
    import numpy as np
    def parseheader(line):
        import time
        keys=[
                'Original file start/end time',
                'Exported file start/end time',
                'Units',
                "Patient's Name",
                'Sampling Rate'
                ]
        for ki,k in enumerate(keys):
            pat=re.compile(r'% *{}:?\t+(?P<key>.*)'.format(k))
            m=pat.match(line)
            if m:
                val=m.group('key').strip()
                if ki==0 or ki==1:
                    val=re.split(r'\t+',val)
                    for vi,v in enumerate(val):
                        val[vi]=time.strptime(v,'%m/%d/%Y %H:%M:%S')
                if ki==4:
                    val=float(val.split(' ')[0])
                return k,val
        return None,None

    def parsefile(filename,ign_col=[],ampsat=10.0):
        f=open(filename,'r')
        for i,l in enumerate(f):
            pass
        numlines=i+1

        header=dict()
        times=[]
        content=[]
        tstr=[]
        f.seek(0)
        print "Parsing {}".format(filename)

        from progressbar import ProgressBar
        pbar=ProgressBar(maxval=numlines).start()
        for lidx,line in enumerate(f):
            if lidx % 100 ==0:
                pbar.update(lidx)
            if line.strip()[0]=='%':
                key,value=parseheader(line)
                if key:
                    header[key]=value
            else:
                issat=False
                if 'AMPSAT' in line:
                    line=line.replace('AMPSAT','inf')
                    issat=True
                line_list=re.split('\t+',line)
                if not tstr or lidx==numlines-1:
                    tstr.append(line_list[0].strip())
                times.append(int(line_list[1]))
                del line_list[:2]
                if line_list[-1].strip().upper() in ["OFF","ON"]:
                    del line_list[-1]
                if len(line_list)>19:
                    for iidx, col in enumerate(ign_col):
                        del line_list[col-iidx]
                if len(line_list)>19:
                    line_list=line_list[:19]
                num_list = np.array(line_list,dtype='|S4').astype(np.float)
                if issat:
                    if content:
                        c1=content[-1]
                        cols=np.isinf(num_list)
                        num_list[cols]=np.sign(c1[cols])*ampsat
                content.append(num_list)
        pbar.finish()
        f.close()
        if 'Sampling Rate' not in header:
            from time import strptime,mktime
            tstr = [ mktime(strptime(ts,"%m/%d/%Y %H:%M:%S")) for ts in tstr ]
            header['Sampling Rate']=float(len(times))/float(tstr[1]-tstr[0])
        return times,content,header

    def cmpFileNames(f1,f2):
        fn1,ext1=os.path.splitext(f1)
        fn2,ext2=os.path.splitext(f2)
        v=cmp(ext1.upper(),ext2.upper())
        if v:
            return v

        fn1=fn1.upper()
        fn2=fn2.upper()
        
        p=re.compile(r'(?P<idx>\d{1,})\Z')
        fnv1=re.sub(p,'',fn1)
        fnv2=re.sub(p,'',fn2)

        v=cmp(fnv1,fnv2)
        if v:
            return v

        m1=p.search(fn1)
        m2=p.search(fn2)
        if m1 and m2:
            return cmp(int(m1.group('idx')),int(m2.group('idx')))
        elif m1:
            return -1
        elif m2:
            return 1
        else:
            return cmp(f1,f2)

    if os.path.exists(path):
        if os.path.isdir(path):
            files=os.listdir(path)
            p=re.compile(r'(?P<idx>\d{1,}).TXT\Z')
            files = [ f for f in files if p.search(f.upper()) ]
            files.sort(cmp=cmpFileNames)
            files = [ os.path.join(path,f) for f in files ]
        elif os.path.isfile(path):
            files=[path]
    times=[]
    contents=[]
    t,c,header = parsefile(files[0],ign_col,ampsat)
    t_init=t[0]
    fs=header['Sampling Rate']
    times.append((np.array(t)-t_init)/fs)
    contents.append(np.array(c).T)
    for filepath in files[1:]:
        t,c,h=parsefile(filepath,ign_col,ampsat)
        times.append((np.array(t)-t_init)/fs)
        contents.append(np.array(c).T)
    return times,contents,fs,header

def showSpectrum(data,fs):
    import matplotlib.pyplot as plt
    import scipy.signal as sigpack
    import numpy as np
    fig=plt.figure()
    freq,spec = sigpack.welch(np.concatenate(data),fs=fs)
    spec=np.log(spec)
    plt.plot(freq,spec)
    plt.show()

def LoadCateEEG(recname, t0=None, t1=None,
        amplitude=None, offset=None, filter_ext = 10.0, segsep=-1,basepath=None):
    import numpy as np
    if not basepath:
        basepath="../SharedData/Hoda"
    chan_file=os.path.join(basepath,'CATE_Channels.txt')
    if not os.path.exists(chan_file):
        raise Exception("Cannot find channel file {}".format(chan_file))

    folder=os.path.join(basepath,recname)
    if not os.path.exists(folder):
        raise Exception("Cannot find record folder {}".format(folder))

    f=open(chan_file,'r')
    chan_names=f.readlines()
    f.close()
    ign_list=[ chidx for chidx,chan in enumerate(chan_names) if chan.strip()=='IGN' ]
    chan_names = [ chan.strip() for chan in chan_names if chan.strip()!='IGN' ]

    times,contents,fs,header=LoadTxtEEG(folder,ign_list)

    import SignalTools

    t_init=times[0][0]

    for idx in range(len(contents)):
        t=times[idx]
        data=contents[idx]

        # elongate the data
        N_ext=int(float(filter_ext)*fs)
        if N_ext>len(t):
            N_ext=len(t)
        d_ext = np.take(data,range(0, N_ext), axis=data.ndim-1)
        data_long = np.concatenate((d_ext,data,d_ext),axis=data.ndim-1)

        ##Convert to mean value reference
        data_long = data_long - np.mean(data_long,axis=0)

        # notch filter 60Hz
        center_freq=60.0
        notch_filter=(np.zeros(2,dtype=float)+center_freq+np.array([-1.0,1.0]))/fs*2.0
        data_long = SignalTools.filter(data_long,notch_filter,3,'bandstop')

        # bandpass filter
        pass_band=[0.5,70]
        bandpass_filter=np.array(pass_band)/fs*2.0
        data_long = SignalTools.filter(data_long,bandpass_filter,5,'bandpass')

        # cut the data
        contents[idx] = np.take(data_long, range(N_ext,data_long.shape[data_long.ndim-1]-N_ext),
                axis=data.ndim-1)

        # fix time
        if idx>0:
            if segsep==0:
                t0=times[idx-1][-1]+1.0/fs
                times[idx]=t-t[0]+t0
            elif segsep>0:
                t0=times[idx-1][-1]+segsep
                times[idx]=t-t[0]+t0
            if segsep!=0:
                contents[idx-1][:,-1]=0
                contents[idx][:,0]=0

    # Join all segments
    t=np.concatenate(times)
    data=np.concatenate(contents,axis=1)

    # offset and amplitude adjusting
    if offset is None:
        data = data - np.atleast_2d(np.mean(data,-1)).T
    if amplitude is None:
        dmean,drms=SignalTools.SignalStat(data)
        data = data / drms * 0.02

    EEG=dict()
    EEG['t']=t
    EEG['data']=data
    EEG['names']=chan_names
    EEG['fs']=fs
    EEG['file']=folder
    #EEG['bandpass_filter']=bandpass_filter
    #EEG['notch_filter']=notch_filter
    return EEG

def LoadEDF(filename, t0=None, t1=None, chains = None,
        filter_fc = [0.02, 0.30], filter_order = 6):
    """Load EDF file from t0 to t1 (in sec)
    if t0 and t1 are floats, they indicates seconds
    chains is a list of duple defines the reference chain
    filter defines the cut-off frequency of a bend-pass filter
    filter_order is the order of this IIR filter

    Output: tref, chan_data, chan_names
    """
    import eegtools.io.edfplus as edfplus
    from bisect import bisect
    import numpy as np
    import SignalTools as SigProc
    from Utility import DataSegment

    data = edfplus.load_edf(filename)
    times = np.array(data.time,dtype=float)
    if t0 is None:
        t0i = None
    else:
        t0i = bisect(times,float(t0))

    if t1 is None:
        t1i = None
    else:
        t1i = bisect(times,float(t1))


    chan_names = data.chan_lab
    strip_keys = [ 'EEG', '-REF', '- REF' ]
    for idx, name in enumerate(chan_names):
        for k in strip_keys:
            up_name = name.upper()
            v = up_name.find(k)
            if v>=0:
                name = name[:v] + name[v+len(k):]
        chan_names[idx]=name.strip()

    if chains is not None:
        chan_data,chan_names = SigProc.Abs2ChainRef(data.X, chan_names, chains, t0i, t1i)
    else:
        chan_data = DataSegment(data.X,t0i,t1i)

    if filter_fc is not None:
        chan_data = SigProc.filter(chan_data,filter_fc,filter_order)

    tref = DataSegment(times,t0i,t1i)

    return tref, chan_data, chan_names

def LoadVoiceData(filepath):
    import os.path
    import numpy as np
    realpath = os.path.join("..","SharedData","an4","wav",filepath)
    data = np.memmap(realpath,dtype='>i2',mode='r')
    VMAX = float(2**15)
    data = np.asarray(data,dtype=float)/VMAX
    t = np.arange(len(data),dtype=float)/16000.0
    return t,data


def getEPIRecord(rec_name,basepath=None):
    if basepath is None:
        basepath='../SharedData/PublicEEG'
    infoxml = os.path.join(basepath,'Descriptions.xml')
    if not os.path.exists(infoxml):
        raise Errors.FileIOFileNotFound(infoxml)
    from xml.etree import ElementTree
    import re
    xmldoc = ElementTree.parse(infoxml)
    record = xmldoc.find(rec_name)
    if record is None:
        raise Errors.FileIORecordNotFound(rec_name)
    fields = record.getchildren()
    rec_dict = dict()
    for fd in fields:
        if fd.tag in (
                'Num_seizures',
                'Age_at_recording',
                'Age_at_surgery',
                'Num_inter-ictal_events'):
            val = int(fd.text)
        elif fd.tag == 'Hardware_filters':
            ft_txt = re.split('-|Hz',fd.text)
            val = [ float(ft_txt[0]), float(ft_txt[1]) ]
        elif fd.tag == 'Software_filters':
            ft_txt = re.split('Hz',fd.text)
            val = float(ft_txt[0])
        else:
            val = fd.text
        rec_dict[fd.tag]=val
    rec_dict['EEG_file']=os.path.join(basepath,rec_name,rec_name+'_EEG_DATA.edf')
    rec_dict['basepath']=basepath
    return rec_dict

def getEPIRecordNames(basepath=None):
    if basepath is None:
        basepath = '../SharedData/PublicEEG'
    infoxml = os.path.join(basepath,'Descriptions.xml')
    if not os.path.exists(infoxml):
        raise Errors.FileIOFileNotFound(infoxml)
    from xml.etree import ElementTree
    recs = ElementTree.parse(infoxml).getroot().getchildren()
    rec_names = [ r.tag for r in recs ]
    return rec_names

def LoadEPI(epirec, t0=None, t1=None,
        amplitude=None, offset=None, filter_ext = 5.0):
    import numpy as np

    if epirec['Hardware_reference']=='Fpz':
        chains = EEG_chain_bipolar
    elif epirec['Hardware_reference']=='bipolar':
        chains = None
    else:
        raise Exception('Unknown reference')
    t, data, names = LoadEDF(epirec['EEG_file'], t0, t1, chains, None, None)

    fs=float(len(t)-1)/float(t[-1]-t[0])

    #elongate the data
    t_init_idx = int(float(filter_ext)*fs)
    if t_init_idx>len(t):
        t_init_idx=len(t)
    d_init = np.take(data,range(0, t_init_idx), axis=data.ndim-1)
    data_long = np.concatenate((d_init,data),axis=data.ndim-1)

    import SignalTools

    # notch filter
    notch_filter=(np.zeros(2,dtype=float)+epirec['Software_filters']+np.array([-1.0,1.0]))/fs*2.0
    data_long = SignalTools.filter(data_long,notch_filter,3,'bandstop')

    # bandpass filter
    bandpass_filter=np.array(epirec['Hardware_filters'])/fs*2.0
    data_long = SignalTools.filter(data_long,bandpass_filter,5,'bandpass')

    # cut the data
    data = np.take(data_long, range(t_init_idx,data_long.shape[data_long.ndim-1]),
            axis=data.ndim-1)

    # offset and amplitude adjusting
    if offset is None:
        data = data - np.atleast_2d(np.mean(data,-1)).T
    if amplitude is None:
        data = SignalTools.AutoGain(t,data,100.0,showProg=False)

    EEG=dict()
    EEG['t']=t
    EEG['data']=data
    EEG['names']=names
    EEG['fs']=fs
    EEG['bandpass_filter']=bandpass_filter
    EEG['notch_filter']=notch_filter
    EEG['file']=epirec['EEG_file']
    return EEG


if __name__ == "__main__01":
    from pylab import *
    from Utility import Time
    from PlotTools import PlotEEG
    filename="../SharedData/PublicEEG/CHIMIC_EEG_DATA.edf"
    t0 = Time(minutes=0, seconds=0.0)
    t1 = Time(minutes=3, seconds=0.0)
    offset = 1.0
    amplitude = 0.01
    Wn = [0.02,0.30]
    filter_order = 6

    #chan_names=[("Fp1","F7"),("F7","T3"),("T3","T5"),("T5","O1"),("C3","P3")];
    times,data,chan_labels = LoadEDF(filename,t0,t1,EEG_chain_bipolar,Wn,filter_order)
    PlotEEG(times,data,channels=chan_labels,offset=offset,amplitude=amplitude)

if __name__ == "__main__":
    from Utility import Time
    #rec_name = 'Baptist/PAT16_12_27_2013.txt'
    rec_name = 'ABNORMAL/MK'
    t0 = Time(minutes=0)
    t1 = Time(minutes=1)
    EEG=LoadCateEEG(rec_name,segsep=-1)
    from pylab import *
    from PlotTools import PlotEEG
    #PlotTools.PlotEEG(EEG['t'],EEG['data'],channels=EEG['names'])
    #basepath="../SharedData/Hoda/ABNORMAL"
    #fil="../SharedData/Hoda/Abnormal filtered/Abnormal_filtered.mat"
    #from scipy.io import loadmat
    #matlab=loadmat(fil)
    #P=np.array(matlab['P5'])
    #t=np.arange(0,P.shape[1],dtype=float)/200
    #fig,axes=subplots(2,1)
    PlotEEG(EEG['t'],EEG['data'])
    #PlotEEG(t,P,figaxes=axes[1])
    show()

