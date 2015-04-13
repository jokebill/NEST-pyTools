import numpy as np
import Errors

default_plot_font = {
        'family': 'FreeSerif',
        'size'  : 10
        }

def GetSubplotPos(fig):
    p = fig.subplotpars
    print "left:\t{0}".format(p.left)
    print "right:\t{0}".format(p.right)
    print "bottom:\t{0}".format(p.bottom)
    print "top:\t{0}".format(p.top)
    print "hspace:\t{0}".format(p.hspace)
    print "wspace:\t{0}".format(p.wspace)

def PlotEEG(time, data, **kwargs):
    from Utility import TemplateDict,Time
    import SignalTools as SigProc
    from types import NoneType

    import matplotlib
    import matplotlib.pyplot as plt

    # **kwargs includes:
    kwTemplate = {
            'channels'  : ([],    (list,tuple), None, str,  None ),
            'offset'    : (1.0,   (int, float), None, None, float),
            'amplitude' : (3.0,   (int, float), None, None, float),
            'normalize' : (False, bool,         None, None, None ),
            'font'      : (None,  (NoneType,dict), None, None, None),
            'figsize'   : (None,  (NoneType,tuple,list), None, None, None),
            'figaxes'   : (None,  None, None, None, None)
            }

    if not isinstance(time,np.ndarray):
        time = np.array(time,dtype=float)
    if time.ndim<>1:
        raise Errors.ParamDimError('time',1)
    Nx = len(time)

    if not isinstance(data,np.ndarray):
        data = np.array(data,dtype=float)

    if data.ndim>2 or data.ndim<1:
        raise Errors.ParamDimError('data','1 or 2')
    if data.ndim==1:
        if len(data) <> Nx:
            raise Errors.ParamSizeError('data',(Nx,1))
        Ny = 1
    else:
        Ny = data.shape[0]
        if data.shape[1] <> Nx:
            raise Errors.ParamSizeError('data',(Nx,Ny))

    opts = TemplateDict(kwTemplate, **kwargs)

    if Ny>1:
        if bool(opts.channels):
            if len(opts.channels)<>Ny:
                raise Errors.ParamSizeError('channels',Ny)
        offset_vec = np.arange(Ny,0,-1,dtype=float)*opts.offset;
        if opts.normalize:
            data=SigProc.normalize(data,opts.amplitude)
        else:
            data=data*opts.amplitude
        data=data+np.atleast_2d(offset_vec).T

    if opts.font is not None:
        matplotlib.rc('font',**opts.font)

    if opts.figaxes is None:
        if opts.figsize is None:
            fig, ax = plt.subplots()
        else:
            fig, ax = plt.subplots(figsize=opts.figsize)
    else:
        ax = opts.figaxes

    lines=ax.plot(time,data.T)

    if Ny>1:
        ax.set_ylim(offset_vec[-1]-opts.offset,offset_vec[0]+opts.offset)
        ax.set_yticks(offset_vec)
        if bool(opts.channels):
            ax.set_yticklabels(opts.channels)

    ax.set_xlim(time[0],time[-1])


    def second_to_timestr(seconds,loc):
        tmptime=Time(seconds)
        return repr(tmptime)

    ax.get_xaxis().set_major_formatter(matplotlib.ticker.FuncFormatter(second_to_timestr))
    if opts.figaxes is None:
        plt.show()
        return fig
    else:
        if Ny>1:
            return lines,offset_vec
        else:
            return lines

def raster(times, senders, sigmas=None,
        barlen=None, Ax=None, xylim=True,**plotspec):
    import matplotlib.pyplot as plt


    if sigmas is None:
        Y = senders
    else:
        Y = [ sigmas[n] for n in senders ]

    minY = min(Y)
    maxY = max(Y)
    if minY == maxY:
        dY_half = 0.5
    else:
        dY_half = float(maxY-minY)/float(len(set(Y))-1)/2.0

    if Ax is None:
        fig, Ax = plt.subplots(1,1)

    if barlen is None:
        if 'color' not in plotspec:
            plotspec['color']=(0,0,1)
        if 's' not in plotspec:
            plotspec['s']=10
        Ax.scatter(times, Y, **plotspec)
    else:
        if 'color' not in plotspec:
            plotspec['color']=(0,1,0)
        if 'linewidth' not in plotspec:
            plotspec['linewidth']=0.5

        dY_half = dY_half * barlen
        xx = np.ravel([ [ t, t, np.NaN ] for t in times ])
        yy = np.ravel([ [y-dY_half, y+dY_half, np.NaN ] for y in Y])
        Ax.plot(xx,yy,**plotspec)

    if xylim:
        minT=min(times)
        maxT=max(times)
        dT = float(maxT-minT)/float(len(set(times))-1)
        Ax.set_xlim(minT-dT*20,maxT+dT*20)
        Ax.set_ylim(minY-dY_half*2, maxY+dY_half*2)

def pcolor(x,y,C,
        x_label='time / sec', y_label='Sigma / ms',
        x_range=None, y_range=None,
        y_ticks = None, y_ticklabels = None, Ax=None,
        **spec):
    import matplotlib.cm as colormaps
    import matplotlib.pyplot as plt
    x = np.append(x, 2.0*x[-1]-x[-2])
    y = np.append(y, 2.0*y[-1]-y[-2])
    masked_C = np.ma.array(C,mask=np.isnan(C))
    
    if 'cmap' not in spec:
        cmap = colormaps.hot
        cmap.set_bad('w', 1.0)
        spec['cmap'] = cmap

    if Ax is None:
        newFig = True
        fig, Ax = plt.subplots(1,1)
    else:
        newFig = False

    meshes = Ax.pcolor(x,y,masked_C,**spec)

    if x_range is None:
        x_min = min(x)
        x_max = max(x)
    else:
        x_min, x_max = x_range

    if y_range is None:
        y_min = min(y)
        y_max = max(y)
    else:
        y_min, y_max = y_range

    Ax.set_xlim(x_min,x_max)
    Ax.set_ylim(y_min,y_max)

    if x_label is not None:
        Ax.set_xlabel(x_label)
    if y_label is not None:
        Ax.set_ylabel(y_label)
    if y_ticks is not None:
        if y_ticks is True:
            y_ticks = (y[1:]+y[:-1])/2.0
        Ax.set_yticks(y_ticks)
    if y_ticklabels is not None:
        Ax.set_yticklabels(y_ticklabels)

    if newFig:
        fig.colorbar(meshes)
    else:
        return meshes

def EEG_Comb(t, data, names, network, t0=None, t1=None, amplitude=0.5, iis=None, **spec):
    ## t0, t1 unit: second if is type float
    import Utility as utl
    from bisect import bisect

    if t0 is None:
        t0i = 0
        t0 = t[0]
    else:
        t0i = bisect(t, float(t0))

    if t1 is None:
        t1i = len(t)
        t1 = t[-1]
    else:
        t1i = bisect(t, float(t1))


    t_plot = utl.DataSegment(t, t0i, t1i)
    d_plot = utl.DataSegment(data,t0i,t1i)

    t0_nest = utl.Time(seconds=float(t0)-t[0])
    t1_nest = utl.Time(seconds=float(t1)-t[0])

    clks, phase = network.get_WSN_spike_phases(t0=t0_nest, t1=t1_nest, group=True, inverse=False)
    clks_sec = clks/1000.0+t[0]
    #N_ch,N_sig,N_delay,N_clk = phase.shape
    #N_p_ch = N_sig * N_delay
    #y_values = np.arange(0, N_ch, 1.0/N_sig, dtype=float)[::-1]

    #phase = np.reshape(phase, (N_ch*N_sig, N_delay, N_clk))

    #np.place(phase,np.isnan(phase),np.inf)
    #phase = np.amin(phase, -2)
    #np.place(phase,np.isinf(phase),np.nan)

    N_ch, N_delay, N_clk = phase.shape
    y_values = np.arange(0, N_ch, 1.0/N_delay, dtype=float)[::-1]
    y_ticks = np.arange(0, N_ch, dtype=float)+0.5

    phase = np.reshape(phase, (N_ch*N_delay,N_clk))


    if 'font' in spec:
        font = spec['font']
    else:
        font = default_plot_font

    if 'figsize' in spec:
        figsize = spec['figsize']
    else:
        figsize = [3.45*1.5,5.5*1.5]

    if 'dpi' in spec:
        dpi = spec['dpi']
    else:
        dpi = 100

    if 'plt_pos' in spec:
        plt_pos = spec['plt_pos']
    else:
        plt_pos={
                'left': 0.1,
                'right': 0.8,
                'bottom': 0.1,
                'top' : 0.9,
                'hspace': 0.2
                }

    if 'cbar_pos' in spec:
        cbar_pos = spec['cbar_pos']
    else:
        cbar_pos = [0.85,0.1,0.03,0.35]

    if 'c_ticks' in spec:
        c_ticks = spec['c_ticks']
    else:
        c_ticks = np.arange(0,100.1,20)

    import matplotlib
    matplotlib.rc('font', **font)

    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(nrows=2, ncols=1,figsize=figsize, dpi=dpi)

    PlotEEG(t_plot, d_plot, channels = names, amplitude=amplitude,figaxes=axes[0]) 

    pShapes=pcolor(clks_sec, y_values, phase, 
            x_label='', x_range = [float(t0),float(t1)],
            y_label='', y_ticks = y_ticks, y_ticklabels=names[::-1],
            Ax=axes[1])
    if iis is not None:
        from matplotlib.patches import Rectangle as Rect
        for iis_t0,iis_t1 in iis:
            if (iis_t0>=t0 and iis_t0<=t1) or (iis_t1>=t0 and iis_t1<=t1):
                rec_width = float(iis_t1)-float(iis_t0)
                for ax in axes:
                    ylim = ax.get_ylim()
                    rec_loc = (float(iis_t0),ylim[0])
                    rec_height = ylim[1]-ylim[0]
                    ax.add_patch(Rect(rec_loc,rec_width,rec_height,
                        edgecolor='none', facecolor = 'red', alpha=0.3
                        ))

    from Utility import Time
    def second_to_timestr(seconds,loc):
        tmptime=Time(seconds)
        return repr(tmptime)
    from matplotlib.ticker import FuncFormatter
    axes[1].get_xaxis().set_major_formatter(FuncFormatter(second_to_timestr))
    fig.subplots_adjust(**plt_pos)
    cbar_ax = fig.add_axes(cbar_pos)
    cb=fig.colorbar(pShapes,cax=cbar_ax)
    cb.set_ticks(c_ticks)
    plt.show()
    return fig

def PlotNode(data,key=None,t0=None,t1=None):
    from Utility import Time
    from bisect import bisect
    import matplotlib.pyplot as plt
    t = data['times']
    if t0 is None:
        t0i = 0
    else:
        t0i = bisect(t,Time(t0).as_msec())
    if t1 is None:
        t1i = len(t)
    else:
        t1i = bisect(t,Time(t1).as_msec())+1

    plt.figure()
    if key is None:
        N=len(data)-2
        i=1
        for key,value in data.iteritems():
            if key not in ['times','senders']:
                plt.subplot(N,1,i)
                plt.plot(t[t0i:t1i],value[t0i:t1i])
                plt.xlim(t[t0i],t[t1i-1])
                plt.xlabel(key)
                i=i+1
    else:
        assert(key in data.keys())
        plt.plot(t,data[key])
    plt.show()

def PlotWavelet(t,y,sigmas,clk=None,wavelet=None,T_Int=45.0,Offset=1.0,NewFigure=True):
    from scipy import signal
    import matplotlib.pyplot as plt
    from bisect import bisect
    import numpy as np
    h = t[1]-t[0]
    slist = sigmas/h
    if not wavelet:
        wavelet = signal.ricker
    dwt = signal.cwt(y,wavelet,slist)
    if clk is None:
        tlist = np.append(t,t[-1]+h)
        dwt = np.abs(dwt)
    else:
        tlist = clk + T_Int / 2.0 + Offset 
        tids = [ bisect(t, tr) for tr in tlist ]
        dwt = np.abs(np.array([ dwt[:,i] for i in tids ])).T
        tlist = np.append(tlist,tlist[-1]+tlist[1]-tlist[0])
    sigmas = np.append(sigmas,sigmas[-1]+sigmas[1]-sigmas[0])
    if NewFigure:
        plt.figure()
    plt.pcolormesh(tlist,sigmas,dwt,cmap=plt.get_cmap('hot'))
    if clk is None:
        plt.xlim(tlist[0],tlist[-1])
    else:
        plt.xlim(clk[0],clk[-1])
    plt.ylim(sigmas[0],sigmas[-1])
    #plt.xlabel("Time / ms")
    #plt.ylabel("Sigma / ms")
    if NewFigure:
        plt.colorbar()
        plt.show()

def PlotVoice(data):
    pass
