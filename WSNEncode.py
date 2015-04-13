import os,sys
import numpy as np
from Utility import *
import Errors
import nest

class Models:

    _ModelTemplate = {
            'Init': {
                # field_name    : default  , acceptable  , len, content_type, convert_type
                'spike_times'   : ( [1.0]  , (list,tuple), 1  , float       , float      ),
                'spike_weights' : ( [10.0] , (list,tuple), 1  , float       , float      )
                },
            'Clk' : {
                # field_name : default  , acceptable  , len , content_type, convert_type
                'C_m'        : ( 1.0    , (int, float), None, None        , float       ),
                'E_L'        : ( 0.0    , (int, float), None, None        , float       ),
                'I_e'        : ( 0.0    , (int, float), None, None        , float       ),
                'V_m'        : ( 0.0    , (int, float), None, None        , float       ),
                'V_reset'    : ( 0.0    , (int, float), None, None        , float       ),
                'V_th'       : ( 1.0    , (int, float), None, None        , float       ),
                't_ref'      : ( 20.0   , (int, float), None, None        , float       ),
                'tau_m'      : ( 1.0    , (int, float), None, None        , float       )
                },
            'WSN' : {
                # field_name : default  , acceptable  , len , content_type, convert_type
                'C_m'        : ( 1.0    , (int, float), None, None        , float       ),
                'E_L'        : ( 0.0    , (int, float), None, None        , float       ),
                'I_e'        : ( 0.0    , (int, float), None, None        , float       ),
                'Sigma'      : ( ['const',1.0], (list,tuple,Distribution), 
                    None, None, Distribution),
                'D_Int'         : ( ['const',50.0], (list,tuple,Distribution),
                    None, None, Distribution),
                'V_m'        : ( 0.0    , (int, float), None, None        , float       ),
                'V_reset'    : ( -20.0  , (int, float), None, None        , float       ),
                'V_th'       : ( -1.0   , (int, float), None, None        , float       ),
                't_ref'      : ( 1.0    , (int, float), None, None        , float       ),
                'tau_m'      : ( 30.0   , (int, float), None, None        , float       ), 
            },
            'Clk_Clk' : {
                # field_name: default  , acceptable  , len , content_type, convert_type
                'weight'    : (10.0    , (int, float), None, None        , float       ),
                'delay'     : (100.0   , (int, float), None, None        , float       )
                },
            'Clk_WSN' : {
                # field_name   : default  , acceptable  , len , content_type, convert_type
                'weight'       : (1.0     , (int, float), None, None        , float      ),
                'delay'        : (1.0     , (int, float), None, None        , float      )
                #'receptor_type': (1       , int         , None, None        , None       ),
                },
            'Clk_Spk' : {
                # field_name   : default  , acceptable  , len , content_type, convert_type
                'weight'       : (1.0     , (int, float), None, None        , float      ),
                'delay'        : (1.0    , (int, float), None, None        , float      ),
                }
            }

    def __init__(self,**model):
        # __init__(model=allmodes, **special_model)
        if 'model' in model.keys():
            self.__init__(**model['model'])
            del model['model']

        dummyModelIn = dict.fromkeys(model)
        for key,item in Models._ModelTemplate.iteritems():
            if key in model:
                self.__dict__[key] = TemplateDict(Models._ModelTemplate[key],**model[key])
                del dummyModelIn[key]
            else:
                self.__dict__[key] = TemplateDict(Models._ModelTemplate[key])

        if bool(dummyModelIn):
            raise Errors.ModelUndefinedError(', '.join(dummyModelIn.keys()))

    def setParam(self,**model):
        mt = Models._ModelTemplate
        allModelNames = mt.keys()
        for key in model.keys():
            if key not in allModelNames:
                raise Errors.ModelUndefinedError(key)
            if model[key] is None:
                # reset to default
                self.__dict__[key] = TemplateDict(mt[key])
            else:
                self.__dict__[key].setItems(**model[key])

    def D(self,name):
        # return the nest compatible dict
        if not isinstance(name,str):
            raise Errors.ParamTypeError('name','str')
        if name not in self.__dict__.keys():
            raise Errors.ModelUndefinedError(name)
        return self.__dict__[name].genNestDict()

class InputSignal:
    def __init__(self,t,I,I_size):
        # t is the time table, y is the signal amplitude
        if not isinstance(t,(list,tuple,np.ndarray)):
            raise Errors.ParamTypeError('t', 'list, tuple or numpy array')
        if not isinstance(I,(list,tuple,np.ndarray)):
            raise Errors.ParamTypeError('I', 'list, tuple or numpy array')
        if not isinstance(I_size,(list,tuple,np.ndarray)):
            raise Errors.ParamTypeError('I_size', 'list, tuple or numpy array')

        if not isinstance(t, np.ndarray):
            self.t = np.array(t,dtype=float)
        else:
            self.t = t

        if not isinstance(I, np.ndarray):
            self.I = np.array(I,dtype=float)
        else:
            self.I = I

        if not isinstance(I_size,np.ndarray):
            self.I_size = np.array(I_size,dtype=int)
        else:
            self.I_size = I_size

        if self.t.ndim != 1:
            raise Errors.ParamDimError('t', 1)
        else:
            self.N_sample = len(self.t)

        if self.I_size.ndim != 1:
            raise Errors.ParamDimError('I_size', 1)

        if np.prod(self.I_size) == 1:
            if self.I.ndim != 1:
                raise Errors.ParamDimError('I', 1)
            if len(self.I) != self.N_sample:
                raise Errors.ParamSizeError('I', [self.N_sample])
        else:
            if self.I.ndim != len(self.I_size) + 1:
                raise Errors.ParamSizeError('I', len(self.I_size) +1)
            if self.I.shape[:-1] != self.I_size:
                from copy import deepcopy
                corrSize=deepcopy(I_size)
                corrSize.append(self.N_sample)
                raise Errors.ParamSizeError('I', corrSize)

        self.h = t[1]-t[0]
        self.length = t[-1]-t[0]+self.h
        self.freq = 1.0 / self.h
        self.I = self.I.reshape( (np.prod(I_size),self.N_sample) )

class WSNnet:
    _NetParamTemplate = {
            #'name_of_param': (default, (allowed types), size, content, conv_type)
            'sim_h'         : ( 0.01        , (float,int)    , None, None   , float        ),
            'save_h'        : ( 0.01        , (float,int)    , None, None   , float        ),
            'Size_In'       : ( [1]         , (list,tuple)   , None, int    , None         ),
            'N_Sig_Per_Ch'  : ( 100         , int            , None, None   , None         ),
            'N_Per_Sig'     : ( 1           , int            , None, None   , None         ),
            'Vth'           : ( -1.0        , (float,int)    , None, None   , float        ),
            'Vc'            : ( -5.0        , (float,int)    , None, None   , float        ),
            'Tau'           : ( 40.0        , (float,int)    , None, None   , float        ),
            'Sigma'         : (['const',1.0], (list,tuple,Distribution)
                , None   , None   , Distribution ),
            'T_Clk'         : ( 100.0       , (float,int)    , None, None   , float        ),
            #'T_Enc'         : ( 45.0        , (float,int)    , None, None   , float        ),
            #'D_Enc'         : ( 50.0        , (float,int)    , None, None   , float        ),
            'D_Int'         : (['const',1.0], (list,tuple,Distribution)
                , None   , None   , Distribution ),
            'En_Trace'      : ( False       , bool           , None, None   , None         ),
            'TraceNodes'    : ( []          , (list,tuple)   , None, int    , None         )
            }
    _NestParamTemplate = {
            #'WSN_Model_Name' : ( 'iaf_freq_sensor_v2', str , None, None, None ),
            'WSN_Model_Name' : ( 'wsn_hermitian_2', str , None, None, None ),
            'MyModule_Name'  : ( 'mymodule'       , str , None, None, None ),
            #'Syn_Model_Name' : ( ['CLK_WSN_synapse', 'SRC_WSN_synapse'],
                                                    #list, None, None, None ),
            'verbosity'      : ( 'M_ERROR'        , str , None, None, None ),
            'print_time'     : ( True             , bool, None, None, None )
            }

    def __init__(self, **kwargs):
        self.NetP = TemplateDict(WSNnet._NetParamTemplate,**kwargs)
        self.NestP = TemplateDict(WSNnet._NestParamTemplate)
        self.Mod = Models()

        # set Vth, Vc, Tau, Sigma, T_Clk, T_Enc, T_Delay if needed
        # update wsn model
        P = self.NetP
        updates = P.genDict() 
        if 'save_h' not in kwargs.keys():
            updates['save_h'] = P.sim_h

        if 'Sigma' not in kwargs.keys():
            updates['Sigma'] = ['linear',P.sim_h*2.0,P.T_Enc]

        if 'Vc' not in kwargs.keys():
            updates['Vth'] = P.Vth

        self.setParam(**updates)

    def setParam(self,**param):
        self.NetP.setItems(**param)
        P = self.NetP
        M = self.Mod

        U = {} 

        if ('Vc' not in param) and (
                'Vth' in param or 'T_Clk' in param.keys()):
            if P.Vth < 0:
                P.Vc = np.asscalar(P.Vth * np.exp(P.T_Clk/P.Tau))
            elif P.Vth > 0:
                P.Vc = 0.0
            else:
                raise Errors.NetParamError("Vth should not equal to zero")

        if 'Sigma' in param or 'N_Sig_Per_Ch' in param:
            P.Sigma.rand(P.N_Sig_Per_Ch)
            U['Sigma']=P.Sigma
        if 'D_Int' in param or 'N_Per_Sig' in param:
            P.D_Int.rand(P.N_Per_Sig)
            U['D_Int'] = P.D_Int
        if 'Vth' in param:
            U['V_th']=P.Vth
        if 'Tau' in param:
            U['tau_m'] = P.Tau
        if ('Vth' in param) or ('T_Clk' in param) or ('Vc' in param):
            U['V_reset'] = P.Vc

        if bool(U):
            M.setParam(WSN=U)

        if 'T_Clk' in param:
            M.setParam(Clk_Clk = {'delay': P.T_Clk})

    def CreateNet(self,rebuild=True,SaveFile=None):
        from copy import deepcopy
        from os import path
        if SaveFile is not None:
            self.SaveFile = path.abspath(SaveFile)
            rec_dict = {
                    'to_file'        : True,
                    'to_memory'      : False,
                    }
        else:
            self.SaveFile = None
            rec_dict = {
                    'to_file'        : False,
                    'to_memory'      : True,
                    }

        #setup nest kernel
        nest.set_verbosity(self.NestP.verbosity)
        if self.NestP.WSN_Model_Name not in nest.Models():
            nest.Install(self.NestP.MyModule_Name)
        if rebuild:
            nest.ResetKernel()
        nest.SetStatus([0],{'resolution':self.NetP.sim_h})

        Mod = self.Mod
        NetP = self.NetP
        NestP = self.NestP

        #Create Nodes
        NetP.N_Per_Channel = NetP.N_Sig_Per_Ch * NetP.N_Per_Sig
        self.N_In = np.asscalar(np.prod(NetP.Size_In))
        self.WSN_Size = list(deepcopy(NetP.Size_In))
        self.WSN_Size.append(NetP.N_Per_Channel)
        self.N_WSN = self.N_In * NetP.N_Per_Channel

        self.node_init = nest.Create("spike_generator", 1 ,params = Mod.D('Init'))
        self.node_clk  = nest.Create("iaf_psc_alpha"  , 1 ,params = Mod.D('Clk'))
        self.node_src = nest.Create("step_current_generator", self.N_In)
        self.node_src_grps = np.reshape(self.node_src, NetP.Size_In)
        self.node_wsn = nest.Create(NestP.WSN_Model_Name, self.N_WSN, params = Mod.D('WSN'))
        self.node_wsn_chan_grps = np.reshape(self.node_wsn,(self.N_In,NetP.N_Sig_Per_Ch,NetP.N_Per_Sig))
        self.node_wsn_grps = np.reshape(self.node_wsn, self.WSN_Size)

        self.spk_wsn = nest.Create("spike_detector",1,params=rec_dict)
        self.spk_clk = nest.Create("spike_detector",1,params=rec_dict)
        if SaveFile is not None:
            nest.SetStatus(self.spk_wsn,{'label':'_'.join([self.SaveFile,'spk_wsn'])})
            nest.SetStatus(self.spk_clk,{'label':'_'.join([self.SaveFile,'spk_clk'])})

        #Create multimeter nodes for Tracing
        if NetP.En_Trace:
            if not NetP.TraceNodes:
                rec_nodes = self.node_wsn
            else:
                rec_nodes = NetP.TraceNodes
            rec_fields = nest.GetStatus([rec_nodes[0]],"recordables")[0]
            self.rec_wsn = nest.Create("multimeter",len(rec_nodes),params={
                'record_from'   : rec_fields,
                'interval'      : NetP.save_h
                })
            self.rec_clk = nest.Create("multimeter",1,params={
                'record_from'   : ['V_m'],
                'interval'      : NetP.save_h
                })
            if SaveFile is not None:
                nest.SetStatus(self.rec_wsn,rec_dict)
                nest.SetStatus(self.rec_clk, rec_dict)
                nest.SetStatus(self.rec_wsn,{'label':'_'.join([self.SaveFile,'rec_wsn'])})
                nest.SetStatus(self.rec_clk,{'label':'_'.join([self.SaveFile,'rec_clk'])})

        #Connect Nodes
        nest.Connect(self.node_init,self.node_clk)
        nest.Connect(self.node_clk,self.node_clk,syn_spec=Mod.D('Clk_Clk'))
        nest.Connect(self.node_clk, self.spk_clk,syn_spec=Mod.D('Clk_Spk'))
        nest.Connect(self.node_wsn, self.spk_wsn,'all_to_all')

        #Connect Clk to WSN
        nest.Connect(self.node_clk, self.node_wsn,'all_to_all',syn_spec=Mod.D('Clk_WSN'))

        #Connect source nodes to WSN nodes and setup the sigmas and d_ints
        for src, wsn_sig_grps in zip(self.node_src, self.node_wsn_chan_grps):
            # setup sigma and d_int
            sigDist = Mod.WSN.Sigma
            delayDist = Mod.WSN.D_Int
            for w_sig, s in zip(wsn_sig_grps, sigDist.data):
                for w,d in zip(w_sig, delayDist.data):
                    nest.SetStatus([w],{'Sigma':s, 'D_Int': d})

            # Connect Source to WSN
            nest.Connect([src],np.ravel(wsn_sig_grps).tolist(),'all_to_all')

        if NetP.En_Trace:
            # Connect multimeter nodes for tracing
            nest.Connect(self.rec_wsn,rec_nodes)
            nest.Connect(self.rec_clk,self.node_clk)

    def Sim(self, t_in, I_in, reset=True):
        if 'node_wsn' not in self.__dict__:
            raise Errors.NetNotInitError()
        NetP = self.NetP
        #self.signal=InputSignal(t_in,I_in,NetP.Size_In)
        #S = self.signal
        S = InputSignal(t_in, I_in, NetP.Size_In)
        nest.set_verbosity(self.NestP.verbosity)
        nest.SetStatus([0],{
            'print_time': self.NestP.print_time,
            'overwrite_files' : True
            })
        # apply signal to source nodes
        for s,Is in zip(self.node_src,S.I):
            nest.SetStatus([s],{
                'amplitude_times' : S.t-t_in[0],
                'amplitude_values' : Is
                })
        nest.Simulate(S.length)

    def get_sigmas(self,nodes=None):
        if nodes is None:
            nodes = self.node_wsn

        sigs = nest.GetStatus(nodes,"Sigma")
        self.sigs = dict(zip(nodes,sigs))
        return sigs

    def get_delays(self,nodes=None):
        if nodes is None:
            nodes = self.node_wsn
        delays = nest.GetStatus(nodes,'D_Int')
        self.delays = delays
        return delays

        #syns = nest.GetConnections(
                #target=self.node_wsn,
                #synapse_model=self.NestP.Syn_Model_Name[0] )
        #self.int_delays = nest.GetStatus(syns,'delay')
        #syns = nest.GetConnections(
                #target=self.node_wsn,
                #synapse_model=self.NestP.Syn_Model_Name[1] )
        #self.enc_delays = nest.GetStatus(syns,'delay')
        #return nest.GetStatus(syns,'delay')

    def get_WSN_spikes(self,t0=None, t1=None, nodes=None):
        if nodes is None:
            nodes = self.node_wsn
        if self.SaveFile is None:
            wsn_events = nest.GetStatus(self.spk_wsn,'events')[0]
            times = wsn_events['times']
            senders = wsn_events['senders'].astype(int)
            if t0 is None:
                t0msec = min(times)
            else:
                if isinstance(t0,Time):
                    t0msec=t0.as_msec()
                else:
                    t0msec=t0*1000.0
            if t1 is None:
                t1msec = max(times)
            else:
                if isinstance(t1,Time):
                    t1msec=t1.as_msec()
                else:
                    t1msec=t1*1000.0
            use_idx = np.array([ t>=t0msec and t<=t1msec and s in nodes
                    for t,s in zip(times,senders) ])
            return times[use_idx],senders[use_idx]
        else:
            spkfile = self.get_spk_filename('wsn')
            senders, times = self.load_spk_fromfile(spkfile,t0,t1)
            return times,senders
            #wsn_events = np.loadtxt(spkfile)
            #senders = wsn_events[:,0].astype(int)
            #times = wsn_events[:,1]

    def get_clk_times(self,t0=None,t1=None):
        if self.SaveFile is None:
            clks = nest.GetStatus(self.spk_clk, 'events')[0]['times']
            if t0 is None:
                t0msec = min(clks)
            else:
                if isinstance(t0,Time):
                    t0msec = t0.as_msec()
                else:
                    t0msec = t0*1000.0
            if t1 is None:
                t1msec = max(clks)
            else:
                if isinstance(t1,Time):
                    t1msec = t1.as_msec()
                else:
                    t1msec = t1*1000.0
            return clks[np.logical_and(clks>=t0msec,clks<=t1msec)]
        else:
            spkfile = self.get_spk_filename('clk')
            senders, clks = self.load_spk_fromfile(spkfile,t0,t1)
            return clks

    def get_spk_filename(self,ftype='wsn'):
        if self.SaveFile is not None:
            import os
            from fnmatch import fnmatch
            fd,fn = os.path.split(self.SaveFile)
            flist = os.listdir(fd)
            if ftype=='wsn':
                fpat = '_'.join([fn,'spk_wsn*'])
            elif ftype=='clk':
                fpat = '_'.join([fn,'spk_clk*'])
            else:
                raise Errors.NetPathInexistError

            spk_f = [ f for f in flist if fnmatch(f, fpat) ]
            if not spk_f:
                raise Errors.NetPathInexistError
            else:
                return os.path.join(fd,spk_f[-1])
        else:
            raise Errors.NetPathInexistError

    def load_spk_fromfile(self,filename,t0=None,t1=None):
        from Utility import Time
        f = open(filename,'r')
        import csv
        reader = csv.reader(f,delimiter='\t')
        senders = list()
        times = list()
        if t0 is not None:
            if not isinstance(t0,Time):
                t0 = Time(seconds = t0)
            t0 = t0.as_msec()
        if t1 is not None:
            if not isinstance(t1,Time):
                t1 = Time( seconds=t1 )
            t1 = t1.as_msec()

        for r in reader:
            t = float(r[1])
            sv_this = True
            if t0 is not None:
                sv_this = t>=t0
            if t1 is not None:
                sv_this = sv_this and (t<=t1)
            if sv_this:
                senders.append(int(r[0]))
                times.append(t)
        return np.array(senders,dtype=int),np.array(times,dtype=float)

    def get_WSN_spike_phases(self,t0=None,t1=None,nodes=None,inverse=False,group=False):
        from bisect import bisect_right
        clk = self.get_clk_times(t0,t1)
        if nodes is None:
            nodes = self.node_wsn
            grp_shape = np.array(self.node_wsn_chan_grps.shape)
        else:
            grp_shape = np.array(np.array(nodes).shape)
            nodes = np.ravel(nodes)

        times, senders = self.get_WSN_spikes(t0,t1,nodes)

        if self.SaveFile is None:
            self.get_delays(nodes) 

        #int_delays = self.int_delays
        #enc_delays = self.enc_delays
        T_Clk = self.NetP.T_Clk

        nodes_idx = dict(zip(nodes,range(len(nodes))))

        spk_phase = np.full( (len(nodes), len(clk)), np.NaN )
        for t, s in zip(times,senders):
            nid = nodes_idx[s]

            #enc_duration = int_delays[nid] + enc_delays[nid] + T_Enc
            if t>=clk[0] and t<=clk[-1]+T_Clk:
                cid = bisect_right(clk,t)-1
                if inverse:
                    phase = clk[cid] + T_Clk - t
                else:
                    phase = t - clk[cid]
                spk_phase[nid,cid] = phase
        if group:
            return clk[:-1], np.reshape(spk_phase[:,:-1], np.append(grp_shape,len(clk)-1))
        else:
            return clk[:-1], spk_phase[:,:-1]

if __name__=='__main__':
    from pylab import *
    from nest.raster_plot import from_device as raster
    import scipy.signal as sigpack
    sys.path.append("../Tools")
    import FileIO
    import PlotTools

    t_in, I_in = FileIO.LoadVoiceData("an4_clstk/fash/an253-fash-b.raw")
    #t_in, I_in = PF.LoadVoiceData("an4_clstk/fjam/cen1-fjam-b.raw")
    t_in = t_in * 10000.0
    In_Amp = 5.0
    I_in = I_in * In_Amp

    network = WSNnet(
             Sigma=['linear', 0.2, 10.0],
             Size_In=[1],
             N_Sig_Per_Ch=10,
             Vth=-1.0,
             Tau=45.0,
             D_Int=['const', 25.0],
             En_Trace = False
             )
    network.CreateNet()
    network.Sim(t_in,I_in)
    #times, senders = network.get_WSN_spikes(3000,5000)
    #events = nest.GetStatus([network.rec_wsn[3]],'events')[0]
    clk, phase = network.get_WSN_spike_phases(group=False)
    sigmas = network.get_sigmas()
    #PlotTools.pcolor(clk,sigmas,phase,y_ticks=True,y_ticklabels=['a','b','c'])
    #PlotTools.raster(times,nodes)

    show()

    #OF.PlotVoiceResult(u_events,clk,sigmas,u_nodes,t_in,I_in,P.syn.Clk_Enc['delay'],P.neu.U['Ti'])

    ion()

