import Errors

class TemplateDict:
    def __init__(self,templateDict,**items):
        #Template should always have the following format:
        #field_name(str) : default value, acceptable form, length, content type, convert type

        #ToDo: check template format

        self._template = templateDict
        for key,tmpl in templateDict.iteritems():
            if key not in items:
                items[key] = tmpl[0]
            self.__dict__[key] = None
        self.setItems(**items)
    def __getitem__(self,key):
        return self.__dict__[key]
    def __setitem__(self,key,val):
        self.__dict__[key]=val

    def getItems(self):
        return [ v for v in self.__dict__.keys() if not str.startswith(v,'_') ]
    def setItems(self,**items):
        # modify an old dictionary
        dummyNew = dict.fromkeys(items)
        selfKeys = self.getItems() 
        for key,val in items.iteritems():
            if key not in selfKeys:
                raise Errors.TemplUndefinedFieldError(key)
            tmpl = self._template[key]
            if tmpl[1] is not None:
                if not isinstance(val,tmpl[1]):
                    raise Errors.TemplFieldTypeError(key)
            if tmpl[2] is not None:
                if len(val)!=tmpl[2]:
                    raise Errors.TemplFieldLengthError(key,tmpl[2])
                for v in val:
                    if not isinstance(v,tmpl[3]):
                        raise Errors.TemplFieldContentError(key)
            del dummyNew[key]
            if tmpl[4] is None:
                self.__dict__[key]=val
            else:
                if tmpl[2] is None:
                    self.__dict__[key]=tmpl[4](val)
                else:
                    self.__dict__[key]=[ tmpl[4](v) for v in val ]

    def genDict(self):
        # Generate a regular dictionary
        ret = dict()
        for key in self.getItems():
            ret[key] = self.__dict__[key]
        return ret

    def genNestDict(self):
        # Generate a nest compatible dictionary.
        ret = dict()
        for key in self.getItems():
            if isinstance(self.__dict__[key],Distribution):
                ret[key] = self.__dict__[key].genNestDict()
            else:
                ret[key] = self.__dict__[key]
        return ret

class Distribution:
    _DistTypes = {
            'const'                           : ( 'val',                          ),
            'linear'                          : ( 'start', 'end'                  ),
            'normal'                          : ( 'mu', 'sigma'                   ),
            'normal_clipped'                  : ( 'mu', 'sigma', 'low ', 'high'   ),
            'normal_clipped_to_boundary'      : ( 'mu', 'sigma', 'low ', 'high'   ),
            'lognormal'                       : ( 'mu', 'sigma'                   ),
            'lognormal_clipped'               : ( 'mu', 'sigma', 'low', 'high'    ),
            'lognormal_clipped_to_boundary'   : ( 'mu', 'sigma', 'low', 'high'    ),
            'uniform'                         : ( 'low', 'high'                   ),
            'uniform_int'                     : ( 'low', 'high'                   ),
            'binomial'                        : ( 'n', 'p'                        ),
            'binomial_clipped'                : ( 'n', 'p', 'low', 'high'         ),
            'binomial_clipped_to_boundary'    : ( 'n', 'p', 'low', 'high'         ),
            'gsl_binomial'                    : ( 'n', 'p'                        ),
            'exponential'                     : ( 'lambda',                       ),
            'exponential_clipped'             : ( 'lambda', 'low', 'high'         ),
            'exponential_clipped_to_boundary' : ( 'lambda', 'low', 'high'         ),
            'gamma'                           : ( 'order', 'scale'                ),
            'gamma_clipped'                   : ( 'order', 'scale', 'low', 'high' ),
            'gamma_clipped_to_boundary'       : ( 'order', 'scale', 'low', 'high' ),
            'poisson'                         : ( 'lambda',                       ),
            'poisson_clipped'                 : ( 'lambda', 'low', 'high'         ),
            'poisson_clipped_to_boundary'     : ( 'lambda', 'low', 'high'         )
            }

    def __init__(self,args):
        # __init__(Distribution object)
        # __init__(distType,param1, param2,...):
        self.dist = None
        self.params = None
        self.data = None

        from copy import deepcopy

        if isinstance(args,Distribution):
            # copy constructor
            self.dist = deepcopy(args.dist)
            self.params = deepcopy(args.params)
            self.data = deepcopy(args.data)
        else:
            if not isinstance(args,(list,tuple)):
                raise Errors.ParamTypeError('args','Distribution, list or tuple')
            if len(args)<2:
                raise Errors.ParamSizeError('args','>= 2')
            dist = args[0]
            params = args[1:]

            if not isinstance(dist,str):
                raise Errors.ParamTypeError("dist","string")

            distTypes = Distribution._DistTypes

            if dist not in distTypes.keys():
                raise Errors.ParamError("dist","value", \
                    "Undefined distribution type {}".format(dist))
            self.dist = dist

            if len(distTypes[dist])!=len(params):
                raise Errors.ParamSizeError("params",len(distTypes[dist]))
            for v in params:
                if not isinstance(v,(int,float)):
                    raise Errors.ParamTypeError("params content","int or float")
            self.params = [ float(v) for v in params ]

    def genNestDict(self):
        return self.params[0]
        #if self.dist == 'const' or self.dist == 'linear':
            #return self.params[0]
        #else:
            #ret = {'distribution' : self.dist}
            #for i,p in enumerate(Distribution._DistTypes[self.dist]):
                #ret[p]=self.params[i]
            #return ret

    def rand(self,N):
        if self.dist == 'const':
            self.data = [ self.params[0] for i in range(N) ]
        if self.dist == 'linear':
            import numpy as np
            self.data = np.linspace(self.params[0],self.params[1],N).tolist()
        if self.dist == 'uniform':
            from numpy.random import uniform
            self.data = uniform(low=self.params[0],high=self.params[1],size=N).tolist()


def DataSegment(data,t0i=None,t1i=None):
    import numpy as np
    if not isinstance(data,np.ndarray):
        data = np.array(data)

    if t0i is None:
        t0i = 0
    if t1i is None:
        t1i = data.shape[-1]

    return np.take(data,range(t0i,t1i),data.ndim-1)

class Time:
    def __init__(self,seconds=0.0, msec=None, minutes=None, hours=None, days=None):
        if isinstance(seconds,str):
            #convert time string to Time object
            res=seconds
            timekeys = ('d','h','\'','"')
            timelist = [ 0.0 , 0.0 , 0.0  , 0.0 ]
            for i,k in enumerate(timekeys):
                if k in res:
                    splits = res.split(k)
                    if len(splits)!=2:
                        raise Errors.TimeStrError
                    keystr=splits[0].strip()
                    if not keystr:
                        raise Errors.TimeStrError
                    timelist[i] = float(keystr)
                    res=splits[1]
            days, hours, minutes, seconds = timelist
            if res.strip():
                msec = float(res)
        self.sec = float(seconds)
        if msec is not None:
            self.sec = self.sec + float(msec) / 1000.0
        if minutes is not None:
            self.sec = self.sec + float(minutes) * 60.0
        if hours is not None:
            self.sec = self.sec + float(hours) * 3600.0
        if days is not None:
            self.sec = self.sec + float(days) * 3600.0 * 24.0

    def as_msec(self):
        return self.sec * 1000.0

    def __float__(self):
        # return the msec value
        return self.sec

    def __add__(self,x):
        return Time(seconds = self.sec+float(x))

    def __sub__(self,x):
        return Time(seconds = self.sec-float(x))

    def __cmp__(self,x):
        return cmp(float(self),float(x))

    def __repr__(self):
        import numpy as np
        seconds = self.sec
        sign = "" if seconds>=0 else "-"
        seconds = abs(seconds)

        msec = int(seconds % 1.0 * 100)
        sec = int(seconds % 60)
        minutes = int(seconds / 60)

        timestr = "\""
        if msec > 0:
            timestr = ".{:02d}".format(msec)+timestr
        if minutes>0:
            timestr = "{:02d}".format(sec)+timestr
        else:
            timestr = "{:d}".format(sec)+timestr
        if minutes > 0:
            timestr = "{:d}'".format(minutes)+timestr
        timestr = sign+timestr
        return timestr

class fgcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'
