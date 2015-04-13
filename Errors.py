##    This module defines all error types

## The base class for all errors, 
#  takes only string msg as input
class Error(Exception):
    ## @brief Error constructor
    #  @param msg the error message string
    def __init__(self,msg=None):
        self.msg = msg

    ## @return an expression of message string
    def __str__(self):
        return repr(self.msg)

class ArgNumError(Error):
    def __init__(self,corrNum=None):
        if corrNum is None:
            self.msg = "Error number of arguments"
        else:
            self.msg = "This function should take {} arguments".format(corrNum)

class ParamError(Error):
    def __init__(self,pname,errtype,pcorr=None):
        if pcorr is None:
            self.msg = "Parameter {0} {1} error.".format(pname,errtype)
        else:
            self.msg = "Parameter {0} {1} error. {2}".format( pname,errtype,pcorr)

class ParamTypeError(ParamError):
    def __init__(self,pname,corrType):
        super(ParamTypeError,self).__init__(pname,'type',
                'Should be of type {}'.format(corrType))

class ParamDimError(ParamError):
    def __init__(self,pname,corrDim):
        super(ParamDimError,self).__init__(pname,'dimension',
                'Should be {}-D vector'.format(corrDim))

class ParamSizeError(ParamError):
    def __init__(self,pname,corrSize):
        if isinstance(corrSize,int):
            sizeStr = str(corrSize)
        else:
            sizeStr = ' x '.join(map(str,corrSize))
        super(ParamSizeError,self).__init__(pname,'size',
                'Should have size {}.'.format(sizeStr))


class ModelError(Error):
    def __init__(self,model,errmsg=None):
        if errmsg is None:
            self.msg = "Model {} error".format(model)
        else:
            self.msg = "Model {} error: {}".format(model,errmsg)

class ModelUndefinedError(ModelError):
    def __init__(self,model):
        super(ModelUndefinedError,self).__init__(model,"Undefined model")

class ModelParamTypeError(ModelError):
    def __init__(self,model):
        super(ModelParamTypeError,self).__init__(model,
                "model parameters should be a dictionary")

class ModelIncompleteError(ModelError):
    def __init__(self,model):
        super(ModelIncompleteError,self).__init__(model,
                'parameters have not initialized.')

class TemplateError(Error):
    def __init__(self,field,errmsg=None):
        if errmsg is None:
            self.msg = "Field {} error.".format(field)
        else:
            self.msg = "Field {} error: {}".format(field,errmsg)

class TemplUndefinedFieldError(TemplateError):
    def __init__(self,field):
        super(TemplUndefinedFieldError,self).__init__(field,
                'doesnot exist in the template.')

class TemplFieldTypeError(TemplateError):
    def __init__(self,field,corrtype=None):
        if corrtype is None:
            msg = "field type mismatch template."
        else:
            msg = "field type should be {}.".format(corrtype)
        super(TemplFieldTypeError,self).__init__(field,msg)

class TemplFieldLengthError(TemplateError):
    def __init__(self,field,length=None):
        if length is None:
            msg = "field length mismatch template."
        else:
            msg = "field length should be {}".format(length)
        super(TemplFieldLengthError,self).__init__(field,msg)

class TemplFieldContentError(TemplateError):
    def __init__(self,field):
        super(TemplFieldContentError,self).__init__(field,
                "content type mismatch template")

class NetError(Error):
    def __init__(self,errtype,errmsg=None):
        if errmsg is None:
            msg = "WSN net {} error.".format(errtype)
        else:
            msg = "WSN net {} error: {}.".format(errtype,errmsg)
        super(NetError,self).__init__(msg)

class NetParamError(NetError):
    def __init__(self,errmsg=None):
        super(NetParamError,self).__init__('parameter',errmsg)

class NetNotInitError(Error):
    def __init__(self):
        super(NetNotInitError,self).__init__('Network hasnot been initialized') 

class NetPathInexistError(NetError):
    def __init__(self):
        super(NetPathInexistError,self).__init__('SaveFile inexist')

class FileIOError(Error):
    def __init__(self,errmsg=None):
        if errmsg is None:
            msg = "IO error"
        else:
            msg = "IO Error: {}.".format(errmsg)
        super(FileIOError,self).__init__(msg)

class FileIOFileNotFound(FileIOError):
    def __init__(self,filename):
        errmsg = "File {} does not exist".format(filename)
        super(FileIOFileNotFound,self).__init__(errmsg)

class FileIORecordNotFound(FileIOError):
    def __init__(self,rec_name):
        errmsg = "Record {} does not exist".format(rec_name)
        super(FileIORecordNotFound,self).__init__(errmsg)

class TimeError(Error):
    def __init__(self,errmsg=None):
        if errmsg is None:
            msg = "Time object error"
        else:
            msg = "Time object error: {}".format(errmsg)
            super(TimeError,self).__init__(msg)

class TimeStrError(TimeError):
    def __init__(self):
        errmsg = "Corrupted time string format"
        super(TimeStrError,self).__init__(errmsg)

class TimeParaError(TimeError):
    def __init__(self,par):
        errmsg = "Incorrect parameter {}".format(par)
        super(TimeParaError,self).__init__(errmsg)
