_myConf = {}

def getConf(name,default=None):
    global _myConf
    return _myConf[name] if name in _myConf else default

def setConf(name,value=True):
    global _myConf
    _myConf[name] = value
