"""
===================================================
A class for VCF format 
===================================================
"""

class VCFHeader :

    def __init__ ( self, hInfo = None ) : 

        self.header = {}
        if hInfo and (type( hInfo ) is not dict ) : raise ValueError ('The data type should be "dict" in class of "VCFHeader", but found %s' % str( type(hInfo) ))
        if hInfo : self.header = hInfo
        
    def Add ( self, key, context ) :
        self.header[key] = context
        return self

class VCFInfo : 

    def __init__ (self, info= None ) :
        self.info = {}
        if info and ( type(info) is not dict ) : raise ValueError ('The data type should be "dict" in class of "VCFInfo", but found %s' % str( type(info) ))
        if info : self.info = info

    def Add ( self, key, context ) :
        self.info[key] = context
        return self

class VCFContext : 

    def __init__ (self) :

        chrom = None
        pos   = None  
        Id    = None
        ref   = None  
        alt   = None 
        qual  = None    
        filters= None  
        info   = None   
        formats = None 
        sample  = []


        

