"""
================================================
My own Gaussion Mixture Model for SV genotyping.
Learn form scikit-learn
================================================

Author : Shujia Huang
Date   : 2014-01-06 14:33:45

"""
import sys
import string
import re
import numpy as np
# My own class
import VariantRecalibratorArgumentCollection as VRAC
import VariantDatum as vd
import VCF

class VariantDataManager :

    def __init__ ( self, data=None ) :
        self.VRAC = VRAC.VariantRecalibratorArgumentCollection()
        self.annotationMean = None
        self.annotationSTD  = None

        self.data = [] # list < VariantDatum >
        if data : # data is not None
            if not isinstance(data[0],vd.VariantDatum): raise ValueError('[ERROR] The data type should be "VariantDatum" in VariantDataManager(),but found %s'% str(type(data[0])) )
            self.data = data
            for i, d in enumerate( self.data ) : self.data[i].annotations = np.array( self.data[i].annotations )

    def SetData (self, data) :
        if not data or not isinstance(data[0],vd.VariantDatum): raise ValueError ('[ERROR] The data type should be "VariantDatum" in VariantDataManager(),but found %s' % str(type(data[0])) )
        self.data = data
        for i, d in enumerate( self.data ) : self.data[i].annotations = np.array( self.data[i].annotations )

    def NormalizeData (self) :
        data = np.array( [ d.annotations for d in self.data ], dtype=float )
        mean = data.mean(axis=0); self.annotationMean = mean
        std  = data.std(axis=0) ; self.annotationSTD  = std

        # foundZeroVarianceAnnotation
        if any( std < 1e-5 ) : raise ValueError ( '[ERROR] Found annotations with zero variance. They must be excluded before proceeding.' )
        
        # Each data point is now (x - mean) / standard deviation
        for i, d in enumerate( data ) :
            self.data[i].annotations = (d - mean) / std 
            # trim data by standard deviation threshold and mark failing data for exclusion later
            remove = False
            if any( np.abs(self.data[i].annotations) > self.VRAC.STD_THRESHOLD ) : remove = True
            self.data[i].failingSTDThreshold = remove

    def GetTrainingData (self) :

        trainingData = [ d for d in self.data if (not d.failingSTDThreshold) and d.atTrainingSite ]
        print >> sys.stderr, '[INFO] Training with %d variants after standard deviation thresholding.' % len(trainingData)
        if len(trainingData) < self.VRAC.MIN_NUM_BAD_VARIANTS :
            print >> sys.stderr,'[WARNING] Training with very few variant sites! Please check the model reporting PDF to ensure the quality of the model is reliable.'
        if len(trainingData) > self.VRAC.MAX_NUM_TRAINING_DATA:
            print >> sys.stderr, '[WARING] Very large training set detected. Downsampling to %d training variants.' % self.VRAC.MAX_NUM_TRAINING_DATA
            np.random.shuffle(trainingData) # Random shuffling
            return trainingData[ range(self.VRAC.MAX_NUM_TRAINING_DATA) ]
        return trainingData 

    def SelectWorstVariants(self ) :

        trainingData = [ d for d in self.data if (d.lod < self.VRAC.BAD_LOD_CUTOFF) and (not d.failingSTDThreshold) ]
        print >> sys.stderr, '[INFO] "Training with worst %d scoring variants --> variants with LOD <= %.2f.' % ( len(trainingData), self.VRAC.BAD_LOD_CUTOFF )
        return trainingData


def SampleFaLen ( faLenFile ) :

    if faLenFile[-3:] == '.gz' : I = os.popen( 'gzip -dc %s' % faLenFile )
    else                       : I = open( faLenFile )
    data = {}
    while 1 :
        lines = I.readlines ( 100000 )
        if not lines : break
        for line in lines :
            col = line.strip('\n').split()
            data[col[0]] = string.atoi( col[1] )
    I.close()
    return data

def LoadFaLen ( faLenLstFile ) :

    data = {}
    I    = open (faLenLstFile)
    for line in I.readlines() :

        if len( line.strip('\n').split() ) != 2: raise ValueError('[ERROR] The format of Fa length list maybe not right. It could just be : "sample FalenghtFile", but found',line)
        sampleId, fileName = line.strip('\n').split()
        if sampleId not in data : data[sampleId] = {}
        data[sampleId] = SampleFaLen( fileName )
    I.close()
    return data

def LoadTrainingSiteFromVCF( vcffile ) :

    if vcffile[-3:] == '.gz' :
        I = os.popen( 'gzip -dc %s' % vcffile )
    else :
        I = open ( vcffile )

    dataSet = set()
    while 1 :

        lines = I.readlines( 100000 )
        if not lines : break

        for line in lines :

            if re.search(r'^#', line) : continue
            col  = line.strip('\n').split()
            if not re.search(r'^PASS', col[6] ) : continue
            dataSet.add( col[0] + ':' + col[1] )
    I.close()

    return dataSet

def LoadDataSet ( vcfInfile, traningSet, qFaLen ) :

    if len( traningSet ) == 0 : raise ValueError ('[ERROR] No Training Data found')
    if vcfInfile[-3:] == '.gz' :
        I = os.popen( 'gzip -dc %s' % vcfInfile )
    else :
        I = open ( vcfInfile )

    data, hInfo = [], VCF.VCFHeader()
    while 1 : # VCF format

        lines = I.readlines( 100000 )
        if not lines : break
        for line in lines :

            col = line.strip('\n').split()
            if re.search (r'^##fileformat', line ) : hInfo.Add( '###', line.strip('\n') )
            elif re.search( r'^#CHROM', line ) :
                col2sam = { i+9:sam for i,sam in enumerate(col[9:]) }
                hInfo.Add( '#CHROM', line.strip('\n') )
            elif re.search(r'^#', line) :
                tag = line.strip('\n').split(',')[0]
                hInfo.Add( tag , line.strip('\n') )
            if re.search( r'^#', line) : continue

            fmat = { k:i for i,k in enumerate( col[8].split(':') ) } # Get Format
            for tag in ['AA', 'QR', 'FN'] :
                if tag not in fmat : raise ValueError ('[ERROR] The "Format" fields did not contian %s in VCF %s' %( tag, opt.vcfInfile) )


            annotations = []
            for i, sample in enumerate ( col[9:] ) : 
                sampleId = col2sam[9+i]
                if len( sample.split(':')[fmat['AA']].split(',') ) != 4 :
                    print >> sys.stderr,'[WARNING] %s\n%s' % (line, sample.split(':')[fmat['AA']])
                    continue

                qr  = sample.split(':')[fmat['QR']].split(',')[-1]
                qId, qSta, qEnd = qr.split('-')
                qSta = string.atoi(qSta)
                qEnd = string.atoi(qEnd)

                if sampleId not in qFaLen           : raise ValueError ('[ERROR] The sample name $s(in vcf) is not in the name of Fa list.' % sampleId )
                if      qId not in qFaLen[sampleId] : raise ValueError ('[ERROR]', qId, 'is not been found in file', opt.qFalen, '\n' )
                qSta= int( qSta * 100 / qFaLen[sampleId][qId] + 0.5 )
                qEnd= int( qEnd * 100 / qFaLen[sampleId][qId] + 0.5 )
                if qSta > 100 or qEnd > 100 : raise ValueError ('[ERROR] Query size Overflow! sample : %s; scaffold : %s' % (sampleId, qId) )

                leg = min(qSta, 100 - qEnd)
                nn  = string.atof(sample.split(':')[fmat['FN']])
                n   = int(1000 * nn + 0.5) / 10.0
                alt = string.atoi( sample.split(':')[fmat['AA']].split(',')[1] ) # Alternate perfect
                bot = string.atoi( sample.split(':')[fmat['AA']].split(',')[3] ) # Both imperfect
                annotations.append( [leg, n , alt, bot] )

            datum                = vd.VariantDatum()
            datum.annotations    = np.median( annotations, axis = 0 )
            datum.variantContext = col
            key                  = col[0] + ':' + col[1]
            if key in traningSet : datum.atTrainingSite = True
            data.append( datum )

    I.close()
    return hInfo, np.array(data)

