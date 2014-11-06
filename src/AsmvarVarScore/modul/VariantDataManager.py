"""
================================================
My own Gaussion Mixture Model for SV genotyping.
Learn form scikit-learn
================================================

Author: Shujia Huang
Date  : 2014-01-06 14:33:45

"""
import sys
import string
import re
import os
import numpy as np
from sklearn.metrics import roc_curve
# My own class
import VariantRecalibratorArgumentCollection as VRAC
import VariantDatum as vd
import VCF

class VariantDataManager:

    def __init__(self, data=None):
        self.VRAC = VRAC.VariantRecalibratorArgumentCollection()
        self.annotationMean = None
        self.annotationSTD  = None
        self.annoTexts      = [['AllelicNum', 'Float', 'Allelic number: bi- or mult-allelic'], \
                               ['InbCoeff',   'Float', 'Inbreeding coefficient: 1.0 - hetCount/Expected_hetCount'], \
                               ['Position',   'Float', 'Relative position on Assembly Scaffold'], \
                               ['NRatio'  ,   'Float', 'N normal ratio of the query sequences'], \
                               ['AlternatePerfect', 'Float', 'Support Alt_Perfect'  ], \
                               ['BothImperfect',    'Float', 'Support Both_Imperfect']]

        self.data = [] # list < VariantDatum >
        if data: # data is not None
            if not isinstance(data[0],vd.VariantDatum): 
                raise ValueError('[ERROR] The data type should be "VariantDatum" in VariantDataManager(),but found %s'% str(type(data[0])))
            self.data = data
            for i, d in enumerate(self.data): self.data[i].annotations = np.array(self.data[i].annotations)

    def SetData(self, data):

        if not isinstance(data[0],vd.VariantDatum): 
            raise ValueError('[ERROR] The data type should be "VariantDatum" in VariantDataManager(),but found %s' % str(type(data[0])))
        self.data = data
        for i, d in enumerate(self.data): self.data[i].annotations = np.array(d.annotations)

    def NormalizeData(self):

        data = np.array([d.annotations for d in self.data], dtype=float)
        mean = data.mean(axis=0); self.annotationMean = mean
        std  = data.std(axis=0) ; self.annotationSTD  = std

        # foundZeroVarianceAnnotation
        if any(std < 1e-5): raise ValueError('[ERROR] Found annotations with zero variance. They must be excluded before proceeding.')
        
        # Each data point now is (x - mean)/sd
        for i, d in enumerate(data):
            self.data[i].annotations = (d - mean) / std 
            # trim data by standard deviation threshold and mark failing data for exclusion later
            self.data[i].failingSTDThreshold = False
            if any(np.abs(self.data[i].annotations) > self.VRAC.STD_THRESHOLD):
                self.data[i].failingSTDThreshold = True

    def GetTrainingData(self):

        trainingData = [d for d in self.data if(not d.failingSTDThreshold) and d.atTrainingSite]
        print >> sys.stderr, '[INFO] Training with %d variants after standard deviation thresholding.' % len(trainingData)
        if len(trainingData) < self.VRAC.MIN_NUM_BAD_VARIANTS:
            print >> sys.stderr,'[WARNING] Training with very few variant sites! Please check the model reporting PDF to ensure the quality of the model is reliable.'
        if len(trainingData) > self.VRAC.MAX_NUM_TRAINING_DATA:
            print >> sys.stderr, '[WARING] Very large training set detected. Downsampling to %d training variants.' % self.VRAC.MAX_NUM_TRAINING_DATA
            np.random.shuffle(trainingData) # Random shuffling
            return [trainingData[i] for i in range(self.VRAC.MAX_NUM_TRAINING_DATA)]
        return trainingData 

    def SelectWorstVariants(self, badLod):

        trainingData = []
        for i,d in enumerate(self.data):
            if(d.lod < badLod) and(not d.failingSTDThreshold):
                trainingData.append(d)
                self.data[i].atAntiTrainingSite = True # I need the i order must be the same size as self.data 
        print >> sys.stderr, '[INFO] Training with worst %d scoring variants --> variants with LOD < %.2f.' %(len(trainingData), badLod)
        return trainingData

    def CalculateWorstLodCutoff(self):

        lodThreshold, lodCum = None, []
        if len(self.data) > 0:

            lodDist = np.array([[d.atTrainingSite, d.lod] for d in self.data if(not d.failingSTDThreshold)])

            # I just use the 'roc_curve' function to calculate the worst LOD threshold, not use it to draw ROC curve
            # And 'roc_curve' function will output the increse order, so that I don't have to sort it again
            _, tpr, thresholds = roc_curve(lodDist[:,0], lodDist[:,1]) 
            lodCum = [[thresholds[i], 1.0 - r] for i, r in enumerate(tpr)]

            for i, r in enumerate(tpr):
                if r > 1.0 - self.VRAC.POSITIVE_TO_NEGATIVE_RATE: 
                    lodThreshold = round(thresholds[i])
                    break

        return lodThreshold, np.array(lodCum)

def SampleFaLen(faLenFile):

    if faLenFile[-3:] == '.gz': I = os.popen('gzip -dc %s' % faLenFile)
    else                      : I = open(faLenFile)
    data = {}
    while 1:
        lines = I.readlines(100000)
        if not lines: break
        for line in lines:
            col = line.strip('\n').split()
            data[col[0]] = string.atoi(col[1])
    I.close()
    return data

def LoadFaLen(faLenLstFile):

    data = {}
    I    = open(faLenLstFile)
    for line in I.readlines():

        if len(line.strip('\n').split()) != 2: 
            raise ValueError('[ERROR] The format of Fa length list maybe not right. It could just be: "sample FalenghtFile", but found',line)
        sampleId, fileName = line.strip('\n').split()
        if sampleId not in data: data[sampleId] = {}
        data[sampleId] = SampleFaLen(fileName)
    I.close()
    return data

def LoadTrainingSiteFromVCF(vcffile):
    """
    Just record the training site positions
    """

    if vcffile[-3:] == '.gz':
        I = os.popen('gzip -dc %s' % vcffile)
    else:
        I = open(vcffile)

    dataSet = set()
    while 1:

        lines = I.readlines(100000)
        if not lines: break

        for line in lines:

            if re.search(r'^#', line): continue
            col = line.strip('\n').split()
            dataSet.add(col[0] + ':' + col[1])
    I.close()

    return dataSet

def LoadDataSet(vcfInfile, traningSet, qFaLen):

    if len(traningSet) == 0: raise ValueError('[ERROR] No Training Data found')
    if vcfInfile[-3:] == '.gz':
        I = os.popen('gzip -dc %s' % vcfInfile)
    else:
        I = open(vcfInfile)

    data, hInfo = [], VCF.VCFHeader()
    while 1: # VCF format

        lines = I.readlines(100000)
        if not lines: break
        for line in lines:

            col = line.strip('\n').split()
            if re.search(r'^#CHROM', line): col2sam = {i+9:sam for i,sam in enumerate(col[9:])}

            # Record the header information
            if re.search(r'^#', line):
                hInfo.Record(line.strip('\n'))
                continue

            fmat = {k:i for i,k in enumerate(col[8].split(':'))} # Get Format
            if 'QR' not in fmat: continue # Cause by INTERGAP. But We'd better delete this statment, because the error is cause by the USER 
            for tag in ['AA', 'QR', 'NR']:
                if tag not in fmat: raise ValueError('[ERROR] The "Format" fields did not contian "%s" in VCF: %s\nAT: %s\n' %(tag, vcfInfile, line))

            isBiallelic = True
            if len(col[4].split(',')) > 1: isBiallelic = False
            # Get inbreeding coefficient
            # It's calculated like: 1.0 - hetCount/Expected_hetCount in VCF
            #inbCoeff = re.search(r';?InbCoeff=([^;]+)', col[7])
            inbCoeff = re.search(r';F=([^;]+)', col[7])
            if not inbCoeff:
                print >> sys.stderr, '[ERROR] No inbreeding coefficient "InbCoeff=..." in INFO field in vcf:\n%s\n' % vcfInfile
            inbCoeff = float('%.2f' % float(inbCoeff.group(1)))

            annotations = []
            atleastOne  = False
            for i, sample in enumerate(col[9:]): 

                sampleId  = col2sam[9+i]
                if sample == './.': continue
                field = sample.split(':')
                if len(field[fmat['AA']].split(',')) != 4: continue

                if len(field) < fmat['QR'] + 1: continue
                qr    = field[fmat['QR']].split(',')[-1]
                if qr == '.': continue

                atleastOne = True
                qregion    = np.array(qr.split('-'))
                if len(qregion) > 3: qId = qregion[0] + '-' + qregion[1]
                else               : qId = qregion[0]
                qSta = string.atoi(qregion[-2])
                qEnd = string.atoi(qregion[-1])

                if sampleId not in qFaLen          : raise ValueError('[ERROR] The sample name $s(in vcf) is not in the name of Fa list.' % sampleId)
                if      qId not in qFaLen[sampleId]: raise ValueError('[ERROR]', qId, 'is not been found in fa file\n')
                qSta = int(qSta * 100 / qFaLen[sampleId][qId] + 0.5)
                qEnd = int(qEnd * 100 / qFaLen[sampleId][qId] + 0.5)
                if qSta > 100 or qEnd > 100: raise ValueError('[ERROR] Query size Overflow! sample: %s; scaffold: %s' %(sampleId, qId))

                leg = min(qSta, 100 - qEnd)
                nn  = string.atof(sample.split(':')[fmat['NR']])
                n   = int(1000 * nn + 0.5) / 10.0 # n ratio range: [0, 100]
                alt = string.atoi(sample.split(':')[fmat['AA']].split(',')[1]) # Alternate perfect
                bot = string.atoi(sample.split(':')[fmat['AA']].split(',')[3]) # Both imperfect
                annotations.append([isBiallelic, inbCoeff, leg, n , alt, bot])

            if not atleastOne: raise ValueError('[ERROR] All the samples don\'t contain this variant.', col)
            datum                = vd.VariantDatum()
            datum.annotations    = np.median(annotations, axis = 0)
            pos                  = col[0] + ':' + col[1]
            datum.variantOrder   = pos
            if pos in traningSet: datum.atTrainingSite = True
            data.append(datum)

    I.close()

    return hInfo, np.array(data)

