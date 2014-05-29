"""
===============================================
===============================================
Author : Shujia Huang
Date   : 2014-05-23 11:21:53
"""

import sys
import re
import os
import optparse
import string
# My own class
import model.VariantDataManager as vdm
import model.VariantRecalibratorEngine as vre
import model.VariantRecalibratorArgumentCollection as VRAC
import model.VariantRecalibrator as vror

def main ( opt ) :

    traningSet     = vdm.LoadTrainingSiteFromVCF( opt.trainData ) # Just record the sites of training data
    hInfo, dataSet = vdm.LoadDataSet( opt.vcfInfile, traningSet, vdm.LoadFaLen(opt.qFalen) )
    vr             = vror.VariantRecalibrator()
    vr.OnTraversalDone( dataSet )

    # Outputting the result as VCF format
    hInfo.Add ('##INFO=<ID=VQ', '##INFO=<ID=VQ,Number=1,Type=String,Description="Variant Quality">')
    hInfo.Add ('##INFO=<ID=culprit', '##INFO=<ID=culprit,Number=1,Type=String,Description="The annotation which was the worst performing in the Gaussian mixture model, likely the reason why the variant was filtered out">')
    for k,v in sorted (hInfo.header.items(), key = lambda d : d[0] ) : print v

    for d in dataSet :
        # Deal with the INFO line
        vcfinfo = {}
        for info in d.variantContext[7].split(';') : 
            k = info.split('=')[0]
            if k in vcfinfo: raise ValueError('[ERROR] The tag: %s double hits in the INFO column at %s'%(k, info))
            vcfinfo[k] = info
        vcfinfo['VQ']       = 'VQ=' + str( int(d.lod+0.5) )
        vcfinfo['culprit']  = 'culprit=' + d.annoTexts[d.worstAnnotation]
        d.variantContext[7] = ';'.join( sorted(vcfinfo.values()) )
        print '\t'.join( d.variantContext )

if __name__ == '__main__' :

    usage = '\nUsage : %prog [--Train Training data set] [-i SampleVcfInfile] > Output'
    optp  = optparse.OptionParser(usage=usage)
    optp.add_option('-i', '--InVcf' , dest='vcfInfile', metavar='VCF', help='VCF for predict.', default=[] )
    optp.add_option('-T', '--Train' , dest='trainData', metavar='TRU', help='Traing data set at true  category', default=[] )
    optp.add_option('-q', '--qFalen', dest='qFalen'   , metavar='LEN', help='The length list of query sequence', default=[] )

    opt, _ = optp.parse_args()
    if len( opt.vcfInfile ) == 0 : optp.error( "Required [-i vcfInfile]\n"             )
    if len( opt.trainData ) == 0 : optp.error( "Required [-T trainData. VCF Format]\n" )
    if len( opt.qFalen    ) == 0 : optp.error( "Required [-q Query fa length list ]\n" )
    print >> sys.stderr, '[INFO] Parameters: python' , sys.argv[0], '\n\t-i', opt.vcfInfile, \
          '\n\t-T', opt.trainData, '\n\t-q', opt.qFalen, '\n'

    main( opt )
    print >> sys.stderr, '*********************** ALL DONE ***********************\n'

