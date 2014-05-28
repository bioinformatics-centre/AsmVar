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
    hInfo.Add ('##FORMAT=<ID=VQ', '##FORMAT=<ID=VQ,Number=1,Type=String,Description="Variant Quality">')
    for k,v in sorted (hInfo.header.items(), key = lambda d : d[0] ) : print v
    for d in dataSet :
        fmat = d.variantContext[8].split(':')
        if 'VQ' not in fmat :
            d.variantContext[8] += ':VQ' # Variant Quality
            for i,samp in enumerate( d.variantContext[9:] ) :
                d.variantContext[i+9] += ':' + str( int(d.lod+0.5) )
        else :
            for i,samp in enumerate( d.variantContext[9:] ) :
                fi             = samp.split(':')
                fi[fmat['VQ']] = str( int(d.lod+0.5) )
                d.variantContext[i+9] = ':'.join( fi )
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

