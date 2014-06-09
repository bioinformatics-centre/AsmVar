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
    hInfo, dataSet = vdm.LoadDataSet(opt.vcfInfile,traningSet,vdm.LoadFaLen(opt.qFalen)) # Identify the traning sites
    vr             = vror.VariantRecalibrator() # init VariantRecalibrator object
    vr.OnTraversalDone( dataSet ) # Traning model and calculate the VQ for all the dataSet

    # Outputting the result as VCF format
    hInfo.Add ('##INFO=<ID=VQ', '##INFO=<ID=VQ,Number=1,Type=Float,Description="Variant Quality">')
    hInfo.Add ('##INFO=<ID=CU', '##INFO=<ID=CU,Number=1,Type=String,Description="The annotation which was the worst performing in the Gaussian mixture model, likely the reason why the variant was filtered out. It\'s the same tag as "culprit" in GATK">')
    hInfo.Add ('##INFO=<ID=NEGATIVE_TRAIN_SITE', '##INFO=<ID=NEGATIVE_TRAIN_SITE,Number=0,Type=Flag,Description="This variant was used to build the negative training set of bad variants">')
    hInfo.Add ('##INFO=<ID=POSITIVE_TRAIN_SITE', '##INFO=<ID=POSITIVE_TRAIN_SITE,Number=0,Type=Flag,Description="This variant was used to build the positive training set of good variants">')
    # For Record the Annnotations' values
    hInfo.Add ('##INFO=<ID=Position', '##INFO=<ID=Position,Number=1,Type=Float,Description="The median of Position on Assembly Scaffold">')
    hInfo.Add ('##INFO=<ID=NRatio', '##INFO=<ID=NRatio,Number=1,Type=Float,Description="The median of N ratio of the query sequences">')
    hInfo.Add ('##INFO=<ID=AlternatePerfect', '##INFO=<ID=AlternatePerfect,Number=1,Type=Float,Description="The median of Depth of Alt_Perfect">')
    hInfo.Add ('##INFO=<ID=BothImperfect', '##INFO=<ID=BothImperfect,Number=1,Type=Float,Description="The median of Depth of Both_Imperfect">')
    
    culprit, tot = {}, 0.0
    for k,v in sorted (hInfo.header.items(), key = lambda d : d[0] ) : print v
    idx = {c:i for i,c in enumerate(['Position', 'NRatio', 'AlternatePerfect', 'BothImperfect']) }
    for d in dataSet :
        # Deal with the INFO line
        vcfinfo = {}
        for info in d.variantContext[7].split(';') : 
            k = info.split('=')[0]
            if k in vcfinfo: raise ValueError('[WARNING] The tag: %s double hits in the INFO column at %s'%(k, info))
            vcfinfo[k] = info

        tot += 1.0 # Record For summary
        culprit[d.annoTexts[d.worstAnnotation]] = culprit.get( d.annoTexts[d.worstAnnotation], 0.0 ) + 1.0 # For summary

        vcfinfo['VQ'] = 'VQ=' + str(d.lod)
        #vcfinfo['VQ'] = 'VQ=' + str( int(d.lod+0.5) )
        vcfinfo['CU'] = 'CU=' + d.annoTexts[d.worstAnnotation]
        vcfinfo['Position']         = 'Position=' + str(d.annotations[ idx['Position'] ])
        vcfinfo['NRatio'  ]         = 'NRatio=' + str(d.annotations[ idx['NRatio'] ])
        vcfinfo['AlternatePerfect'] = 'AlternatePerfect=' + str(d.annotations[ idx['AlternatePerfect'] ])
        vcfinfo['BothImperfect']    = 'BothImperfect=' + str(d.annotations[ idx['BothImperfect'] ])
        if d.atTrainingSite     : vcfinfo['POSITIVE_TRAIN_SITE'] = 'POSITIVE_TRAIN_SITE'
        if d.atAntiTrainingSite : vcfinfo['NEGATIVE_TRAIN_SITE'] = 'NEGATIVE_TRAIN_SITE'
        d.variantContext[7] = ';'.join( sorted(vcfinfo.values()) )

        print '\t'.join( d.variantContext )

    ## Output Summary
    print >> sys.stderr, '\n[Summmary] Here is the cupprit summary information:'
    for k, v in sorted( culprit.items(), key = lambda k:k[0] ) :
        print >> sys.stderr, '  ** %s: %.1f  %.2f%' % ( k, v, v * 100.0 / tot )

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

