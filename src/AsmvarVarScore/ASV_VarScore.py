"""
===============================================
===============================================
Author: Shujia Huang & Siyang Liu
Date  : 2014-05-23 11:21:53
"""

import sys
import re
import os
import optparse
import string
# My own class
import modul.VariantDataManager as vdm
import modul.VariantRecalibratorEngine as vre
import modul.VariantRecalibratorArgumentCollection as VRAC
import modul.VariantRecalibrator as vror

def main(opt):

    traningSet     = vdm.LoadTrainingSiteFromVCF(opt.trainData) # Just record the sites of training data
    hInfo, dataSet = vdm.LoadDataSet(opt.vcfInfile, traningSet, vdm.LoadFaLen(opt.qFalen)) #Identify the traning sites
    vr             = vror.VariantRecalibrator() # init VariantRecalibrator object
    vr.OnTraversalDone(dataSet) # Traning modul and calculate the VQ for all the dataSet
    vr.VisualizationLodVStrainingSet(opt.figure + '.BadLodSelectInTraining')

    # Outputting the result as VCF format
    hInfo.Add('INFO', 'VQ', 1, 'Float' , 'Variant Quality')
    hInfo.Add('INFO', 'CU', 1, 'String', 'The annotation which was the worst performing in the Gaussian mixture modul, likely the reason why the variant was filtered out. It\'s the same tag as <culprit> in GATK')
    hInfo.Add('INFO', 'NEGATIVE_TRAIN_SITE', 0, 'Flag', 'This variant was used to build the negative training set of bad variants')
    hInfo.Add('INFO', 'POSITIVE_TRAIN_SITE', 0, 'Flag', 'This variant was used to build the positive training set of good variants')
    # For Record the Annnotations' values
    for d in vr.dataManager.annoTexts: hInfo.Add('INFO', d[0], 1, d[1], d[2])
 
    culprit, good, tot = {}, {}, 0.0
    annoTexts = [d[0] for d in vr.dataManager.annoTexts]
    idx = {c:i for i,c in enumerate(annoTexts)}

    for k,v in sorted(hInfo.header.items(), key = lambda d: d[0]): print v
    if opt.vcfInfile[-3:] == '.gz':
        I = os.popen('gzip -dc %s' % opt.vcfInfile)
    else:
        I = open(opt.vcfInfile)

    j, monitor = 0, True
    while 1:

       lines = I.readlines(100000)
       if not lines: break

       for line in lines:

           if re.search(r'^#', line): continue

           col  = line.strip('\n').split()
           fmat = {k:i for i,k in enumerate(col[8].split(':'))} # Get Format
           if 'QR' not in fmat: continue # Cause by INTERGAP. But We'd better delete this statment, because the error is cause by the USER

           order = col[0] + ':' + col[1]
           d     = dataSet[j]
           j    += 1 # Increase the index of dataSet for the next cycle
           if d.variantOrder != order: 
               raise ValueError('[BUG] The order(%s) must be the same as dataSet(%s)' %(order, d.variantOrder))

           vcfinfo = {}
           # Deal with the INFO line
           for info in col[7].split(';'): 
               k = info.split('=')[0]
               if monitor and k in vcfinfo: 
                   monitor = False
                   print >>sys.stderr,'[WARNING] The tag: %s double hits in the INFO column at %s'%(k, opt.vcfInfile)
               vcfinfo[k] = info

           tot += 1.0 # Record For summary
           culprit[annoTexts[d.worstAnnotation]] = culprit.get(annoTexts[d.worstAnnotation], 0.0) + 1.0 #For summary
           d.lod = float('%.2f' % d.lod)
           for lod in [0, 1, 2, 3, 4]:
               if d.lod >= lod: good[lod] = good.get(lod, 0.0) + 1.0
        
           if d.atTrainingSite    : vcfinfo['POSITIVE_TRAIN_SITE'] = 'POSITIVE_TRAIN_SITE'
           if d.atAntiTrainingSite: vcfinfo['NEGATIVE_TRAIN_SITE'] = 'NEGATIVE_TRAIN_SITE'
           vcfinfo['VQ'] = 'VQ=' + str(d.lod)
           vcfinfo['CU'] = 'CU=' + annoTexts[d.worstAnnotation]
           for text in annoTexts: 
               vcfinfo[text] = text + '=' + str(d.annotations[idx[text]])
           col[7] = ';'.join(sorted(vcfinfo.values()))
           if d.lod < 0: d.lod = 0 # QUAL: donot allow value below 0
           col[5] = str(d.lod) # In fact QUAL field should use phred scala

           print '\t'.join(col)

    I.close()

    ## Output Summary
    print >> sys.stderr, '\n[Summmary] Here is the summary information:\n'
    for k, v in sorted(good.items(), key = lambda k:k[0]): 
        print >> sys.stderr, '  ** Variant Site score >= %d: %d\t%0.2f' % (k, v, v * 100 / tot)
    for k, v in sorted(culprit.items(), key = lambda k:k[0]):
        print >> sys.stderr, '  ** Culprit by %s: %d\t%.2f' % (k, v, v * 100.0 / tot)

if __name__ == '__main__':

    usage = '\nUsage: %prog [--Train Training data set] [-i SampleVcfInfile] > Output'
    optp  = optparse.OptionParser(usage=usage)
    optp.add_option('-i', '--InVcf' , dest='vcfInfile', metavar='VCF', help='VCF for predict.', default=[])
    optp.add_option('-T', '--Train' , dest='trainData', metavar='TRU', help='Traing data set at true  category', default=[])
    optp.add_option('-q', '--qFalen', dest='qFalen'   , metavar='LEN', help='The length list of query sequence', default=[])
    optp.add_option("-f", "--fig",    dest="figure",    metavar="FIG", help="The prefix of figure.",             default='figtest')

    opt, _ = optp.parse_args()
    if len(opt.vcfInfile) == 0: optp.error("Required [-i vcfInfile]\n"            )
    if len(opt.trainData) == 0: optp.error("Required [-T trainData. VCF Format]\n")
    if len(opt.qFalen   ) == 0: optp.error("Required [-q Query fa length list]\n")
    print >> sys.stderr, '[INFO] Parameters: python' , sys.argv[0], '\n\t-i', opt.vcfInfile, \
          '\n\t-T', opt.trainData, '\n\t-q', opt.qFalen, '\n\t-f', opt.figure, '\n'

    main(opt)
    print >> sys.stderr, '*********************** ALL DONE ***********************'

