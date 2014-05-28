"""
==============================
Draw ROC Line
==============================
Author : Shujia Huang
Date   : 2014-04-22 19:12:54
"""
import sys
import re
import os
import optparse
import string
import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics import roc_curve, auc


def main ( opt ) :

    sampleID = opt.sample

    if opt.vcfInfile[-3:] == '.gz' :
        I = os.popen( 'gzip -dc %s' % opt.vcfInfile )
    else :
        I = open ( opt.vcfInfile )
    data = []
    while 1 :
        lines = I.readlines( 100000 )
        if not lines : break

        for line in lines :

            line = line.strip('\n')
            col  = line.split()
            if re.search(r'^#CHROM', line) :
                if len(col) < 10 :
                    print >> sys.stderr, '# [ERROR] The input vcf file (%s) does not contain the "Sample" fields!\n' % vcfInfile
                    sys.exit(1)
                sam2col = { sam:i+9 for i,sam in enumerate(col[9:]) }
                if sampleID not in sam2col :
                    print >> sys.stderr, '# [ERROR] The input sample id (%s) is not match in Vcf file(%s)\n' % ( sampleID, vcfInfile)
                    sys.exit(1)
            if re.search(r'^#', line) : continue
            idx  = sam2col[sampleID]
            
            flag = re.search(r';Novel=\(([^\|]+)\|', col[7]).group(1)
            if flag != 'Novel' and flag != 'ESame' and flag != 'Same' : continue

            fmat = { k:i for i,k in enumerate( col[8].split(':') ) } # Get Format
            for type in ['VQ'] :
                if type not in fmat :
                    print >> sys.stderr, '# [ERROR] The "Format" fields did not contian %s in VCF %s' %( 'AA', opt.vcfInfile )
                    sys.exit(1)

            vq = string.atoi(col[idx].split(':')[ fmat['VQ'] ])
            v  = 0
            if 'Same' in flag : v = 1
            data.append([vq, v])
                
    I.close()
    DrawFig( opt.figure, np.array(data) )

def DrawFig ( figPrefix, data ) :

    fpr, tpr, thresholds = roc_curve(data[:, 1], data[:, 0])
    roc_auc = auc(fpr, tpr)

    fig = plt.figure()
    plt.title('ROC', fontsize = 14)
    plt.plot(fpr, tpr, label='ROC curve (area = %0.2f)' % roc_auc)
    plt.plot([0, 1], [0, 1], 'k--')
    #plt.xlabel('Count of novel(not in HuRef dataset) events', fontsize = 14)
    #plt.ylabel('Fraction of known(in HuRef dataset) events', fontsize = 14)
    plt.legend(loc="lower right")
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')

    plt.show()
    fig.savefig(figPrefix + '.png')
    fig.savefig(figPrefix + '.pdf')

    for i in range( len(thresholds) ) :
        print '%f\t%f\t%f' % (thresholds[i], fpr[i], tpr[i])


if __name__ == '__main__' :

    usage = '\nUsage : %prog [-i vcfInfile] [-f figPrefix] > Output'
    optp  = optparse.OptionParser(usage=usage)
    optp.add_option('-i', '--vcf' , dest='vcfInfile', metavar='VCF', help='Input vcf file.'      , default=[] )
    optp.add_option('-f', '--fig' , dest='figure'   , metavar='FIG', help='Figure output prefix.', default='test' )
    optp.add_option('-s', '--smp' , dest='sample'   , metavar='SMP', help='Sample ID.'           , default=[] )

    opt, _ = optp.parse_args()
    print >> sys.stderr, '# [INFO] Paraeters: python' , sys.argv[0], '\n\t-i', opt.vcfInfile, \
          '\n\t-s', opt.sample, '\n\t-f', opt.figure, '\n'
    main( opt )
    print >> sys.stderr, '*********************** ALL DONE ***********************\n'



