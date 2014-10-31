"""
==============================
Draw ROC Line
==============================
Author: Shujia Huang
Date  : 2014-04-22 19:12:54
"""
import sys
import re
import os
import optparse
import string
import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics import roc_curve, auc


def main(opt):

    if opt.vcfInfile[-3:] == '.gz':
        I = os.popen('gzip -dc %s' % opt.vcfInfile)
    else:
        I = open (opt.vcfInfile)
    data = []
    while 1:
        lines = I.readlines(100000)
        if not lines: break

        for line in lines:

            if re.search(r'^#', line): continue
            col  = line.strip('\n').split()
            """
            flag     = re.search(r';Novel=\(([^\|]+)\|', col[7]).group(1)
            if flag != 'Novel' and flag != 'ESame' and flag != 'Same': continue
            vq = string.atof(re.search(r';VQ=([^;]+)', col[7]).group(1))
            v  = 0
            if 'Same' in flag: v = 1
            """
            if not re.search(r'_TRAIN_SITE', col[7]): continue
            vq = string.atof(re.search(r';VQ=([^;]+)', col[7]).group(1))
            vq = int(vq + 0.5)
            v  = 1
            if 'NEGATIVE_TRAIN_SITE' in col[7]: v = 0
            data.append([vq, v])
    I.close()
    DrawFig(opt.figure, np.array(data))

def DrawFig(figPrefix, data):

    fpr, tpr, thresholds = roc_curve(data[:, 1], data[:, 0])
    roc_auc = auc(fpr, tpr)

    fig = plt.figure()
    plt.title('ROC', fontsize = 14)
    plt.plot(fpr, tpr, 'ro-', label='ROC curve (area = %0.2f)' % roc_auc)
    plt.plot([0, 1], [0, 1], 'k--')
    plt.legend(loc="lower right")
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')

    #plt.show()
    fig.savefig(figPrefix + '.png')
    fig.savefig(figPrefix + '.pdf')

    for i in range(len(thresholds)):
        print '%f\t%f\t%f' % (thresholds[i], fpr[i], tpr[i])


if __name__ == '__main__':

    usage = '\nUsage: %prog [-i vcfInfile] [-f figPrefix] > Output'
    optp  = optparse.OptionParser(usage=usage)
    optp.add_option('-i', '--vcf' , dest='vcfInfile', metavar='VCF', help='Input vcf file.'      , default=[])
    optp.add_option('-f', '--fig' , dest='figure'   , metavar='FIG', help='Figure output prefix.', default='test')

    opt, _ = optp.parse_args()
    print >> sys.stderr, '# [INFO] Paraeters: python' , sys.argv[0], '\n\t-i', opt.vcfInfile, \
         '\n\t-f', opt.figure, '\n'
    main(opt)
    print >> sys.stderr, '*********************** ALL DONE ***********************\n'



