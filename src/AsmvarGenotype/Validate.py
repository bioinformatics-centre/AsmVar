"""
========================================
========================================
"""
import os
import sys
import re
import string
import numpy as np
import matplotlib.pyplot as plt

def DrawFig ( alleleCount ) :

    fig = plt.figure()
    if len( alleleCount ) > 0 :
        plt.title('AC Number', fontsize=12)
        plt.hist(alleleCount , 60, histtype='bar', normed=1, facecolor = 'c', color=['c'] )
        plt.ylabel('#', fontsize=12)

    figureFile = 'AC'
    fig.savefig(figureFile + '.png')
    fig.savefig(figureFile + '.pdf')

def main ( argv ) :

    vcfInfile = argv[0]
    if vcfInfile[-3:] == '.gz' :
        I = os.popen( 'gzip -dc %s' % vcfInfile )
    else :
        I = open ( vcfInfile )

    alleleCount = []
    while 1 :
        lines = I.readlines(100000)
        if not lines : break;

        for line in lines :
        # S
            line = line.strip( '\n' )
            if re.search(r'^#', line) : continue
            col  = line.split()
            if not re.search( r'^PASS', col[6]) : continue

            format = {}
            for i,t in enumerate( col[8].split(':') ) : format[t] = i # Get Format

            for type in ['RR', 'VT', 'VS'] :
                if type not in format.keys() :
                    print >> sys.stderr, '[ERROR] The format of VCF file is not right which you input, it did not contian %s field' % type
                    sys.exit(1)

            svsize = re.search(r';SVSIZE=([^;]+)', col[7])
            svtype = re.search(r';SVTYPE=([^;]+)', col[7])
            factor = re.search(r';F=([^;]+)', col[7])
            factor = string.atof( factor.group(1) )

            ac     = 0
            for sample in col[9:] : 
                gt = sample.split(':')[0]
                gnt= gt.split('/')
                if len(gnt) < 2 : gnt = gt.split('|')
                ac += ( string.atoi(gnt[0]) + string.atoi(gnt[1]) )
            print col[0],'\t',col[1], '\t', col[6], '\t', ac
            alleleCount.append(ac)


    I.close()
    DrawFig ( alleleCount )


if __name__ == '__main__' :

    main(sys.argv[1:])

