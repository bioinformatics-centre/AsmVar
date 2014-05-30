"""
==============================
==============================
Author : Shujia Huang
Date   : 2014-05-30 17:29:39
"""
import os
import sys
import re
import string
import numpy as np
import matplotlib.pyplot as plt

def Draw ( know, novel, figname ) :

    numbins = 100
    fig     = plt.figure ()
    if len(know) : plt.hist( know , numbins, normed = 1, facecolor = 'b', label = 'Know' )
    if len(novel): plt.hist( novel, numbins, normed = 1, facecolor = 'r', label = 'Novel')
    plt.legend(loc='upper left')
    plt.xlabel('VQ Distribution')
    plt.ylabel('Number')
    fig.savefig( figname + '.png' )

def main ( opt ) :

    vcfInfile = opt[0]
    if vcfInfile[-3:] == '.gz' :
        I = os.popen( 'gzip -dc %s' % vcfInfile )
    else :
        I = open ( opt.vcfInfile )

    know, novel = [], []
    while 1 :

        lines = I.readlines( 100000 )
        if not lines : break

        for line in lines :
            if re.search(r'^#', line) : continue
            col= line.strip('\n').split()
            vq = string.atoi( re.search(r';VQ=([^;]+)', col[7]).group(1) )
            if col[2][:2] == 'rs' : 
                know.append(vq)
                print '1\t',vq
            else : 
                novel.append(vq)
                print '0\t',vq

    I.close()

    figName = 'fig'
    if len(opt) > 1 : figName = opt[1]
    Draw ( np.array(know), np.array(novel), figName )

if __name__ == '__main__' :

    usage = '\nUsage : %prog '
    main ( sys.argv[1:])

