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

def DrawOld ( know, novel, figname ) :

    numbins = 100
    fig     = plt.figure ()
    plt.subplot(2,1,1)
    if len(know) : plt.hist( know , numbins, normed = 1, facecolor = 'b', label = 'Know' )
    plt.legend(loc='upper left')
    plt.xlim (-25, 10)
    plt.ylabel('Number')

    plt.subplot(2,1,2)
    if len(novel): plt.hist( novel, numbins, normed = 1, facecolor = 'r', label = 'Novel')
    plt.legend(loc='upper left')
    plt.xlim (-25, 10)
    plt.xlabel('VQ Distribution')
    plt.ylabel('Number')

    fig.savefig( figname + '.png' )
def Draw ( newData, trainData, figname ) :

    for k in newData.keys() :

        newData[k]   = np.array( newData[k]   )
        trainData[k] = np.array( trainData[k] )

        fig = plt.figure ()
        plt.subplot(2,1,1)
        if len( trainData[k] ) : plt.scatter(trainData[k][:,0], trainData[k][:,1], c='r', marker='.', label = 'Training Data')
        plt.legend(loc='upper left')
        plt.ylabel( 'Evidence for ' + k )

        plt.subplot(2,1,2)
        if len(newData[k] ) : plt.scatter(newData[k][:,0], newData[k][:,1], c='g', marker='.', label = 'New Data')
        plt.legend(loc='upper left')
        plt.xlabel( 'VQ' )
        plt.ylabel( 'Evidence for ' + k )

        fig.savefig( figname + '.' + k + '.png' )
    

def main ( opt ) :

    vcfInfile = opt[0]
    if vcfInfile[-3:] == '.gz' :
        I = os.popen( 'gzip -dc %s' % vcfInfile )
    else :
        I = open ( opt.vcfInfile )

    newdata,trainData = {}, {}
    while 1 :

        lines = I.readlines( 100000 )
        if not lines : break

        for line in lines :
            if re.search(r'^#', line) : continue
            col= line.strip('\n').split()

            vcfinfo = { d.split('=')[0] : d.split('=')[1] for d in col[7].split(';') if len(d.split('=')) == 2 }
            vq      = string.atoi( vcfinfo['VQ'] )
            culprit = string.atof( vcfinfo[ vcfinfo['CU'] ] )

            if vcfinfo['CU'] not in newdata   : newdata[vcfinfo['CU']]   = []
            if vcfinfo['CU'] not in trainData : trainData[vcfinfo['CU']] = []

            if re.search(r'_TRAIN_SITE', col[7]) : 
                trainData[vcfinfo['CU']].append( [vq, culprit] )
            else :
                newdata[vcfinfo['CU']].append( [vq, culprit] )

    I.close()

    figName = 'fig'
    if len(opt) > 1 : figName = opt[1]
    Draw ( newdata, trainData, figName )

if __name__ == '__main__' :

    usage = '\nUsage : %prog '
    main ( sys.argv[1:])

