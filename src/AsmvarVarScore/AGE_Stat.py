"""
========================================================
Statistic the SV Stat after AGE Process
========================================================
Author : Shujia Huang
Date   : 2014-03-07 0idx:54:15
"""
import sys
import re
import os
import string
import numpy as np
import matplotlib.pyplot as plt

def DrawFig( figureFile, distance, leftIden, rigthIden, nr, aa, bb ) : 

    fig = plt.figure( num=None, figsize=(16, 18), facecolor='w', edgecolor='k' )
    plt.subplot(321)
    plt.title('Distance distribution', fontsize=16)
    P = distance[:,0] == 1; N = distance[:,0] == 2; X = distance[:,0] == 3
    plt.scatter(distance[:,1][P] , distance[:,2][P] )
    plt.xlabel('The breakpoints of varints span on assemble sequence(%)', fontsize=16)
    plt.ylabel('% of Number', fontsize=16)

    fig.savefig('test.png')
    #fig.savefig(figureFile + '.pdf')

def Accum ( data, isBig = False) :

    tmpD= data
    k   = sorted( tmpD.keys(), key = lambda d: float(d) )
    dat = []
    for i in range( len(k) ) :
        if isBig : 
            for j in range(i,len(k)) : tmpD[k[i]][1] += tmpD[k[j]][0]
        else : 
            for j in range(i+1) : tmpD[k[i]][1] += tmpD[k[j]][0]
        dat.append( [ float(k[i]), float(tmpD[k[i]][0]), float(tmpD[k[i]][1]) ] )
    
    return dat

def SampleFaLen ( faLenFile ) :

    if faLenFile[-3:] == '.gz' : I = os.popen( 'gzip -dc %s' % faLenFile )
    else                       : I = open( faLenFile )
    data = {}
    while 1 :
        lines = I.readlines ( 100000 )
        if not lines : break
        for line in lines :
            col = line.strip('\n').split()
            data[col[0]] = string.atoi( col[1] )
    I.close()
    return data

def LoadFaLen ( faLenLstFile ) :

    data = {}
    I    = open (faLenLstFile)
    for line in I.readlines() :

        if len( line.strip('\n').split() ) != 2: raise ValueError('[ERROR] The format of Fa length list maybe not right. It could just be : "sample FalenghtFile", but found',line)
        sampleId, fileName = line.strip('\n').split()
        if sampleId not in data : data[sampleId] = {}
        data[sampleId] = SampleFaLen( fileName )
    I.close()
    return data

def main ( argv ) :

    qFaLen    = LoadFaLen( argv[1] )
    figPrefix = 'test'
    if len(argv) > 2 : figPrefix = argv[2]
    if argv[0][-3:] == '.gz' :
        I = os.popen( 'gzip -dc %s' % argv[0] )
    else :
        I = open ( argv[0] )

    s, annotations, mark = set(), [], []
    print '#Chr\tPosition\tDistance\tLeftIden\tRightIden\tAveIden\tN-Ratio\tAA'
    while 1 : # VCF format

        lines = I.readlines( 100000 )
        if not lines : break

        for line in lines :

            col = line.strip('\n').split()
            if re.search( r'^#CHROM', line ) : col2sam = { i+9:sam for i,sam in enumerate(col[9:]) }
            if re.search(r'^#', line) : continue
            key = col[0] + ':' + col[1]
            if key in s : continue
            s.add(key)

            #if re.search(r'^PASS', col[6] ) : continue
            #if not re.search(r'_TRAIN_SITE', col[7]) : continue
            if not re.search(r'^PASS', col[6] ) : continue

            fmat = { k:i for i,k in enumerate( col[8].split(':') ) }
            if 'VS' not in fmat or 'QR' not in fmat: continue
            if len(annotations) == 0 : annotations = [ [] for _ in col[9:] ]

            vcfinfo = { d.split('=')[0] : d.split('=')[1] for d in col[7].split(';') if len(d.split('=')) == 2 }
            vq      = string.atof( vcfinfo['VQ'] )
            if 'POSITIVE_TRAIN_SITE' in col[7] and 'NEGATIVE_TRAIN_SITE' in col[7] :
                mark.append( [3, vq] )
            elif 'POSITIVE_TRAIN_SITE' in col[7] : 
                mark.append( [1, vq] )
            elif 'NEGATIVE_TRAIN_SITE' in col[7] : 
                mark.append( [2, vq] )
            else :
                mark.append( [0, vq] )

            for i, sample in enumerate ( col[9:] ) :
                sampleId = col2sam[9+i]
                qr = sample.split(':')[fmat['QR']].split(',')[-1]
                if qr == '.' : 
                    annotations[i].append( [0, 0, 0, 0, 0, 0] )
                    continue

                qId, qSta, qEnd = qr.split('-')
                qSta = string.atoi(qSta)
                qEnd = string.atoi(qEnd)

                if sampleId not in qFaLen : raise ValueError ('[ERROR] The sample name $s(in vcf) is not in the name of Fa list.' % sampleId )
                if qId not in qFaLen[sampleId] : raise ValueError ('[ERROR]', qId, 'is not been found in file', opt.qFalen, '\n' )
                qSta= int( qSta * 100 / qFaLen[sampleId][qId] + 0.5 )
                qEnd= int( qEnd * 100 / qFaLen[sampleId][qId] + 0.5 )
                if qSta > 100 or qEnd > 100 : raise ValueError ('[ERROR] Query size Overflow! sample : %s; scaffold : %s' % (sampleId, qId) )
                leg = qSta
                if 100 - qEnd < qSta : leg = qEnd
                nn  = string.atof(sample.split(':')[fmat['FN']])
                n   = int(1000 * nn + 0.5) / 10.0
                alt = string.atoi( sample.split(':')[fmat['AA']].split(',')[1] ) # Alternate perfect
                bot = string.atoi( sample.split(':')[fmat['AA']].split(',')[3] ) # Both imperfect
                pro = string.atoi( sample.split(':')[fmat['RP']].split(',')[0] ) # Proper Pair
                ipr = string.atoi( sample.split(':')[fmat['RP']].split(',')[1] ) # ImProper Pair
                annotations[i].append( [leg, n, alt, bot, pro, ipr] )
    I.close()
    if len( mark ) != len( annotations[0] ) : raise ValueError ('[ERROR] The size is not match!')
    annotations = np.array( annotations );

    sampleNum = len( annotations )
    for i in range( sampleNum ) : 
        if np.sum(annotations[i]) == 0: continue
        mean = np.array( [ d for d in annotations[i] if np.sum(d) > 0 ] ).mean(axis=0)
        std  = np.array( [ d for d in annotations[i] if np.sum(d) > 0 ] ).std (axis=0)
        annotations[i] = np.array( np.round( (annotations[i] - mean)/std) )  # Normalization Per sample
        print >> sys.stderr, '# Sample NO.', i + 1,'\n', mean,'\n', std

    data, distance, leftIdn, rightIdn, nr, aa, bb = [],[],[],[],[],[],[]
    for i in range( len(annotations[0]) ) : 

        if not mark[i][0] : continue # Next if makr[i] == 0

        anno = np.array( [ annotations[s][i] for s in range( sampleNum ) if len(annotations[s][i][annotations[s][i]!=0]) > 0 ] ) # each person in the same position
        if len( anno ) == 0 : continue
        leg, n, alt, bot, leftIden, rightIden = np.median( anno, axis=0 )

        distance.append( [ mark[i][0], mark[i][1], leg       ] )
        leftIdn.append ( [ mark[i][0], mark[i][1], leftIden  ] )
        rightIdn.append( [ mark[i][0], mark[i][1], rightIden ] )
        nr.append      ( [ mark[i][0], mark[i][1], n   ] )
        aa.append      ( [ mark[i][0], mark[i][1], alt ] )
        bb.append      ( [ mark[i][0], mark[i][1], bot ] )
 
        data.append([leg, alt, leftIden, rightIden, n, bot])
        print mark[i][0], mark[i][1], '\t', leg, '\t', leftIden, '\t', rightIden,'\t', n, '\t', alt, '\t', bot

    data = np.array(data)
    print >> sys.stderr, '\nPosition\tALTernatePerfect\tLeftIdentity\tRightIdentity\tAveIden\tNRatio\tBothImperfect'
    print >> sys.stderr, 'Means: ', data.mean(axis=0), '\nstd  : ', data.std(axis=0), '\nMedian: ', np.median( data, axis=0 ) 
    print >> sys.stderr, '25 Percentile:', np.percentile(data, 25,axis=0), '\n50 Percentile:', np.percentile(data, 50,axis=0), '\n75 Percentile:', np.percentile(data, 75,axis=0)

    DrawFig( figPrefix, \
             np.array (distance ), \
             np.array (leftIdn  ), \
             np.array (rightIdn ), \
             np.array (nr       ), \
             np.array (aa       ), \
             np.array (bb       )  )

if __name__ == '__main__' :

    main(sys.argv[1:])

