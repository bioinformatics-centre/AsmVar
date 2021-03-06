"""
========================================================
Statistic the SV Stat after AGE Process
========================================================
Author : Shujia Huang & Siyang Liu
Date   : 2014-03-07 0idx:54:15
"""
import sys
import re
import os
import string
import numpy as np
import matplotlib.pyplot as plt

def DrawFig( figureFile, distance, leftIden, rigthIden, aveIden, nr, aa, bb, test ) : 

    fig = plt.figure( num=None, figsize=(16, 18), facecolor='w', edgecolor='k' )
    plt.subplot(321)

    """
    from matplotlib.colors import LogNorm
    plt.hist2d(test[:,4], test[:,5], bins=50, norm=LogNorm())
    plt.plot(test[:,0], test[:,1], 'co')
    """
    plt.title('Distance distribution', fontsize=16)
    plt.plot(distance[:,0] , 100 * distance[:,1]/np.sum(distance[:,1])  , 'ro-' )
    plt.xlabel('The breakpoints of varints span on assemble sequence(%)', fontsize=16)
    plt.ylabel('% of Number', fontsize=16)

    plt.subplot(322)
    plt.title('Left Side', fontsize=16)
    plt.plot(leftIden[:,0] , leftIden[:,2]/np.sum(leftIden[:,1])  , 'go-' )
    plt.axis([0,100,0.0,1.0])
    plt.xlabel('Left Side Identity of varints(<=%)', fontsize=16)
    plt.ylabel('% of Accumulate', fontsize=16)

    plt.subplot(323)
    plt.title('Right Side', fontsize=16)
    plt.plot(rigthIden[:,0], rigthIden[:,2]/np.sum(rigthIden[:,1]), 'bo-' )
    plt.axis([0,100,0.0,1.0])
    plt.xlabel('Right Side Identity of varints(<=%)', fontsize=16)
    plt.ylabel('% of Accumulate', fontsize=16)

    plt.subplot(324)
    plt.title('Averge', fontsize=16)
    plt.plot(aveIden[:,0]  , aveIden[:,2]/np.sum(aveIden[:,1])    , 'co-' )
    plt.axis([0,100,0.0,1.0])
    plt.xlabel('Averge Identity of varints(<=%)', fontsize=16)
    plt.ylabel('% of Accumulate', fontsize=16)

    plt.subplot(325)
    plt.title('N Ratio', fontsize=16)
    plt.plot(nr[:,0], nr[:,2]/np.sum(nr[:,1]), 'yo-' )
    plt.axis([0,5,0.0,1.0])
    plt.xlabel('N Ratio of varints\' regions(>=%)', fontsize=16)
    plt.ylabel('% of Accumulate', fontsize=16)

    plt.subplot(6,2,10)
    plt.plot(aa[:,0], aa[:,2]/np.sum(aa[:,1]), 'mo-' )
    plt.axis([0,100,0.0,1.0])
    plt.xlabel('Perfect Depth(<=)', fontsize=12)
    plt.ylabel('% of Accumulate', fontsize=16)

    plt.subplot(6,2,12)
    plt.plot(bb[:,0], bb[:,2]/np.sum(bb[:,1]), 'ko-' )
    plt.axis([0,100,0.0,1.0])
    plt.xlabel('Both ImPerfect Depth(<=)', fontsize=12)
    plt.ylabel('% of Accumulate', fontsize=16)

    fig.savefig(figureFile + '.png')
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

def LoadFaLen ( faLenFile ) : 

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


def main ( argv ) :

    qFaLen    = LoadFaLen( argv[1] )
    figPrefix = 'test'
    if len(argv) > 2 : figPrefix = argv[2]
    if argv[0][-3:] == '.gz' :
        I = os.popen( 'gzip -dc %s' % argv[0] )
    else :
        I = open ( argv[0] )

    data, distance, leftIdn, rightIdn, aveIdn, nr, aa, bb = [],{},{},{},{},{},{},{}
    s    = set()
    print '#Chr\tPosition\tDistance\tLeftIden\tRightIden\tAveIden\tN-Ratio\tAA'
    while 1 : # VCF format

        lines = I.readlines( 100000 )
        if not lines : break

        for line in lines :

            if re.search(r'^#', line) : continue
            col = line.strip('\n').split()
            key = col[0] + ':' + col[1]
            if key in s : continue
            s.add(key)

            #if re.search(r'^PASS', col[6] ) : continue
            if not re.search(r'^PASS', col[6] ) : continue
            #if svsize < 500 : continue
            mage = re.search(r';?AGE=([^;]+)'  , col[7])
            if mage :
                age = mage.group(1).split(':')
            else    :
                age = ['0','0','0','0','0','0','0']

            fmat = { k:i for i,k in enumerate( col[8].split(':') ) }
            if 'VS' not in fmat or 'QR' not in fmat: continue

        
            idx = 9 
            if col[idx].split(':')[fmat['QR']] == '.': continue
            for type in ['AA', 'QR', 'VS'] :
                if type not in fmat :
                    print >> sys.stderr, '# [ERROR] The "Format" fields did not contian %s in VCF %s' %( 'AA', argv[0] )
                    sys.exit(1)
            svsize = string.atoi( col[idx].split(':')[fmat['VS']] )
            qId    = col[idx].split(':')[fmat['QR']].split(',')[-1].split('-')[0]
            qSta   = string.atof(col[idx].split(':')[fmat['QR']].split(',')[-1].split('-')[1])
            qEnd   = string.atof(col[idx].split(':')[fmat['QR']].split(',')[-1].split('-')[2])

            if qId not in qFaLen : 
                print >> sys.stderr, '#[ERROR]', qId, 'is not been found in file', argv[1], '\n'
                sys.exit(1)
            qSta= int( qSta * 100 / qFaLen[qId] + 0.5 )
            qEnd= int( qEnd * 100 / qFaLen[qId] + 0.5 )
            if qSta > 100 or qEnd > 100 : 
                print >> sys.stderr, '#[ERROR] Query size Overflow!', col[7]
                sys.exit(1)
            leg = min(qSta, 100 - qEnd)
            #if 100 - qEnd < qSta : leg = qEnd

            leftIden, rightIden = string.atoi(age[2]), string.atoi(age[4])
            avelen  , aveIden   = string.atoi(age[5]), string.atoi(age[6])
            nn   = re.search(r';?NR=([^;]+)' , col[7])
            n    = int(1000 * string.atof(nn.group(1)) + 0.5) / 10.0

            if len( col[idx].split(':')[fmat['AA']].split(',') ) < 4 :
                print >> sys.stderr, col[0], col[1], fmat['AA'], col[idx].split(':')[fmat['AA']], col[idx]
                continue

            alt = string.atoi( col[idx].split(':')[fmat['AA']].split(',')[1] ) # Alternate perfect
            bot = string.atoi( col[idx].split(':')[fmat['AA']].split(',')[3] ) # Both imperfect
            if leg       not in distance : distance[leg]       = [0,0] 
            if leftIden  not in leftIdn  : leftIdn[leftIden]   = [0,0] 
            if rightIden not in rightIdn : rightIdn[rightIden] = [0,0] 
            if aveIden   not in aveIdn   : aveIdn[aveIden]     = [0,0] 
            if n         not in nr       : nr[n]               = [0,0]
            if alt       not in aa       : aa[alt]             = [0,0]
            if bot       not in bb       : bb[bot]             = [0,0]
            distance[leg][0]       += 1
            leftIdn[leftIden][0]   += 1
            rightIdn[rightIden][0] += 1
            aveIdn[aveIden][0]     += 1
            nr[n][0]               += 1
            aa[alt][0]             += 1
            bb[bot][0]             += 1
 
            data.append([leg, alt, leftIden, rightIden, aveIden, n, bot])
            
            print col[0], '\t', col[1], '\t', leg, '\t', leftIden, '\t', rightIden, '\t', aveIden, '\t', n, '\t', alt, '\t', bot, '\t',svsize,qId, qFaLen[qId]
    I.close()

    data = np.array(data)
    print >> sys.stderr, '\nPosition\tALTernatePerfect\tLeftIdentity\tRightIdentity\tAveIden\tNRatio\tBothImperfect'
    print >> sys.stderr, 'Means: ', data.mean(axis=0), '\nstd  : ', data.std(axis=0), '\nMedian: ', np.median( data, axis=0 ) 
    print >> sys.stderr, '25 Percentile:', np.percentile(data, 25,axis=0), '\n50 Percentile:', np.percentile(data, 50,axis=0), '\n75 Percentile:', np.percentile(data, 75,axis=0)

    DrawFig( figPrefix, \
             np.array (Accum( distance )), \
             np.array (Accum( leftIdn  )), \
             np.array (Accum( rightIdn )), \
             np.array (Accum( aveIdn   )), \
             np.array (Accum( nr, True )), \
             np.array (Accum( aa       )), \
             np.array (Accum( bb       )), \
             np.array ( data ) )


if __name__ == '__main__' :

    main(sys.argv[1:])


