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

def DrawFig ( r, g, b, figureFile ) :

    fig = plt.figure()

    numbins = 80
    plt.subplot(211)
    plt.title('RR/(AA + RR) After alternate align', fontsize=14)
    if len( b ) > 0 : plt.hist(b[:,0], numbins, normed=1, facecolor='blue'  , alpha= .8)
    if len( g ) > 0 : plt.hist(g[:,0], numbins, normed=1, facecolor='green' , alpha= .8)
    if len( r ) > 0 : plt.hist(r[:,0], numbins, normed=1, facecolor='red'   , alpha= .8)
    plt.ylim(0, 10)
    plt.ylabel('Normed Number', fontsize=14)
    """
    if len( g ) > 0 : plt.plot(g[:,0], g[:,1], 'go')
    if len( r ) > 0 : plt.plot(r[:,0], r[:,1], 'ro')
    if len( b ) > 0 : plt.hist(b[:,0], b[:,1], 'bo')
    plt.xlabel('Ref_perfect', fontsize=14)
    plt.ylabel('Alt_perfect', fontsize=14)
    """

    plt.subplot(212)
    numbins = 160 
    plt.title('PEP Ratio', fontsize=14)
    if len( b ) > 0 : plt.hist(b[:,1], numbins, normed=1, facecolor='blue'  , alpha= .8)
    if len( g ) > 0 : plt.hist(g[:,1], numbins, normed=1, facecolor='green' , alpha= .8)
    if len( r ) > 0 : plt.hist(r[:,1], numbins, normed=1, facecolor='red'   , alpha= .8)
    plt.ylabel('Normed Number', fontsize=14)
    """
    if len( g ) > 0 : plt.plot(g[:,2], g[:,3], 'go')
    if len( r ) > 0 : plt.plot(r[:,2], r[:,3], 'ro')
    if len( b ) > 0 : plt.hist(b[:,2], b[:,3], 'bo')
    plt.xlabel('Proper', fontsize=14)
    plt.ylabel('Expect Proper', fontsize=14)
    """
    plt.xlim(0, 1.5)

    fig.savefig(figureFile + '.png')
    fig.savefig(figureFile + '.pdf')

def main ( argv ) :

    vcfInfile = argv[0]
    if vcfInfile[-3:] == '.gz' :
        I = os.popen( 'gzip -dc %s' % vcfInfile )
    else :
        I = open ( vcfInfile )

    r, g, b = [], [], []
    tot, err1, err2 = 0,0,0
    while 1 :
        lines = I.readlines(100000)
        if not lines : break;

        for line in lines :

            line = line.strip( '\n' )
            col  = line.split()
            if re.search(r'^#', line) : continue
#            if not re.search( r';Novel=\(Same', col[7] ) and not re.search( r';Novel=\(ESame', col[7] ) : continue
            
            # GT:EQ:ER:FA:FN:FV:QC:QP:QR:RC:RP:RR:VS:VT:AA
            fmat = {}
            for i,t in enumerate( col[8].split(':') ) : fmat[t] = i # Get Format

            for type in ['AA', 'ER', 'RP', 'VT', 'VS'] :
                if type not in fmat.keys() :
                    print >> sys.stderr, '[ERROR] The format of VCF file is not right which you input, it did not contian %s field' % type
                    sys.exit(1)

            fi    = col[9].split(':')
            rr,aa = string.atof( fi[fmat['AA']].split(',')[0] ), string.atof( fi[fmat['AA']].split(',')[1] )
            et,pt = fi[fmat['ER']], fi[fmat['RP']]
            ep,pp = string.atof(et.split(',')[0]), string.atof(pt.split(',')[0])
            svsize= string.atoi( fi[fmat['VS']] )
            if svsize < 10 : continue

            r1, r2 = 1.2, 1.2
            if rr + aa > 0.0 : r1 = rr / (rr + aa)
            if ep      > 0.0 : r2 = pp / ep       

            tot  += 1
            if re.search(r'homozygous_', col[7]) :
                r.append([r1, r2])
                #r.append( [rr/ep, aa/ep, pp, ep] )
                if r1 > 0.4: 
                    print '\t'.join( col )
                    err1 += 1
            elif re.search(r'heterozygous_', col[7]) :
                g.append([r1, r2])
                #g.append([rr/ep, aa/ep, pp, ep])
                if r1 < 0.2 : 
                    print '\t'.join( col )
                    err2 += 1
            else :
                b.append([r1, r2])
                #b.append([rr/ep, aa/ep, pp, ep])
    I.close()
    figureFile = 'AA'
    if len(argv) > 1 : figureFile = argv[1]

    DrawFig ( np.array(r), np.array(g), np.array(b), figureFile )
    print >> sys.stderr, '# Total :', tot, '; Homo Err :', err1, '; Heter Err :', err2, '; Ratio :', float(err1+err2)/tot


if __name__ == '__main__' :

    main(sys.argv[1:])

