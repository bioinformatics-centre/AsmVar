"""
===================================================
Use Guassion Misture Model to estimate SV genotype 
from population data. VCF Format
===================================================

Author : Shujia Huang & Siyang Liu
Date   : 2013-12-14 10:48:02
Modify : 2013-12-16 13:58:57 Debuging

"""
import optparse
import os
import re
import string
import sys

import numpy as np
import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
from matplotlib.colors import LogNorm
import random as rd

import GMM2D as mixture
#from sklearn import mixture

def GetPPrate(fmat, formatInfo):

    trainIndx, ppr, pp = [], [], []
    theta = 1.0

    for i, t in enumerate(formatInfo):
        # 0/1:52,3:6,13,0,0,0,1,0:6,14:20-121512-121512:13:INS

        fi = t.split(':')
        if t != './.':
            rr = string.atof(fi[fmat['AA']].split(',')[0]) 
            aa = string.atof(fi[fmat['AA']].split(',')[1])
        else:
            rr, aa = 0, 0
        r  = theta
        if rr + aa >  0: 
            r = rr/(rr + aa)
            trainIndx.append(i)

        ppr.append([r])
        pp.append([rr, aa])

    #return range(len(ppr)), np.array(ppr), np.array(pp) # All for training!
    return trainIndx, np.array(ppr), np.array(pp)

###################

def UpdateInfoFromGMM(gmm, ppr, grey, red, green, blue, data, sam2col, family): 
# gmm : It's the GMM model
# data: It's the vcf line

    child = []
    for k,v in family.items():
        if k not in sam2col or v[0] not in sam2col or v[1] not in sam2col: continue
        child.append(sam2col[k])
    child = set(child)

    # determine the relationship of predictLabel and the SV genotype by the result of gmm.means_
    c2g = gmm.Label2Genotype()
    g2c = {v:k for k,v in c2g.items()}
    if not gmm.converged_: 
        if data[6] == '.' or data[6] == 'PASS':
            data[6] = 'FALSE_GENOTYPE'
        else:
            data[6] = 'FALSE_GENOTYPE;' + data[6]
    
    predict      = gmm.predict(ppr)
    predictProba = gmm.predict_proba(ppr)
    weights      = gmm.weights_

    #print >> sys.stderr, '### Predict:', predict, '### PPR:\n', ppr

    genotypeQuality = []
    for p,i in enumerate(predict):
        pd = predictProba[p][i]
        if pd > PRECISION: pd = PRECISION
        genotypeQuality.append(int(-10 * np.log10(1.0 - pd) + 0.5))

    logprob, posteriors = gmm.score_samples( ppr )
    for i in range(len(posteriors)): 
        posteriors[i][posteriors[i] == 0.0] = 1.0 - PRECISION

    loglhd = logprob[:, np.newaxis] + np.log( posteriors ) - np.log(weights) # The Log(e) genotype likelihoods 
    lhd    = np.exp(loglhd)
    for i in range(len(lhd)): lhd[i][lhd[i] == 0.0] = 1.0 - PRECISION
    loglhd = np.log10(lhd) # change the radices

    # Add info to the format fields for each individual
    fmat = {t:i for i,t in enumerate(data[8].split(':'))}
    if 'GQ' not in fmat : data[8] += ':GQ'
    if 'PL' not in fmat : data[8] += ':PL'

    first, gnt, ac, iac, N = True, [], 0, 0, 0
    refCount, hetCount, homCount = 1, 1, 1 # Assign 1 to prevent 0 in denominator 
    for i in range(len(data[9:])):

        if data[9+i] == './.':
            for j in range(len(fmat) - 1): data[9+i] += ':.'
            data[9+i] += ':0:65535,65535,65535'
            gnt.append('./.')
            continue
        
        fi = data[9+i].split(':')
        # Change raw genotype
        gt = 0
        if fi[0] != './.' and fi[0] != '0/0' :
            tt = [string.atoi(g) for g in fi[0].split('/')]
            if tt[0] > 1: 
                gt = tt[0]
            elif tt[1] > 1: 
                gt = tt[1]
        
        # Change the genotype
        if gt > 0:
            tt = [string.atoi(g) for g in c2g[predict[i]].split('/')]
            if tt[0] > 0: tt[0] = gt
            if tt[1] > 0: tt[1] = gt
            fi[0] = str(tt[0]) + '/' + str(tt[1])
        else: 
            fi[0] = c2g[predict[i]]
        pl = [int(p+0.5) for p in -10 * loglhd[i]]
        gnt.append(fi[0])

        if fi[0] == '1/1': 
            if i not in child: iac += 2
            ac       += 2
            homCount += 1 
        if fi[0] == '0/1': 
            if i not in child : iac += 1
            ac       += 1
            hetCount += 1
        if fi[0] == '0/0': 
            if i not in child : iac += 0
            ac       += 0
            refCount += 1
        if fi[0] != './.': N += 1

        if '0/0' not in g2c:
            phredScale        = '65535'
            if first : weight = '0'
        else:
            phredScale        = str(pl[g2c['0/0']])
            if first : weight = str(weights[g2c['0/0']])

        if '0/1' not in g2c:
            phredScale        += ',65535'
            if first : weight += ',0'
        else :
            phredScale        += ',' + str(pl[g2c['0/1']])
            if first : weight += ',' + str(weights[g2c['0/1']])

        if '1/1' not in g2c :
            phredScale        += ',65535'
            if first : weight += ',0'
        else :
            phredScale       += ',' + str(pl[g2c['1/1']])
            if first: weight += ',' + str(weights[g2c['1/1']])
        first = False

        if 'GQ' not in fmat: 
            fi.append(str(genotypeQuality[i]))
        else: 
            fi[fmat['GQ']] = str(genotypeQuality[i])

        if 'PL' not in fmat: 
            fi.append(phredScale)
        else: 
            fi[fmat['PL']] = phredScale

        data[9+i] = ':'.join(fi)

    if N == 0: N = 1
    p = float(2.0 * refCount + hetCount) / (2.0 * (refCount + hetCount + homCount)) # expected reference allele frequency
    q = 1.0 - p                              # expected alternative allele frequency
    f = 1.0 - (hetCount / (2.0 * p * q * N)) # The hetCount VS expected of hetCount

    """
    if re.search(r';K=([^;]+)', data[7]): 
        data[7] = re.sub(r';K=([^;]+)', ';K=' + str(len(g2c)), data[7])
    else: 
        data[7] += ';K=' + str(len(g2c))
   
    if re.search(r';W=([^;]+)', data[7]): 
        data[7] = re.sub(r';W=([^;]+)', ';W=' + weight, data[7])
    else: 
        data[7] += ';W=' + weight
    """

    if re.search(r';InbCoeff=([^;]+)', data[7]): 
        data[7] = re.sub(r';InbCoeff=([^;]+)', ';InbCoeff=%.2f' % f, data[7])
    else: 
        data[7] += ';InbCoeff=%.2f' % f

    if homCount == N + 1 or refCount == N + 1: # The all sample are totally 1/1 or 0/0! 
        if data[6] == '.' or data[6] == 'PASS':
            data[6] = 'FALSE_GENOTYPE'
        elif 'FALSE_GENOTYPE' not in data[6]:
            data[6] = 'FALSE_GENOTYPE;' + data[6]
        gmm.converged_ = False

    ng = set([])
    if gmm.n_components > 1 and gmm.converged_: 
        _, _, _, ng = gmm.Mendel(gnt, sam2col, family)

    if gmm.converged_:

        for i,g in enumerate([c2g[c] for c in predict]): # For figure

            grey.append([ppr[i]])
            if gnt[i] == './.': continue
            if g == '0/0':
                if i in ng:
                    green.append([ppr[i]])
                else: 
                    green.append([ppr[i]])
            if g == '0/1':
                if i in ng:
                    red.append([ppr[i]])
                else: 
                    red.append([ppr[i]])
            if g == '1/1':
                if i in ng:
                    blue.append([ppr[i]])
                else: 
                    blue.append([ppr[i]])

    return np.array(genotypeQuality), np.array(gnt), ac, iac, f 

def DrawModel2 ( gmm, ppr ) :

    #fig = plt.figure()

    minv = np.min([-1.5 * np.max(ppr), 1.5 * np.max(ppr)])
    maxv = np.max([-1.5 * np.max(ppr), 1.5 * np.max(ppr)])
    x = np.linspace( minv, maxv )
    y = np.linspace( minv, maxv )
    X, Y = np.meshgrid(x, y)
    XX = np.c_[X.ravel(), Y.ravel()]
    Z  = np.log( np.abs(gmm.score_samples(XX)[0]) ) #np.log(-gmm.score_samples(XX)[0])
    Z  = Z.reshape(X.shape)

    CS = plt.contour(X, Y, Z)
    CB = plt.colorbar(CS, shrink=0.8, extend='both')
    plt.scatter(ppr[:, 0], ppr[:, 1], .8)

    plt.axis('tight')
    # plt.show()

def DrawModel(figureFile, gmm, ppr, pp):

    predict = gmm.predict(ppr)
    mu      = gmm.means_
    sigma   = gmm.covars_
    fig = plt.figure()
    c2g = gmm.Label2Genotype()

    color1   = ['bo', 'ro', 'go']
    labels   = ['Hom','Het','Ref']
    for i in range( len(mu) ) :

        x = []
        for j in range(len(predict)) :
           if predict[j] == i : x.append(pp[j])
        x = np.array(x)
        plt.plot(x[:,0], x[:,1], color1[i], label = labels[i])

    plt.legend()
    plt.xlabel('Reference depth', fontsize=14)
    plt.ylabel('Alternate depth', fontsize=14)

    """
    color    = ['r','g','b']
    for i,x in enumerate( (red, green, blue) ) :

        n, bins, patches = plt.hist(x, num_bins, normed=1, facecolor=color[i], alpha=1.0 )
        x = mu[i][0] + sigma[i][0] * np.random.randn(300)
        n, bins, patches = plt.hist(x, num_bins, normed=1, facecolor=color[i], alpha=0.0 )
        y = mlab.normpdf(bins, mu[i][0], sigma[i][0])
        plt.plot( bins, y, color2[i], lw=2 )
    """

    fig.savefig(figureFile + '.png')
    fig.savefig(figureFile + '.pdf')


def DrawFig1 ( figureFile, title, red, green, blue, grey ) :

    numbins = 50
    fig     = plt.figure(1)
    plt.title(title, fontsize=12)
    
    if len( green ) > 0 : n, bins, patches = plt.hist(green[:,0], numbins, normed=1, facecolor='green' , alpha= .8 )
    if len( blue  ) > 0 : n, bins, patches = plt.hist(blue[:,0] , numbins, normed=1, facecolor='blue'  , alpha= .8 )
    if len( red   ) > 0 : n, bins, patches = plt.hist(red[:,0]  , numbins, normed=1, facecolor='red'   , alpha= .8 )
    plt.xlim( 0, 1.3)
    plt.ylabel('Normed Number', fontsize=16)
    plt.xlabel('PEP-Ratio VS Reference', fontsize=16)
    fig.savefig(figureFile + '.png')
    fig.savefig(figureFile + '.pdf')

def DrawAC ( figPrefix, title, green, blue, red, power, alleleCount, ialleleCount, inbCoff, svsize, numbins = 60 ) :

    highNum = 3
    widNum  = 2
    fig = plt.figure(num=None, figsize=(4 * widNum, 3 * highNum), facecolor='w', edgecolor='k')
    if len( alleleCount ) > 0 :
        plt.subplot(321)
        plt.title('AC Number', fontsize=12)
        plt.hist(alleleCount , numbins, histtype='bar', normed=1, facecolor = 'c', color=['c'] )
        plt.ylabel('#', fontsize=12)
    if len ( ialleleCount ) > 0 :
        plt.subplot(322)
        plt.title('IAC Number', fontsize=12)
        plt.hist(ialleleCount , numbins, histtype='bar', normed=1, facecolor = 'c', color=['c'] )
        plt.ylabel('#', fontsize=12)
    if len( inbCoff ) > 0 : 
        plt.subplot(323)
        plt.hist(inbCoff, 20 ,normed=1, facecolor = 'c', color=['c'] )
        plt.xlabel('1.0-hetCount/Expected_hetCount',fontsize=12)
        plt.ylabel('#', fontsize=14)

    plt.subplot(324)
    plt.title(title, fontsize=12)
    if len( green ) > 0 : plt.hist(green[:,0], numbins, normed=1, facecolor='green' , color=['g'], label=['Het'], alpha= .8 )
    if len( blue  ) > 0 : plt.hist(blue[:,0] , numbins, normed=1, facecolor='blue'  , color=['b'], label=['Hom'], alpha= .8 )
    if len( red   ) > 0 : plt.hist(red[:,0]  , 30     , normed=1, facecolor='red'   , color=['r'], label=['Ref'], alpha= .8 )
    plt.legend()
    plt.xlim( 0, 2.0)
    plt.ylabel('Normed Number', fontsize=12)
    plt.xlabel('PEPE Ratio'   , fontsize=12)

    if len ( power ) > 0 :
        plt.subplot(313)
        plt.bar( np.arange(len(power)), power[:,1], color = 'c' )
        plt.xticks( np.arange(len(power)) + 0.5, power[:,0], rotation = 90 )
        plt.ylim(0.0, 1.2)
        plt.xlabel( 'SV SIZE(bp)', fontsize=12 )
        plt.ylabel( 'Proper of Genotype', fontsize=12 )

    fig.savefig(figPrefix + '.png')
    fig.savefig(figPrefix + '.pdf')
    #plt.show()

def DrawFig ( figureFile, alleleCount, red, green, blue ) :

    numbins = len(set(alleleCount))
    fig = plt.figure( num=None, facecolor='w', edgecolor='k' )
    plt.subplot(211)
    if len( alleleCount ) > 0 :
        plt.subplot(321)
        plt.title('AC Number', fontsize=14)
        plt.hist(alleleCount , numbins, histtype='bar', normed=1, facecolor = 'c', color=['c'] )
        plt.ylabel('Normed Number', fontsize=14)

    plt.subplot(212)
    numbins = 60
    plt.title('Data Distribution', fontsize=14)
    if len( green ) > 0 : plt.hist(green[:,0], numbins, histtype='bar', normed=1, facecolor='green', color=['g'], label=['Ref'], alpha= .8 )
    if len( blue  ) > 0 : plt.hist(blue[:,0] , numbins, histtype='bar', normed=1, facecolor='blue' , color=['b'], label=['Hom'], alpha= .8 )
    if len( red   ) > 0 : plt.hist(red[:,0]  , numbins, histtype='bar', normed=1, facecolor='red'  , color=['r'], label=['Het'], alpha= .8 )
    plt.legend()
    plt.ylabel('Normed Number', fontsize=14)

    fig.savefig(figPrefix + '.png')
    fig.savefig(figPrefix + '.pdf')


##################
def LoadFamily(file):

    if len(file) == 0: return {}

    family = {}
    for line in open(file):
    #1006 1006-05 1006-01 1006-02 0 0
        line = line.strip('\n') # Cut the reture char at the end.
        col  = line.split()
        if col[1] in family.keys(): 
            print >> sys.stderr, 'The key is already in family! Your file may have the duplication sample name : ', col[1], '\n'
            sys.exit(1)

        family[col[1]] = [col[2], col[3]]

    return family


############################# Main Process #############################
def main(opt):

    gqSummary = {1: 0.0, 2: 0.0, 3: 0.0, 10: 0.0, 20: 0.0, 30: 0.0, 'sum': 0.0, 'Yes': 0.0, 'No': 0.0}
    family    = LoadFamily(opt.family)
    outPrefix = opt.outPrefix

    for f in infile: 
    
        outHandle = open(outPrefix + '.vcf', 'w')
        outFailGtyHandle = open(outPrefix + '.false_genotype.vcf', 'w')

        print >> sys.stderr, '# *** Reading File: ', f, ' ***\n'

        if f[-3:] == '.gz': 
            if len(opt.chroms) > 0:
                gzformat, I = True, os.popen("/home/siyang/Bin/software_pip/tabix-0.2.6/tabix -h %s %s " % (f, ' '.join(opt.chroms)))
            else:
                gzformat, I = True, os.popen("gzip -dc %s " % f)
        else:
            gzformat, I = False, open(f) 

        # 'grey', 'red', 'green' and 'blue' are used for draw figure
        grey, red, green, blue, inbCoff, svsize = [], [], [], [], [], []
        power, alleleCount, ialleleCount = {}, [], []
        sam2col = {}
        mdr, mde, mdr_t, mde_t, mde_n = 0.0, 0.0, 0.0, 0.0, 0
        while 1: 
            # Read 100000 lines at one time make effectly
            lines = I.readlines(100000)
            if not lines : break

            for line in lines: 

                line = line.strip('\n') # Cut the reture char at the end.
                col  = line.split()
                if re.search(r'^##FORMAT=<ID=GT', line):
                    outHandle.write(line + '\n')
                    outFailGtyHandle.write(line + '\n')
                    outHandle.write('##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality. -10*log10(1-p), p=w*N/(sigma(W*N)) N is gaussian density (p: The posterior probability of Genotype call is correct)">\n')
                    outFailGtyHandle.write('##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality. -10*log10(1-p), p=w*N/(sigma(W*N)) N is gaussian density (p: The posterior probability of Genotype call is correct)">\n')
                    outHandle.write('##FORMAT=<ID=PL,Number=1,Type=String,Description="Phred-scaled genotype likelihood. Rounded to the closest integer as defined in the VCF specification (The lower the better). The value calculate -10*log10(p), p is the predict posterior probability. And the order is : HOM_REF,HETE_VAR,HOM_VAR">\n')
                    outFailGtyHandle.write('##FORMAT=<ID=PL,Number=1,Type=String,Description="Phred-scaled genotype likelihood. Rounded to the closest integer as defined in the VCF specification (The lower the better). The value calculate -10*log10(p), p is the predict posterior probability. And the order is : HOM_REF,HETE_VAR,HOM_VAR">\n')
                    continue
                elif re.search(r'^##INFO=<ID=AC', line):
                    outHandle.write(line + '\n')
                    outFailGtyHandle.write(line + '\n')
                    outHandle.write('##FILTER=<ID=FALSE_GENOTYPE,Description="False in genotype process">\n')
                    outFailGtyHandle.write('##FILTER=<ID=FALSE_GENOTYPE,Description="False in genotype process">\n')
                    #outHandle.write('##INFO=<ID=K,Number=1,Type=String,Description="Number of genotype stats">\n')
                    #outFailGtyHandle.write('##INFO=<ID=K,Number=1,Type=String,Description="Number of genotype stats">\n')
                    #outHandle.write('##INFO=<ID=W,Number=1,Type=String,Description="The wieghts for each genotype stats. And the order is: HOM_REF,HETE_VAR,HOM_VAR">\n')
                    #outFailGtyHandle.write('##INFO=<ID=W,Number=1,Type=String,Description="The wieghts for each genotype stats. And the order is: HOM_REF,HETE_VAR,HOM_VAR">\n')
                    outHandle.write('##INFO=<ID=InbCoeff,Number=1,Type=Float,Description="Inbreeding coefficient: 1.0 - hetCount/Expected_hetCount">\n')
                    outFailGtyHandle.write('##INFO=<ID=InbCoeff,Number=1,Type=Float,Description="Inbreeding coefficient: 1.0 - hetCount/Expected_hetCount">\n')
                elif re.search(r'^##', line):
                    outHandle.write(line + '\n')
                    outFailGtyHandle.write(line + '\n')
                    continue
                elif re.search(r'^#', line):
                    outHandle.write(line + '\n')
                    outFailGtyHandle.write(line + '\n')
                    sam2col = {sam:i for i,sam in enumerate(col[9:])}
                    continue

                if (len(opt.chroms)>0) and (col[0] not in opt.chroms): continue

                fmat = {t:i for i,t in enumerate(col[8].split(':'))}
                if 'AA' not in fmat: continue

                if col[6] == 'PASS': col[6] == '.'

                trainIndx, ppr, pp = GetPPrate(fmat, col[9:])
                if len(trainIndx) < 10: 

                    gqSummary['No'] += 1.0
                    if col[6] == '.' or col[6] == 'PASS':
                        col[6] = 'FALSE_GENOTYPE'
                    elif 'FALSE_GENOTYPE' not in col[6]:
                        col[6] = 'FALSE_GENOTYPE;' + col[6]
                    
                    if 'GQ' not in fmat: col[8] += ':GQ'
                    if 'PL' not in fmat: col[8] += ':PL'

                    for i in range(len(col[9:])): 

                        if col[9+i] == './.':
                        
                            for j in range(len(fmat) - 1): col[9+i] += ':.'
                            col[9+i] += ':0:65535,65535,65535'
                        else:

                            fi    = col[9+i].split(':')
                            fi[0] = './.'

                            if 'GQ' not in fmat: 
                                fi.append('0')
                            else: 
                                fi[fmat['GQ']] = '0'

                            if 'PL' not in fmat: 
                                fi.append('65535,65535,65535')
                            else: 
                                fi[fmat['PL']] = '65535,65535,65535'
                            col[9+i] = ':'.join(fi)

                    outFailGtyHandle.write('\t'.join(col) + '\n')
                    continue

                nc  = 3
                clf = mixture.GMM(n_components=nc, n_iter=50, n_init=8, covariance_type='full', thresh=0.001, params='wmc')
                clf.fit(ppr[trainIndx])
                if not clf.converged_:
                    print >> sys.stderr, '#+++ Position:', col[0], col[1], "couldn't converge with 3 components in GMM. Now trying 2 components ... "
                    nc  = 2
                    clf = mixture.GMM(n_components=nc, n_iter=50, n_init=8, covariance_type='full', thresh=0.001, params='wmc')
                    clf.fit(ppr[trainIndx])
                if not clf.converged_:
                    print >> sys.stderr, '#+++ Position:', col[0], col[1], "couldn't converge with 2 components in GMM. Now trying 1 components ... "
                    nc  = 1
                    clf = mixture.GMM(n_components=nc, n_iter=50, n_init=8, covariance_type='full', thresh=0.001, params='wmc')
                    clf.fit(ppr[trainIndx])
                print >> sys.stderr, '#--> Position:', col[0], col[1], 'with', clf.n_components,'components. Converge information :', clf.converged_
                print >> sys.stderr, '# Means: \n', clf.means_, '\nCovars: \n', clf.covars_,'\nWeight', clf.weights_, '\n*************'

                genotypeQuality,gnt,ac,iac,ef = UpdateInfoFromGMM(clf, ppr, grey, red, green, blue, col, sam2col, family)
                inbCoff.append(ef)
                if ef > -0.7 and clf.converged_: 
                    alleleCount.append(ac)
                    ialleleCount.append(iac)
                else: 
                    if col[6] == '.' or col[6] == 'PASS':
                        col[6] = 'FALSE_GENOTYPE'
                    elif 'FALSE_GENOTYPE' not in col[6]:
                        col[6] = 'FALSE_GENOTYPE;' + col[6]
                    clf.converged_ = False

                fmat = {t:i for i,t in enumerate(col[8].split(':'))}
                for i in range(len(col[9:])):
                    fi = col[9+i].split(':')
                    if fi[0] != './.': 
                        rr = string.atof(fi[fmat['AA']].split(',')[0])
                        aa = string.atof(fi[fmat['AA']].split(',')[1])
                    else:
                        rr, aa = 0, 0

                    # still keep the genotype info in 'FALSE_GENOTYPE'
                    # if rr + aa == 0 or 'FALSE_GENOTYPE' in col[6]:
                    if rr + aa == 0:
                        fi[fmat['GQ']] = '0'
                        fi[fmat['GT']] = './.'
                        fi[fmat['PL']] = '.'
                        col[9+i] = ':'.join(fi)
                        gnt[i]   = fi[fmat['GT']]
                        genotypeQuality[i] = 0
                    if 'FALSE_GENOTYPE' in col[6]: gnt[i] = './.'

                """
                size = re.search ( r';SVSIZE=([^;]+)', col[7] )
                size = string.atoi(size.group(1))
                if size < 2000 : 
                    if clf.converged_ : svsize.append( size )
                    # Calculate the bins length
                    base = 10 ** int( np.log10(size) )
                    bins = size / base
                    if size % base != 0 : bins += 1 # Can increat 1!
                    length = bins * base
                    if length not in power : power[length] = [0,0]
                    if clf.converged_ : power[length][0] += 1
                    else              : power[length][1] += 1
                    ###
                """

                if clf.converged_ and len(family) > 0: 
                    sm, sn, snum, _ = clf.Mendel(gnt, sam2col, family)
                    mdr_t += sm 
                    mde_t += sn
                    mde_n += snum

                if 'FALSE_GENOTYPE' in col[6]:
                    outFailGtyHandle.write('\t'.join(col) + '\n')
                else:
                    outHandle.write('\t'.join(col) + '\n')
                #DrawModel(figPrefix, clf, ppr, pp)

                if clf.converged_: 
                    gqSummary['Yes'] += 1.0
                    gqSummary['sum'] += len(ppr)
                    gqSummary[10]    += len(genotypeQuality[genotypeQuality>=10])
                    gqSummary[20]    += len(genotypeQuality[genotypeQuality>=20])
                    gqSummary[30]    += len(genotypeQuality[genotypeQuality>=30])
                    gqSummary[clf.n_components] += 1.0
                else: 
                    gqSummary['No']  += 1.0
        I.close()
        outHandle.close()
        outFailGtyHandle.close()

        if gqSummary['sum'] == 0: gqSummary['sum'] = 1.0
        if gqSummary['Yes'] + gqSummary['No'] == 0: 
            gqSummary['Yes'] = 1.0
            gqSummary['No']  = 1e9
        if gqSummary['Yes'] == 0: gqSummary['Yes'] = 1.0
        mderr_t = '-'
        if mdr_t + mde_t > 0.0: mderr_t = mde_t/(mdr_t+mde_t)

        print >> sys.stderr, '\n******** Output Summary information ************************************\n'
        print >> sys.stderr, '# ** The count of positions which can be genotype:',int(gqSummary['Yes']),',',gqSummary['Yes']/(gqSummary['Yes']+gqSummary['No'])
        print >> sys.stderr, '# ** (Just for the genotype positions)The mendelian violation of', f, 'is : ',mderr_t, '\t', mde_n 
        print >> sys.stderr, '# ** (Just for the genotype positions)Proportion of 1 component : ', gqSummary[1] / gqSummary['Yes']
        print >> sys.stderr, '# ** (Just for the genotype positions)Proportion of 2 component : ', gqSummary[2] / gqSummary['Yes']
        print >> sys.stderr, '# ** (Just for the genotype positions)Proportion of 3 component : ', gqSummary[3] / gqSummary['Yes']
        print >> sys.stderr, '# ** (Just for the genotype positions)The ratio of Homo-Ref(0/0): ', len( green ) / gqSummary['sum']
        print >> sys.stderr, '# ** (Just for the genotype positions)The ratio of Hete-Var(0/1): ', len( red   ) / gqSummary['sum']
        print >> sys.stderr, '# ** (Just for the genotype positions)The ratio of Homo-Var(1/1): ', len( blue  ) / gqSummary['sum']
        print >> sys.stderr, '# ** (Just for the genotype positions)Genotype Quality >= 10 : ', gqSummary[10] / gqSummary['sum']
        print >> sys.stderr, '# ** (Just for the genotype positions)Genotype Quality >= 20 : ', gqSummary[20] / gqSummary['sum']
        print >> sys.stderr, '# ** (Just for the genotype positions)Genotype Quality >= 30 : ', gqSummary[30] / gqSummary['sum']

        DrawFig (figPrefix, np.array(alleleCount), np.array(red), np.array(green), np.array(blue))
        """
        DrawAC( figPrefix, 'Mendelian violation = ' + str( mderr ), np.array(red), np.array(green), np.array(blue), 
                np.array([ [k, float(v[0])/(v[0]+v[1])] for k,v in sorted(power.items(), key = lambda d:d[0]) ]  ), 
                np.array(alleleCount), 
                np.array(ialleleCount), 
                np.array( inbCoff ), np.array(svsize), 60 )
        """

################################
################################

if __name__ == '__main__':

    usage = "Usage : %prog [option] [vcfInfile] > Output"
    optp  = optparse.OptionParser(usage=usage)
    optp.add_option("-c", "--chr", dest="chroms", metavar="CHR", help="process only specified chromosomes, separated by ','. [default: all]\nexample: --chroms=chr1,chr2", default=[])
    optp.add_option("-p", "--ped", dest="family", metavar="PED", help="Family information. ", default=[])
    optp.add_option("-f", "--fig", dest="figure", metavar="FIG", help="The prefix of figure about the GMM.", default=[])
    optp.add_option("-o", "--out", dest="outPrefix", metavar="OUT", help="The prefix of output. [out]", default = 'out')

    opt, infile = optp.parse_args()

    figPrefix = 'test'
    if len(infile)     == 0: optp.error("Required at least one [vcfInfile]\n")
    if len(opt.figure) > 0: figPrefix = opt.figure
    if any(opt.chroms): opt.chroms = opt.chroms.split(',')

    COMPONENT_NUM = 3        # The tyoe number of genotype
    PRECISION     = 0.9999999999

    main(opt)
    print >> sys.stderr, '\n# [INFO] Closing the two Ouput files:\n  -- (1/4) %s' % (opt.outPrefix + '.vcf')
    print >> sys.stderr, '  -- (2/4) %s' % (opt.outPrefix + '.false_genotype.vcf')
    print >> sys.stderr, '  -- (3/4) %s' % (opt.figure + '.pdf')
    print >> sys.stderr, '  -- (4/4) %s' % (opt.figure + '.png')
    print >> sys.stderr, '******************************* ALL DONE *******************************'

