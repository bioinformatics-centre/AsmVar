"""
===========================================
===========================================
Author : Shujia Huang
Date   : 2014-05-20 08:50:06
"""
import sys
import numpy as np
from sklearn import mixture
from sklearn.utils.extmath import logsumexp
# My own class
import VariantDatum as vd
import VariantRecalibratorArgumentCollection as VRAC

class VariantRecalibratorEngine :

    def __init__ (self, vrac=None ) :
        self.VRAC = VRAC.VariantRecalibratorArgumentCollection()
        if vrac : self.VRAC = vrac
        self.MIN_PROB_CONVERGENCE     = 2e-3
        self.MIN_ACCEPTABLE_LOD_SCORE = -20000.0
        
    def GenerateModel ( self, data, maxGaussians ) :
        if len(data)    == 0 : raise ValueError ( '[ERROR] No data found. The size is %d' %len(data) )
        if not isinstance(data[0], vd.VariantDatum) : 
            raise ValueError ('[ERROR] The data type should be "VariantDatum" in GenerateModel() of class VariantRecalibratorEngine(), but found %s'% str(type(data[0])) )
        if maxGaussians <= 0 : raise ValueError ( '[ERROR] maxGaussians must be a positive integer but found: %d' % maxGaussians )
        """
        gmm = mixture.GMM( n_components = maxGaussians, covariance_type='full', thresh = self.MIN_PROB_CONVERGENCE, 
                                 n_iter = self.VRAC.NITER , n_init = self.VRAC.NINIT, params='wmc', init_params='wmc' )
        gmm.fit( [ d.annotations for d in data ] )
        """
        gmms = [ mixture.GMM(n_components = n + 1, covariance_type='full', thresh = self.MIN_PROB_CONVERGENCE, 
                                   n_iter = self.VRAC.NITER , n_init = self.VRAC.NINIT, params='wmc', init_params='wmc') for n in range(maxGaussians) ]
        trainingData = np.array([d.annotations for d in data])
        minBIC, bics = np.inf, []
        for g in gmms : 
            print >> sys.stderr, '[INFO] Trying %d gaussian in GMM process ...' % g.n_components
            g.fit( trainingData ); bic = g.bic( trainingData )
            bics.append(bic)
            if bic != float('inf') or (bic < minBIC and g.converged_) : 
                bestgmm, minBIC = g, bic
        print >> sys.stderr, '[INFO] All the BIC:', bics
        print >> sys.stderr, '[INFO] Model Training Done. And take the model with %d gaussiones which with BIC %f.' % ( len(bestgmm.means_), minBIC )
        
        return bestgmm

    def EvaluateData ( self, data, gmm, evaluateContrastively = False ) :

        if not isinstance(data[0], vd.VariantDatum) : 
            raise ValueError ('[ERROR] The data type should be "VariantDatum" in EvaluateData() of class VariantRecalibratorEngine(), but found %s'% str(type(data[0])) )
        print >> sys.stderr, '[INFO] Evaluating full set of', len(data), 'variants...'
        for i,_ in enumerate( data ) : 

            thisLod = gmm.score( data[i].annotations[np.newaxis,:] ) / np.log(10)
            if np.math.isnan( thisLod ) : 
                gmm.converged_ = False
                return
            if evaluateContrastively :
                # contrastive evaluation: (prior + positive model - negative model)
                data[i].lod =  data[i].prior + data[i].lod - thisLod # Some problem to get understand of the data[i].prior
                if thisLod  == float( 'inf' ) : data[i].lod = self.MIN_ACCEPTABLE_LOD_SCORE * ( 1.0 + np.random.rand(1)[0] )
            else :
                data[i].lod = thisLod # positive model only so set the lod and return
        return self

    def CalculateWorstPerformingAnnotation ( self, data, goodModel, badModel ) : 

        for i, d in enumerate(data) :
            probDiff = [self.EvaluateDatumInOneDimension(goodModel, d, k) - self.EvaluateDatumInOneDimension(badModel, d, k) for k in range(len(d.annotations))]
            data[i].worstAnnotation = np.argsort( probDiff )[0] # Get the index of the worst annotations
        return self

    def EvaluateDatumInOneDimension ( self, gmm, datum, iii ) :

        pVarInGaussianLogE = [np.log(w) + NormalDistributionLoge(gmm.means_[k][iii], gmm.covars_[k][iii][iii], datum.annotations[iii]) for k,w in enumerate(gmm.weights_)]
        return logsumexp( np.array(pVarInGaussianLogE) ) / np.log( 10 ) # np.log10( Sum(pi_k * p(v|n,k)) )

def NormalDistributionLoge(mu, sigma, x) :

    if sigma <= 0 : raise ValueError ( '[ERROR] sd: Standard deviation of normal must be > 0 but found : %f' % sigma )
    if mu == float( 'inf' ) or mu == float( '-inf' ) or sigma == float( 'inf' ) or sigma == float( '-inf' ) or \
       x  == float( 'inf' ) or x  == float( '-inf' ) :
        raise ValueError ( '[ERROR] mean, sd, or, x : Normal parameters must be well formatted (non-INF, non-NAN)' )

    a = -1.0 * ( np.log(sigma) + 0.5 * np.log(2 * np.pi) )
    b = -0.5 * ( ( x - mu ) / sigma ) ** 2

    return a + b  # The Natural log

