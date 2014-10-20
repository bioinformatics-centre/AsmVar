"""
================================================
My own Gaussion Mixture Model for SV genotyping.
Learn form scikit-learn
================================================

Author : Shujia Huang
Date   : 2014-01-06 14:33:45

"""
import sys
import numpy as np
from scipy import linalg
from sklearn import cluster
from sklearn.base import BaseEstimator
from sklearn.utils.extmath import logsumexp

EPS = np.finfo(float).eps

class GMM ( BaseEstimator ) :

    """
    Copy from scikit-learn
    """

    def __init__(self, n_components=1, covariance_type='diag', random_state=None, thresh=1e-2, min_covar=1e-3,
                 n_iter=100, n_init=10, params='wmc', init_params='wmc'):
        self.n_components    = n_components
        self.covariance_type = covariance_type
        self.thresh          = thresh
        self.min_covar       = min_covar
        self.random_state    = random_state
        self.n_iter          = n_iter
        self.n_init          = n_init
        self.params          = params
        self.init_params     = init_params

        self.init_means      = []
        self.init_covars     = []
        self.category        = []   # For genotype

        if not covariance_type in ['spherical', 'tied', 'diag', 'full']:
            raise ValueError( 'Invalid value for covariance_type: %s' % covariance_type )

        if n_init < 1: raise ValueError('GMM estimation requires at least one run')

        self.weights_ = np.ones(self.n_components) / self.n_components

        # flag to indicate exit status of fit() method: converged (True) or
        # n_iter reached (False)
    def score_samples(self, X):
        """Return the per-sample likelihood of the data under the model.

        Compute the log probability of X under the model and
        return the posterior distribution (responsibilities) of each
        mixture component for each element of X.

        Parameters
        ----------
        X: array_like, shape (n_samples, n_features)
            List of n_features-dimensional data points. Each row
            corresponds to a single data point.

        Returns
        -------
        logprob : array_like, shape (n_samples,)
            Log probabilities of each data point in X.

        responsibilities : array_like, shape (n_samples, n_components)
            Posterior probabilities of each mixture component for each
            observation
        """
        X = np.asarray(X)
        if X.ndim == 1:
            X = X[:, np.newaxis]
        if X.size == 0:
            return np.array([]), np.empty((0, self.n_components))
        if X.shape[1] != self.means_.shape[1]:
            raise ValueError('The shape of X  is not compatible with self')

        lpr = (log_multivariate_normal_density(X, self.means_, self.covars_,self.covariance_type)
               + np.log(self.weights_))

        logprob = logsumexp(lpr, axis=1)
        responsibilities = np.exp(lpr - logprob[:, np.newaxis])
        return logprob, responsibilities

    def predict(self, X):
        """
        Predict label for data.

        Parameters
        ----------
        X : array-like, shape = [n_samples, n_features]

        Returns
        -------
        C : array, shape = (n_samples,)
        """
        logprob, responsibilities = self.score_samples(X)
        return responsibilities.argmax(axis=1)

    def predict_proba(self, X):
        """
        Predict posterior probability of data under each Gaussian
        in the model.

        Parameters
        ----------
        X : array-like, shape = [n_samples, n_features]

        Returns
        -------
        responsibilities : array-like, shape = (n_samples, n_components)
            Returns the probability of the sample for each Gaussian
            (state) in the model.
        """
        logprob, responsibilities = self.score_samples(X)
        return responsibilities

    def fit(self, X, sam2col, family):
        """
        Copy form scikit-learn: gmm.py
        Estimate model parameters with the expectation-maximization
        algorithm.

        A initialization step is performed before entering the em
        algorithm. If you want to avoid this step, set the keyword
        argument init_params to the empty string '' when creating the
        GMM object. Likewise, if you would like just to do an
        initialization, set n_iter=0.

        Parameters
        ----------
        X : array_like, shape (n, n_features)
            List of n_features-dimensional data points.  Each row
            corresponds to a single data point.
        """

        X = np.asarray(X, dtype=np.float)
        if X.ndim == 1:
            X = X[:, np.newaxis]
        if X.shape[0] < self.n_components:
            raise ValueError(
                'GMM estimation with %s components, but got only %s samples' %
                (self.n_components, X.shape[0]))

        lowest_bias = np.infty

        c1,c2,c3    = '1/1', '0/1', '0/0'
        m1,m2,m3    = 0.001 , 0.5 , 1.0
        v1,v2,v3    = 0.002, 0.002, 0.002
        category    = np.array([ [c1,c2,c3],
                                 [c1,c2], [c1,c3], [c2,c3] ,
                                 [c1]         , [c2]       , [c3] ])
        init_means  = np.array([ [[ m1],[ m2] , [ m3]], 
                                 [[ m1],[ m2]], [[m1],[m3]], [[m2],[m3]], 
                                 [[m1]]       , [[m2]]     , [[m3]] ])
        init_covars = np.array([ [[[ v1]],[[ v2]],[[ v3]]], 
                                 [[[ v1]],[[ v2]]], [[[ v1]],[[ v3]]], [[[ v2]],[[ v3]]], 
                                 [[[ v1]]]        , [[[ v2]]]        , [[[ v3]]] ])

        bestCovars, bestMeans, bestWeights, bestConverged, bestCategory = [], [], [], [], []
        for i, (m,v,c) in enumerate( zip(init_means, init_covars, category) ) : 

            if i == 0 and self.n_components != 3 : continue
            if i < 4  and self.n_components == 1 : continue
            self.init_means  = np.array(m)
            self.init_covars = np.array(v)
            self.category    = np.array(c)
            best_params,bias = self.training(X, sam2col, family)

            if lowest_bias > bias :
                lowest_bias         = bias
                bestCovars    = best_params['covars']
                bestMeans     = best_params['means']
                bestWeights   = best_params['weights']
                bestConverged = best_params['converged']
                bestCategory  = best_params['category']

            if self.n_components == 3 : break
            if self.n_components == 2 and i == 3 : break

        bestWeights = np.tile(1.0 / self.n_components, self.n_components)

        self.covars_     = bestCovars
        self.means_      = bestMeans
        self.weights_    = bestWeights
        self.converged_  = bestConverged
        self.category    = bestCategory

        return self

    ####
    def training(self, X, sam2col, family) :

        max_log_prob        = -np.infty
        lowest_bias         = np.infty

        wmin, wmax = 0.8, 1.2  # Factor intervel [wmin, wmax]
        for w in np.linspace(wmin, wmax, self.n_init) :
            if 'm' in self.init_params or not hasattr(self, 'means_'):
                #self.means_ = cluster.KMeans(n_clusters=self.n_components, random_state=self.random_state).fit(X).cluster_centers_
                self.means_  = w * self.init_means

            if 'w' in self.init_params or not hasattr(self, 'weights_'):
                self.weights_= np.tile(1.0 / self.n_components, self.n_components)

            if 'c' in self.init_params or not hasattr(self, 'covars_'):
                """
                cv = np.cov(X.T) + self.min_covar * np.eye(X.shape[1])
                if not cv.shape :
                    cv.shape = (1, 1)
                self.covars_ = distribute_covar_matrix_to_match_covariance_type(cv, self.covariance_type, self.n_components)
                """
                self.covars_ = self.init_covars

            # EM algorithms
            log_likelihood = []
            # reset self.converged_ to False
            self.converged_= False
            for i in range(self.n_iter):
                # Expectation step
                curr_log_likelihood, responsibilities = self.score_samples(X)
                log_likelihood.append(curr_log_likelihood.sum())

                # Check for convergence.
                if i > 0 and abs(log_likelihood[-1] - log_likelihood[-2]) < self.thresh :
                    self.converged_ = True
                    break
                #Maximization step
                self._do_mstep(X, responsibilities, self.params, self.min_covar)

            if   self.n_components == 3 :
                curr_bias =(self.means_[0][0]-self.init_means[0][0])+np.abs(self.means_[1][0]-self.init_means[1][0])+(self.init_means[2][0]-self.means_[2][0])
            elif self.n_components == 2 :
                curr_bias =np.abs(self.means_[0][0] - self.init_means[0][0]) + np.abs(self.init_means[1][0] - self.means_[1][0])
            elif self.n_components == 1 :
                curr_bias =np.abs (self.means_[0][0] - self.init_means[0][0])
            else :
                 print >> sys.stderr, '[ERROR] The companent could only between [1,3]. But yours is ', self.n_components
                 sys.exit(1)
            
            self.Label2Genotype()
            if w == wmin :
                max_log_prob = log_likelihood[-1]
                best_params  = { 'weights':self.weights_, 
                                 'means':self.means_, 
                                 'covars':self.covars_, 
                                 'converged':self.converged_, 
                                 'category':self.category }
                if self.converged_ :
                    lowest_bias = curr_bias

            if self.converged_ and lowest_bias > curr_bias :
                max_log_prob        = log_likelihood[-1]
                lowest_bias         = curr_bias
                best_params         = { 'weights': self.weights_, 
                                        'means': self.means_, 
                                        'covars': self.covars_,
                                        'converged': self.converged_, 
                                        'category':self.category }

        # check the existence of an init param that was not subject to
        # likelihood computation issue.
        if np.isneginf(max_log_prob) and self.n_iter:
            raise RuntimeError(
                "EM algorithm was never able to compute a valid likelihood " +
                "given initial parameters. Try different init parameters "   +
                "(or increasing n_init) or check for degenerate data." )

       # if neendshift :
       #     self.covars_    = tmp_params['covars']
       #     self.means_     = tmp_params['means']
       #     self.weights_   = tmp_params['weights']
       #     self.converged_ = tmp_params['converged']
       #     self.category   = tmp_params['category']

        return best_params, lowest_bias

    def _do_mstep(self, X, responsibilities, params, min_covar=0):
        """ 
        Perform the Mstep of the EM algorithm and return the class weihgts.
        """
        weights = responsibilities.sum(axis=0)
        weighted_X_sum  = np.dot(responsibilities.T, X)
        inverse_weights = 1.0 / (weights[:, np.newaxis] + 10 * EPS)

        if 'w' in params:
            self.weights_ = (weights / (weights.sum() + 10 * EPS) + EPS)
        if 'm' in params:
            self.means_ = weighted_X_sum * inverse_weights
        if 'c' in params:
            covar_mstep_func = _covar_mstep_funcs[self.covariance_type]
            self.covars_ = covar_mstep_func(self, X, responsibilities, weighted_X_sum, inverse_weights,min_covar)

        return weights

    """
    Here is just for genotyping process
    """
    # Decide the different guassion mu(mean) to seperate the genotype
    def Label2Genotype ( self ) :

        label2genotype = {}
        if self.converged_ :

            if len( self.means_ ) > 3 :
                print >> sys.stderr, 'Do not allow more than 3 components. But you set', len( self.means_ )
                sys.exit(1)

            for label,mu in enumerate( self.means_[:,0] ) :

                best_distance, bestIndx = np.infty, 0
                for i,m in enumerate(self.init_means[:,0]) :
                    distance = np.abs( mu - m )
                    if distance < best_distance :
                        bestIndx      = i
                        best_distance = distance

                label2genotype[label] = self.category[bestIndx]    

            # Put False if there are more than one 'label' points to the same 'genotype'
            g2c = { v:k for k,v in label2genotype.items() }
            if len(label2genotype) != len(g2c) : self.converged_ = False 
        else :
            label2genotype = { label: './.' for label in range( self.n_components ) }

        return label2genotype

    def Mendel(self, genotype, sample2col, family):

        ngIndx = []

        m,n,num = 0.0,0.0,0   # m is match; n is not match
        for k,v in family.items() :

            #if v[0] not in sample2col or v[1] not in sample2col : continue
            if k not in sample2col or v[0] not in sample2col or v[1] not in sample2col : continue
            if k not in sample2col :
                print >> sys.stderr, 'The sample name is not in vcf file! ', k
                sys.exit(1)

            # c1 is son; c2 and c3 are the parents
            c1,c2,c3 = genotype[ sample2col[k] ], genotype[ sample2col[v[0]] ], genotype[ sample2col[v[1]] ]

            if c1 == './.' or c2 == './.' or c3 == './.': continue
            num += 1;
            ng   = False
            if c2 == c3 :
                if c2 == '0/0' or c2 == '1/1' :
                    if c1 == c2 : m += 1
                    else : 
                        n += 1
                        ng = True
                else : # c2 == '0/1' and c3 == '0/1'
                    m += 1
            elif c2 == '0/1' and c3 == '1/1' :

                if c1 == '0/0' : 
                    n += 1
                    ng = True
                else : m += 1
            elif c2 == '0/1' and c3 == '0/0' :

                if   c1 == '1/1' : 
                    n += 1
                    ng = True
                else : m += 1
            elif c2 == '1/1' and c3 == '0/1' :

                if c1 == '0/0' : 
                    n += 1
                    ng = True
                else : m += 1
            elif c2 == '1/1' and c3 == '0/0' :

                if c1 == '1/1' or c1 == '0/0': 
                    n += 1
                    ng = True
                else : m += 1
            elif c2 == '0/0' and c3 == '0/1' :

                if c1 == '1/1' : 
                    n += 1
                    ng = True
                else : m += 1
            elif c2 == '0/0' and c3 == '1/1' :

                if c1 == '0/0' or c1 == '1/1' : 
                    n += 1
                    ng = True
                else : m += 1

            if ng : 
                ngIndx.append(sample2col[k])
                ngIndx.append(sample2col[v[0]])
                ngIndx.append(sample2col[v[1]])

        return m,n,num,set(ngIndx)
###
def log_multivariate_normal_density (X, means, covars, covariance_type='full'):
    """
    Log probability for full covariance matrices.
    """
    X = np.asarray(X)
    if X.ndim == 1:
        X = X[:, np.newaxis]
    if X.size == 0:
        return np.array([])
    if X.shape[1] != means.shape[1]:
        raise ValueError('The shape of X  is not compatible with self')

    log_multivariate_normal_density_dict = {
        'full' : _log_multivariate_normal_density_full
    }

    return log_multivariate_normal_density_dict[covariance_type]( X, means, covars )

def _log_multivariate_normal_density_full(X, means, covars, min_covar=1.e-7):
    """
    Log probability for full covariance matrices.
    """
    if hasattr(linalg, 'solve_triangular'):
        # only in scipy since 0.9
        solve_triangular = linalg.solve_triangular
    else:
        # slower, but works
        solve_triangular = linalg.solve
    n_samples, n_dim = X.shape
    nmix = len(means)
    log_prob = np.empty((n_samples, nmix))
    for c, (mu, cv) in enumerate(zip(means, covars)):
        try:
            cv_chol = linalg.cholesky(cv, lower=True)
        except linalg.LinAlgError:
            # The model is most probabily stuck in a component with too
            # few observations, we need to reinitialize this components
            cv_chol = linalg.cholesky(cv + min_covar * np.eye(n_dim),
                                      lower=True)
        
        cv_log_det = 2 * np.sum(np.log(np.diagonal(cv_chol)))
        cv_sol = solve_triangular(cv_chol, (X - mu).T, lower=True).T
        log_prob[:, c] = - .5 * (np.sum(cv_sol ** 2, axis=1) +
                                 n_dim * np.log(2 * np.pi) + cv_log_det)

    return log_prob


def distribute_covar_matrix_to_match_covariance_type( tied_cv, covariance_type, n_components) :
    """
    Create all the covariance matrices from a given template
    """
    if covariance_type == 'spherical':
        cv = np.tile(tied_cv.mean() * np.ones(tied_cv.shape[1]),
                     (n_components, 1))
    elif covariance_type == 'tied':
        cv = tied_cv
    elif covariance_type == 'diag':
        cv = np.tile(np.diag(tied_cv), (n_components, 1))
    elif covariance_type == 'full':
        cv = np.tile(tied_cv, (n_components, 1, 1))
    else:
        raise ValueError("covariance_type must be one of " +
                         "'spherical', 'tied', 'diag', 'full'")
    return cv

def _covar_mstep_full(gmm, X, responsibilities, weighted_X_sum, norm, min_covar ):
    """Performing the covariance M step for full cases"""
    # Eq. 12 from K. Murphy, "Fitting a Conditional Linear Gaussian
    # Distribution"
    n_features = X.shape[1]
    cv = np.empty((gmm.n_components, n_features, n_features))
    for c in range(gmm.n_components):
        post = responsibilities[:, c]
        # Underflow Errors in doing post * X.T are  not important
        np.seterr(under='ignore')
        avg_cv = np.dot(post * X.T, X) / (post.sum() + 10 * EPS)
        mu = gmm.means_[c][np.newaxis]
        cv[c] = (avg_cv - np.dot(mu.T, mu) + min_covar * np.eye(n_features))
    return cv

_covar_mstep_funcs = { 'full': _covar_mstep_full }


