"""
===============================================
===============================================
Author : Shujia Huang
Date   : 2014-05-23 11:21:53
"""

import sys
import VariantDataManager as vdm
import VariantRecalibratorEngine as vre
import VariantRecalibratorArgumentCollection as VRAC

class VariantRecalibrator : 

    def __init__ (self) :

        self.VRAC        = VRAC.VariantRecalibratorArgumentCollection()
        self.dataManager = vdm.VariantDataManager()
        self.engine      = vre.VariantRecalibratorEngine( self.VRAC )
    def OnTraversalDone (self, data ) :

        self.dataManager.SetData( data )
        self.dataManager.NormalizeData()

        # Generate the positive model using the training data and evaluate each variant
        positiveTrainingData = self.dataManager.GetTrainingData()
        print >> sys.stderr, '[INFO] Training the goodModel ...'
        goodModel            = self.engine.GenerateModel( positiveTrainingData, self.VRAC.MAX_GAUSSIANS )
        print >> sys.stderr, '[INFO] The converged information of goodModel is :', goodModel.converged_
        self.engine.EvaluateData ( self.dataManager.data, goodModel, False)

        # Generate the negative model using the worst performing data and evaluate each variant contrastively
        print >> sys.stderr, '[INFO] Training the badModel ...'
        negativeTrainingData = self.dataManager.SelectWorstVariants()
        badModel             = self.engine.GenerateModel( negativeTrainingData, min(self.VRAC.MAX_GAUSSIANS_FOR_NEGATIVE_MODEL, self.VRAC.MAX_GAUSSIANS))
        print >> sys.stderr, '[INFO] The converged information of badModel is :', badModel.converged_
        self.engine.EvaluateData( self.dataManager.data, badModel, True )

        if (not goodModel.converged_) or (not badModel.converged_) : raise ValueError ( '[ERROR] NaN LOD value assigned. Clustering with this few variants and these annotations is unsafe. Please consider raising the number of variants used to train the negative model (via --minNumBadVariants 5000, for example) or lowering the maximum number of Gaussians allowed for use in the model (via --maxGaussians 4, for example).' )

        # Find the VQSLOD cutoff values which correspond to the various tranches of calls requested by the user
        self.engine.CalculateWorstPerformingAnnotation( self.dataManager.data, goodModel, badModel )



