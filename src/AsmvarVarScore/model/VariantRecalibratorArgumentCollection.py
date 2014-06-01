"""
====================================
====================================
Author : Shujia Huang
Date   : 2014-05-21 18:03:28

"""
class VariantRecalibratorArgumentCollection :

    def __init__ (self) :
        self.NITER = 150
        self.NINIT = 100
        self.STD_THRESHOLD         = 10.0
        self.MIN_NUM_BAD_VARIANTS  = 1000
        self.MAX_NUM_TRAINING_DATA = 20000
        self.BAD_LOD_CUTOFF        = -5.0
        self.MAX_GAUSSIANS         = 8
        self.MAX_GAUSSIANS_FOR_NEGATIVE_MODEL = 4

