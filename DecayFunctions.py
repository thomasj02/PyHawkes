__author__ = 'tjohnson'
import math

class ExponentialDecayFunction:
    """
    Exponential decay function from Liniger thesis p. 32
    """

    def __init__(self,params):
        self.setParams(params)
        self.epsilon=1e-5 #For Q function

    @staticmethod
    def getNumParameters():
        return 1

    @staticmethod
    def getParameterBounds():
        return [[0,None]]

    def setParams(self,params):
        """
        Set parameters with an iterable to support log-likelihood calculation
        This will set alpha=params[0]
        """
        self.alpha=params[0]

    def getW(self,t):
        """
        Get value of decay function

        :param componentIdx: Index of multivariate process component
        :param t: Time
        :return: Value of decay function w
        """
        return self.alpha*math.exp(-self.alpha*t)

    def getWBar(self,t):
        """Get value of cumulative decay function"""

        return 1.0-math.exp(-self.alpha*t)

    def getQ(self):
        """Get value of quantile function"""
        return -math.log(self.epsilon)/self.alpha
