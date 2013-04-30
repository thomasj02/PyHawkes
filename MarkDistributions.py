__author__ = 'tjohnson'
import random

class ParetoMarkDistribution:
    def __init__(self,params):
        """
        Pareto-distributed impact function
        Must have rho>2
        Seems to be different than standard pareto distribution??
        From Liniger thesis p. 37
        """
        self.setParams(params)

    @staticmethod
    def getNumParameters():
        return 5

    @staticmethod
    def getParameterBounds():
        #This specifies that alpha>0, but strictly speaking only one of alpha, beta, gamma have to be greater than zero.
        #See Liniger thesis p. 34
        return [[0,None],[2,None],[0,None],[None,None],[None,None]]

    def setParams(self, params):
        """
        Set parameters with an iterable to support log-likelihood calculation
        This will set:
         mu=params[0]
         rho=params[1]
         alpha=params[2]
         beta=params[3]
         gamma=params[4]
        """

        self.mu = params[0]
        self.rho = params[1]
        self.alpha=params[2]
        self.beta=params[3]
        self.gamma=params[4]

    def getRandomValue(self):
        randomVal=random.random()
        return self.__inverseCumulativeDistribution(randomVal)

    def __inverseCumulativeDistribution(self,x):
        """From wolfram alpha"""
        firstMultiplier=(1.0-x)**(1.0/self.rho)-1.0
        secondMultiplier=(1.0-x)**(-1.0/self.rho)
        retval=-self.mu*firstMultiplier*secondMultiplier
        return retval

    def getDensityFunction(self,x):
        """Liniger f_u,p(x) (p. 21)"""
        numerator=self.rho*self.mu**self.rho
        denominator=(x+self.mu)**(self.rho+1)
        return numerator/denominator

    def getCumulativeDistributionFunction(self,x):
        """Liniger F_u,p(x)"""
        return 1.0-(self.mu/(x+self.mu))**self.rho

    def getImpactFunction(self,x):
        """Impact function for pareto distribution. Liniger g_k(x)"""
        term1Numerator=(self.rho-1.0)*(self.rho-2.0)
        term1Denominator=self.alpha*(self.rho-1.0)*(self.rho-2.0)+self.beta*self.mu*(self.rho-2.0)+2.0*self.gamma*self.mu*self.mu
        term2=self.alpha+self.beta*x+self.gamma*x*x

        return term1Numerator/term1Denominator*term2


class VoidMarkDistribution:
    def __init__(self,params):
        """
        Void impact function
        From Liniger thesis p. 21
        """
        self.setParams(params)

    @staticmethod
    def getNumParameters():
        return 0

    @staticmethod
    def getParameterBounds():
        return []

    def setParams(self, params):
        pass

    def getRandomValue(self):
        return 1.0

    def getDensityFunction(self,x):
        return 1.0

    def getCumulativeDistributionFunction(self,x):
        return 1.0

    def getImpactFunction(self,x):
        return 1.0
