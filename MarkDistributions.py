__author__ = 'tjohnson'
import random
import math

class ParetoMarkDistribution:
    def __init__(self,params):
        """
        Pareto-distributed impact function
        Must have rho>2
        From Liniger thesis p. 37
        Seems to be different than standard pareto distribution??
        """
        self.setParams(params)

    @staticmethod
    def getNumParameters():
        return 5

    @staticmethod
    def getParameterBounds():
        #This specifies that alpha, is positive, but strictly speaking only one of alpha, beta, gamma have to be greater than zero.
        #See Liniger thesis p. 34
        return [[1e-5,None],[2.0+1e-5,None],[1e-5,None],[0,None],[0,None]]

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


class _NonNormalizedExponentialImpactFunction3Params:
    def __init__(self,params):
        """
        Impact function for exponential distribution
        Liniger p.35
        """
        self.setParams(params)

    @staticmethod
    def getNumParameters():
        return 3

    @staticmethod
    def getParameterBounds():
        return [[1e-5,None]]*3

    def setParams(self,params):
        self.alpha=params[0]
        self.beta=params[1]
        self.gamma=params[2]

    def getImpactFunction(self,lambd,x):
        term1Numerator=lambd*lambd
        term1Denominator=self.alpha*lambd*lambd+self.beta*lambd+2.0*self.gamma
        term2=self.alpha+self.beta*x+self.gamma*x*x

        return (term1Numerator/term1Denominator)*term2


class _NonNormalizedExponentialImpactFunction1Params:
    def __init__(self, params):
        """
        Impact function for exponential distribution
        Liniger p.35
        """
        self.setParams(params)

    @staticmethod
    def getNumParameters():
        return 1

    @staticmethod
    def getParameterBounds():
        return [[1e-5, None]]

    def setParams(self, params):
        self.alpha = params[0]

    def getImpactFunction(self, lambd, x):
        term1Numerator=lambd**self.alpha
        term1Denominator=math.gamma(self.alpha+1.0)
        term2=x**self.alpha
        return (term1Numerator/term1Denominator)*term2



class _GenericExponentialMarkDistribution:
    def __init__(self,params,impactFunctionClass):
        """
        Exponential-distributed impact function
        Must have lambda>0
        Pass either _NonNormalizedExponentialImpactFunction3Params or _NonNormalizedExponentialImpactFunction1Param as impactFunctionClass
        See Liniger thesis p.35
        """
        self.impactFunctionClass=impactFunctionClass
        self.setParams(params)


    def getNumParameters(self):
        return 1+self.impactFunctionClass.getNumParameters()

    def getParameterBounds(self):
        parameterBounds=[[1e-5,None]]+self.impactFunctionClass.getParameterBounds()
        return parameterBounds

    def setParams(self,params):
        self.lambd=params[0]
        self.impactFunctionClass.setParams(params[1:])

    def getRandomValue(self):
        return random.expovariate(self.lambd)

    def getDensityFunction(self,x):
        return self.lambd*math.exp(-self.lambd*x)

    def getCumulativeDistributionFunction(self,x):
        return 1.0-math.exp(-self.lambd*x)

    def getImpactFunction(self,x):
        return self.impactFunctionClass.getImpactFunction(self.lambd,x)

class ExponentialMarkDistribution1Param:
    def __init__(self,params):
        self.impactFunction=_NonNormalizedExponentialImpactFunction1Params(params[1:])
        self.delegateMarkDistribution=_GenericExponentialMarkDistribution(params,self.impactFunction)

        self.setParams(params)

    @staticmethod
    def getNumParameters():
        return _NonNormalizedExponentialImpactFunction1Params.getNumParameters()+1

    @staticmethod
    def getParameterBounds():
        return [[1e-5,None]]+_NonNormalizedExponentialImpactFunction1Params.getParameterBounds()

    def setParams(self,params):
        self.delegateMarkDistribution.setParams(params)

    def getRandomValue(self):
        return self.delegateMarkDistribution.getRandomValue()

    def getDensityFunction(self,x):
        return self.delegateMarkDistribution.getDensityFunction(x)

    def getCumulativeDistributionFunction(self,x):
        return self.delegateMarkDistribution.getCumulativeDistributionFunction(x)

    def getImpactFunction(self,x):
        return self.delegateMarkDistribution.getImpactFunction(x)


