__author__ = 'tjohnson'
import random

class ParetoMarkDistribution:
    def __init__(self,mu,rho,alpha,beta,gamma):
        """
        Pareto-distributed impact function
        Must have rho>2
        Seems to be different than standard pareto distribution??
        From Liniger thesis p. 37
        """
        self.mu=mu
        self.rho=rho
        self.alpha=alpha
        self.beta=beta
        self.gamma=gamma

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
        """Liniger f_u,p(x)"""
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
        term2=self.alpha*self.beta*x+self.gamma*x*x

        return term1Numerator/term1Denominator*term2
