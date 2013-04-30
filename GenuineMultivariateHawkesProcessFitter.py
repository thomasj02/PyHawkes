__author__ = 'tjohnson'
import ImmigrationDescendantParameters
import random

class GenuineMultivariateHawkesProcessFitter:
    def __init__(self,hawkesProcess):
        self.hawkesProcess=hawkesProcess

        totalNumParams=0
        parameterBounds=[]

        numComponents=hawkesProcess.numComponents
        self.immigrationDescendantNumParams=ImmigrationDescendantParameters.ImmigrationDescendantParameters.getNumParameters(numComponents)
        totalNumParams+=self.immigrationDescendantNumParams
        parameterBounds.extend(ImmigrationDescendantParameters.ImmigrationDescendantParameters.getParameterBounds(numComponents))

        self.decayFunctionsNumParams=[]
        for decayFunction in self.__getDecayFunctions():
            thisDecayFunctionNumParams=decayFunction.getNumParameters()
            self.decayFunctionsNumParams.append(thisDecayFunctionNumParams)
            totalNumParams+=thisDecayFunctionNumParams
            parameterBounds.extend(decayFunction.getParameterBounds())

        self.markDistributionsNumParams=[]
        for markDistribution in self.__getMarkDistributions():
            thisMarkDistributionNumParams=markDistribution.getNumParameters()
            self.markDistributionsNumParams.append(thisMarkDistributionNumParams)
            totalNumParams+=thisMarkDistributionNumParams
            parameterBounds.extend(markDistribution.getParameterBounds())

        self.numParams=totalNumParams
        self.parameterBounds=parameterBounds

    def __getMarkDistributions(self):
        retval=[]

        lastMarkDistribution=None
        for markDistribution in self.hawkesProcess.markDistributions:
            if markDistribution is lastMarkDistribution:
                continue
            retval.append(markDistribution)
            lastMarkDistribution=markDistribution

        return retval

    def __getDecayFunctions(self):
        retval=[]

        lastDecayFunction=None
        for decayFunction in self.hawkesProcess.decayFunctions:
            if decayFunction is lastDecayFunction:
                continue
            retval.append(decayFunction)
            lastDecayFunction=decayFunction

        return retval

    def getInitialRandomVector(self):
        vector=[]
        for lowerBound,upperBound in self.parameterBounds:
            if lowerBound is None:
                lowerBound=0.0
            randomValue=random.uniform(lowerBound,lowerBound+1.0)
            vector.append(randomValue)

        return vector

    def setParameterValues(self,params):
        parameterStartIndex=0
        self.hawkesProcess.immigrationDescendantParameters.setParameters(params[parameterStartIndex:self.immigrationDescendantNumParams])
        parameterStartIndex+=self.immigrationDescendantNumParams

        for decayFunction,decayFunctionNumParams in zip(self.__getDecayFunctions(),self.decayFunctionsNumParams):
            decayFunctionParams=params[parameterStartIndex:parameterStartIndex+decayFunctionNumParams]
            parameterStartIndex+=decayFunctionNumParams

            decayFunction.setParams(decayFunctionParams)


        for markDistribution,markDistributionNumParams in zip(self.__getMarkDistributions(),self.markDistributionsNumParams):
            markDistributionParams = params[parameterStartIndex:parameterStartIndex + markDistributionNumParams]
            parameterStartIndex+=markDistributionNumParams

            markDistribution.setParams(markDistributionParams)


    def maximizeLikelihood(self,timeComponentMarkTriples,initialGuess):
        #likelihoodFunction=lambda x: self.hawkesProcess.
        pass


if __name__=="__main__":
    import DecayFunctions
    import MarkDistributions
    import GenuineMultivariateHawkesProcess
    import ImmigrationDescendantParameters
    import numpy as np


    q=np.matrix("[0.61 0.16;0.60 0.06]")
    nu=np.array([0.021,0.029])
    alpha=0.015
    rho1=5.6
    rho2=7.2
    mu1=3.6
    mu2=4.2
    phi1=0.47
    phi2=1.1
    psi1=0.22
    psi2=0.0

    immigrationDescendantParameters=\
        ImmigrationDescendantParameters.ImmigrationDescendantParameters(2,[0.021,0.029,0.61,0.16,0.6,0.06])
    print immigrationDescendantParameters.nu
    print immigrationDescendantParameters.q
    decayFunction=DecayFunctions.ExponentialDecayFunction([alpha])
    markDistribution1=MarkDistributions.ParetoMarkDistribution([mu1,rho1,phi1,psi1,0.0])
    markDistribution2=MarkDistributions.ParetoMarkDistribution([mu2,rho2,phi2,psi2,0.0])

    hawkesProcess=GenuineMultivariateHawkesProcess.GenuineMultivariateHawkesProcess(
        immigrationDescendantParameters,
        [decayFunction,decayFunction],
        [markDistribution1,markDistribution2])

    fitter=GenuineMultivariateHawkesProcessFitter(hawkesProcess)
    print fitter.numParams
    print fitter.parameterBounds
    print fitter.getInitialRandomVector()
