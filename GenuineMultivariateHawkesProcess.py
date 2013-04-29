import math

__author__ = 'tjohnson'
import numpy as np
import random

class GenuineMultivariateHawkesProcess:
    def __init__(self,immigrationIntensities,branchingMatrix,decayFunctions,markDistributions):
        self.numComponents=len(immigrationIntensities)
        self.immigrationIntensities=immigrationIntensities #nu_j in Liniger thesis
        self.branchingMatrix=branchingMatrix #script v_jk(x) in Liniger thesis
        self.decayFunctions=decayFunctions #w_j(t) in Liniger thesis
        self.markDistributions=markDistributions

    def getLambda(self,componentIdx,timeComponentMarkTriplets,t,includeT=False):
        """Returns lambda hat as defined in algorithm 1.28 at bottom of Liniger p. 41"""
        lambdaHat=self.immigrationIntensities[componentIdx]

        j=componentIdx
        decayFunction=self.decayFunctions[j]
        quantile=decayFunction.getQ()

        for s,k,x in timeComponentMarkTriplets:
            timeSinceEvent=t-s
            if(timeSinceEvent>quantile or timeSinceEvent<0):
                continue
            if not includeT and timeSinceEvent==0:
                continue

            branchingFactor=self.branchingMatrix[j,k]
            decayFunctionValue=decayFunction.getW(t-s)
            impactFunctionValue=self.markDistributions[k].getImpactFunction(x)

            lambdaHat+=branchingFactor*decayFunctionValue*impactFunctionValue

        return lambdaHat

    def __simulationInnerLoop(self,componentIdx,previousTime,timeComponentMarkTriplets):
        """Linniger thesis bottom p. 30"""

        tau=previousTime
        lambd=self.getLambda(componentIdx,timeComponentMarkTriplets,previousTime,includeT=True)
        while True:

            E=random.expovariate(1.0)
            tau=tau+E/lambd
            lambdaNew=self.getLambda(componentIdx,timeComponentMarkTriplets,tau)

            U=random.random()
            u=U*lambd
            if u<=lambdaNew:
                return tau

    def getLogLikelihood(self,timeComponentMarkTriplets):
        """
        Liniger thesis, Algorithm 1.27, p. 41
        """

        lambdaTermSum=0
        markDensityTermSum=0

        for t,d,x in timeComponentMarkTriplets:
            lambdaTermSum+=math.log(self.getLambda(d,timeComponentMarkTriplets,t,includeT=False))
            markDensityTermSum+=math.log(self.markDistributions[d].getDensityFunction(x))

        compensatorTermSum=0
        for j in range(0,self.numComponents):
            compensatorTermSum+=self.getCompensator(j,timeComponentMarkTriplets)

        print "lambdaTermSum: %s markDensityTermSum: %s compensatorTermSum: %s" % (lambdaTermSum,markDensityTermSum,compensatorTermSum)
        return lambdaTermSum+markDensityTermSum-compensatorTermSum

    def getCompensator(self,componentIdx,timeComponentMarkTriplets):
        """
        Big Lambda from Liniger Thesis
        Algorithm from Liniger thesis p. 44
        """
        #TODO: This can be made a little faster using the approximation at the bottom of Liniger p. 44.
        #TODO: Just don't calculate wBar if lastTime-s > q_j
        firstTime=timeComponentMarkTriplets[0][0]
        lastTime=timeComponentMarkTriplets[-1][0]
        firstTerm=self.immigrationIntensities[componentIdx]*(lastTime-firstTime)

        j=componentIdx
        secondTerm=0
        for s,k,x in timeComponentMarkTriplets:
            branchingFactor=self.branchingMatrix[j,k]
            wBarValue=self.decayFunctions[j].getWBar(lastTime-s)
            gValue=self.markDistributions[k].getImpactFunction(x)
            secondTerm+=branchingFactor*wBarValue*gValue

        retval=firstTerm+secondTerm
        return retval

    def simulate(self,numTimesteps):
        timeComponentMarkTriplets=[]

        currentTime=0
        for timestep in range(0,numTimesteps):
            newTime=float('inf')
            newComponent=0

            for j in range(0,self.numComponents):
                tau_n_j=self.__simulationInnerLoop(j,currentTime,timeComponentMarkTriplets)
                if tau_n_j<newTime:
                    newTime=tau_n_j
                    newComponent=j

            newMark=self.markDistributions[newComponent].getRandomValue()
            newTriple=(newTime,newComponent,newMark)
            timeComponentMarkTriplets.append(newTriple)
            currentTime=newTime

        return timeComponentMarkTriplets





