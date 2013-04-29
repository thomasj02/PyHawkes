__author__ = 'tjohnson'
import numpy as np
class ImmigrationDescendantParameters:
    def __init__(self,numComponents,params):
        self.numComponents=numComponents
        self.nu=np.zeros(numComponents) #Immigration intensities
        self.q=np.zeros([numComponents,numComponents]) #Branching matrix

        self.setParameters(params)

    @staticmethod
    def getNumParameters(numComponents):
        #Square matrix plus vector
        return numComponents*numComponents+numComponents

    def setParameters(self,params):
        for nuIdx,paramVal in zip(range(0,self.numComponents),params[0:self.numComponents]):
            self.nu[nuIdx]=paramVal

        for qIdx,paramVal in zip(range(0,self.numComponents*self.numComponents),params[self.numComponents:]):
            self.q.put(qIdx,paramVal)
