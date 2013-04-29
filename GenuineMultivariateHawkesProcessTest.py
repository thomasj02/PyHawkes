__author__ = 'tjohnson'

import DecayFunctions
import MarkDistributions
import GenuineMultivariateHawkesProcess
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
decayFunction=DecayFunctions.ExponentialDecayFunction(alpha)
markDistribution1=MarkDistributions.ParetoMarkDistribution(mu1,rho1,phi1,psi1,0.0)
markDistribution2=MarkDistributions.ParetoMarkDistribution(mu2,rho2,phi2,psi2,0.0)

hawkesProcess=GenuineMultivariateHawkesProcess.GenuineMultivariateHawkesProcess(
    nu,
    q,
    [decayFunction,decayFunction],
    [markDistribution1,markDistribution2])

timeComponentMarkTriples1000=hawkesProcess.simulate(1000)
secondHalfTimeComponentMarkTriples=timeComponentMarkTriples1000[500:]
goodLogLikelihood=hawkesProcess.getLogLikelihood(secondHalfTimeComponentMarkTriples)
print goodLogLikelihood


badDecayFunction=DecayFunctions.ExponentialDecayFunction(alpha)
badLogLikelihood=GenuineMultivariateHawkesProcess.GenuineMultivariateHawkesProcess(
    np.array([0.021,0.029]),
    np.matrix("[0.61 0.16;0.60 0.06]"),
    [badDecayFunction,badDecayFunction],
    [markDistribution1,markDistribution2]).getLogLikelihood(secondHalfTimeComponentMarkTriples)
print badLogLikelihood