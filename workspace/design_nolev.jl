workspace()
include("../src/UNSflow.jl")
using UNSflow

#Base motion is an Eldredge pitch-up of SD7003 at Re=20k with amp=30 deg, K=0.2, smoothing =0.8
#Plunge is superimposed so that LEV formation is prevented
#Plunge is an integrated Eldredge pitch-up (Derivative of plunge is an Eldrege pitch-up)

alpha_amp = 30*pi/180
K = 0.2
a = 0.8

alphadef = EldUpDef(alpha_amp, K, a)
design_solve(alphadef)
