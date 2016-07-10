workspace()
include("../src/UNSflow.jl")
using UNSflow

init_alpha = 0.
init_alphadot = 0.
init_h = 0.
init_hdot = 0.

pvt = 0.35

lespcrit = [0.11;]


#Provide initial conditions
kinem = KinemPar2DOF(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)

kinem.alpha = kinem.alpha_pr = kinem.alpha_pr2 = init_alpha*pi/180
kinem.h = kinem.h_pr = kinem.h_pr2 = init_h*pi/180

kinem.alphadot = kinem.alphadot_pr = kinem.alphadot_pr2 = init_alphadot*pi/180
kinem.hdot = kinem.hdot_pr = kinem.hdot_pr2 = init_hdot*pi/180

kinem.alphaddot = kinem.alphaddot_pr = kinem.alphaddot_pr2 = 0.
kinem.hddot = kinem.hddot_pr = kinem.hddot_pr2 = 0.

#Provide structural parameters
strpar = TwoDOFPar(0,0,0,0,0,0,0,0,0,0,0)
strpar.x_alpha = 0.05
strpar.r_alpha = 0.5
strpar.w_alpha = 1.0
strpar.w_h = 1.0
strpar.w_alphadot = 0.0
strpar.w_hdot = 0.0
strpar.cubic_alpha_1 = 1.0
strpar.cubic_alpha_3 = 0.0
strpar.cubic_h_1 = 1.0
strpar.cubic_h_3 = 0.0

surf = TwoDSurf_2DOF(1., 1., "sd7003_fine.dat", pvt, 70, 35, strpar, kinem, lespcrit)

curfield = TwoDFlowField()

dt = 0.015
t_tot = 10.
nsteps =round(Int,t_tot/dt)+1

t = 0.

outfile = open("results.dat", "w")

for istep = 1:nsteps
    t = t + dt
    (surf, curfield, cl, cd, cm) = ldvmstep(surf, curfield, dt)

    m11 = (2./surf.c)
    m12 = -surf.strpar.x_alpha*cos(surf.kinem.alpha)
    m21 = -2*surf.strpar.x_alpha*cos(surf.kinem.alpha)/surf.c
    m22=surf.strpar.r_alpha*surf.strpar.r_alpha
    
    R1 = ((4*surf.strpar.kappa*surf.uref*surf.uref*cl)/(pi*surf.c*surf.c)) - (2*surf.strpar.w_h*surf.strpar.w_h*((surf.strpar.cubic_h_1*surf.kinem.h) + (surf.strpar.cubic_h_3*(surf.kinem.h^3)))/surf.c) - (surf.strpar.x_alpha*sin(surf.kinem.alpha)*surf.kinem.alphadot*surf.kinem.alphadot)
    
    R2 = (8*surf.strpar.kappa*surf.uref*surf.uref*cm)/(pi*surf.c*surf.c) - (surf.strpar.w_alpha*surf.strpar.w_alpha*surf.strpar.r_alpha*surf.strpar.r_alpha*((surf.strpar.cubic_alpha_1*surf.kinem.alpha) + (surf.strpar.cubic_alpha_3*(surf.kinem.alpha^3))))
    
    
    #Update all the previous values
    surf.kinem.h_pr3 = surf.kinem.h_pr2
    surf.kinem.alpha_pr3 = surf.kinem.alpha_pr2
    
    surf.kinem.h_pr2 = surf.kinem.h_pr
    surf.kinem.alpha_pr2 = surf.kinem.alpha_pr
    
    surf.kinem.h_pr = surf.kinem.h
    surf.kinem.alpha_pr = surf.kinem.alpha
        
    surf.kinem.hdot_pr3 = surf.kinem.hdot_pr2
    surf.kinem.alphadot_pr3 = surf.kinem.alphadot_pr2
    
    surf.kinem.hdot_pr2 = surf.kinem.hdot_pr
    surf.kinem.alphadot_pr2 = surf.kinem.alphadot_pr
    
    surf.kinem.hdot_pr = surf.kinem.hdot
    surf.kinem.alphadot_pr = surf.kinem.alphadot
    
    surf.kinem.hddot_pr3 = surf.kinem.hddot_pr2
    surf.kinem.alphaddot_pr3 = surf.kinem.alphaddot_pr2

    surf.kinem.hddot_pr2 = surf.kinem.hddot_pr
    surf.kinem.alphaddot_pr2 = surf.kinem.alphaddot_pr

    surf.kinem.hddot_pr = surf.kinem.hddot
    surf.kinem.alphaddot_pr = surf.kinem.alphaddot
    
    #Solve structural equations
    surf.kinem.hddot = (1/(m11*m22 - m21*m12))*(m22*R1 - m12*R2)
    surf.kinem.alphaddot = (1/(m11*m22 - m21*m12))*(-m21*R1 + m11*R2)
    #Continuation/Update kinematics through Adam Bashforth
    surf.kinem.alphadot = surf.kinem.alphadot_pr + ((dt/12.)*(23.*surf.kinem.alphaddot_pr - 16.*surf.kinem.alphaddot_pr2 + 5.*surf.kinem.alphaddot_pr3))
    surf.kinem.hdot = surf.kinem.hdot_pr + ((dt/12.)*(23.*surf.kinem.hddot_pr - 16.*surf.kinem.hddot_pr2 + 5*surf.kinem.hddot_pr3))
    surf.kinem.alpha = surf.kinem.alpha_pr + ((dt/12.)*(23.*surf.kinem.alphadot_pr - 16.*surf.kinem.alphadot_pr2 + 5.*surf.kinem.alphadot_pr3))
    surf.kinem.h = surf.kinem.h_pr + ((dt/12.)*(23.*surf.kinem.hdot_pr - 16.*surf.kinem.hdot_pr2 + 5.*surf.kinem.hdot_pr3))
    
    
    write(outfile, join((t, surf.kinem.alpha, surf.kinem.h, surf.kinem.u, surf.a0[1], cl, cd, cm)," "), "\n")
    
end
close(outfile)
    
