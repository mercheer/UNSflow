using JFVM

Lx = pi
Nx = 64

mesh = createMesh1D(Nx, Lx)

x = mesh.cellcenters.x

ue = createCellVariable(mesh, 2*sin(x))
uex = createCellVariable(mesh, 2*cos(x))
uet = createCellVariable(mesh, 0.)
dt = 0.005
#Inital condition
E = createCellVariable(mesh, 0.4142)
B = createCellVariable(mesh, 131.9*E.value.^3 - 167.32*E.value.^2 + 76.642*E.value - 11.068)
del = createCellVariable(mesh, sqrt(B.value*dt))
F = createCellVariable(mesh, 4.8274*E.value.^4 - 5.9816*E.value.^3 + 4.0274*E.value.^2 + 0.23247*E.value + 0.15174)
S = createCellVariable(mesh,0) #Not required for initial condistion, just initialised
dfde = createCellVariable(mesh, 4*4.8274*E.value.^3 - 3*5.9816*E.value.^2 + 2*4.0274*E.value + 0.23247)
r1 = createCellVariable(mesh, 0.) #Just initialised
r2 = createCellVariable(mesh, 0.) #Just initialised

# BC1 = createBC(mesh) #Newmann by default

# BC1.left.a[:] = 0 #Dirichlet

# BC1.left.b[:] = 1

# BC1.left.c[:] = ue.value[1]*del.value[1]*E.value[1] 

# BC2 = createBC(mesh)

# BC2.left.a[:] = 0 #Dirichlet

# BC2.left.b[:] = 1

# BC2.left.c[:] = ue.value[1]*del.value[1]*F.value[1]

phi1_old = createCellVariable(mesh, del.value[1] , BC1)
phi2_old = createCellVariable(mesh, del.value[1]*(E.value[1]+1), BC2)

phi1 = phi1_old
phi2 = phi2_old

alfa = createCellVariable(mesh,1.0)


FL = fluxLimiter("MUSCL")

final_t = 1.5
t = 0.

while t < final_t
    t = t + dt
    
    println(t)

    #Update variables
    del.value = phi1_old.value
    E.value = phi2_old.value./phi1_old.value - 1


    
    F.value = 4.8274*E.value.^4 - 5.9816*E.value.^3 + 4.0274*E.value.^2 + 0.23247*E.value + 0.15174
    for i = 1:Nx+2
        if (E.value[i] < -0.0616) 
            B.value[i] = -225.86*E.value[i]^3 - 3016.6*E.value[i]^2 - 208.68*E.value[i] - 17.915
        elseif (E.value[i] > -0.0395)
            B.value[i] = 131.9*E.value[i]^3 - 167.32*E.value[i]^2 + 76.642*E.value[i] - 11.068
        else
            B.value[i] = 0.5*(-225.86*E.value[i]^3 - 3016.6*E.value[i]^2 - 208.68*E.value[i] - 17.915 + 131.9*E.value[i]^3 - 167.32*E.value[i]^2 + 76.642*E.value[i] - 11.068)
        end 
        if (E.value[i] < -0.0582)
            S.value[i] = 451.55*E.value[i]^3 + 2010*E.value[i]^2 + 138.96*E.value[i] + 11.296;
        elseif (E.value[i] > -0.042) then
            S.value[i] = -96.739*E.value[i]^3 + 117.74*E.value[i]^2 - 46.432*E.value[i] + 6.8074
        else
            S.value[i] = 0.5*(451.55*E.value[i]^3 + 2010*E.value[i]^2 + 138.96*E.value[i] + 11.296 - 96.739*E.value[i]^3 + 117.74*E.value[i]^2 - 46.432*E.value[i] + 6.8074)
        end 
    end
    dfde.value = 4*4.8274*E.value.^3 - 3*5.9816*E.value.^2 + 2*4.0274*E.value + 0.23247

    (Mbc1, RHSbc1) = boundaryConditionTerm(BC1)
    (Mbc2, RHSbc2) = boundaryConditionTerm(BC2)
        
    (Mt1, RHSt1) = transientTerm(phi1_old, dt, alfa);
    (Mt2, RHSt2) = transientTerm(phi2_old, dt, alfa);
       
    a11 = createFaceVariable(mesh, -1)
    a12 = createFaceVariable(mesh, 1)
    a21 = linearMean(F - (E + 1).*dfde)
    a22 = linearMean(dfde)
    
    (Mconv11, RHSc11) = convectionTvdTerm(a11, phi1_old, FL)
    (Mconv12, RHSc12) = convectionTvdTerm(a12, phi2_old, FL)
    (Mconv21, RHSc21) = convectionTvdTerm(a21, phi1_old, FL)
    (Mconv22, RHSc22) = convectionTvdTerm(a22, phi2_old, FL)
    
    r1.value = B.value./(2*del.value) - del.value.*uet.value./ue.value - (E.value + 1).*del.value.*uex.value - E.value.*del.value.*uex.value
    r2.value = S.value./del.value - 2*E.value.*del.value.*uet.value./ue.value - 2*F.value.*del.value.*uex.value - F.value.*del.value.*uex.value    
    RHSs1 = constantSourceTerm(r1)
    RHSs2 = constantSourceTerm(r2)


    Muw = [Mt1 + Mconv11 + Mbc1 Mconv12; Mconv21 Mt2 + Mconv22 + Mbc2]
    RHS = [RHSt1 + RHSbc1 + RHSs1 + RHSc11 + RHSc12; RHSt2 + RHSbc2 + RHSs2 + RHSc21 + RHSc22]
    X = Muw\RHS

    phi1.value = X[1:Nx+2]
    phi2.value = X[(Nx+2+1):end]
    
    phi1_old = phi1
    phi2_old = phi2

end


crit = diff(del.value[2:Nx+1])./diff(x)
