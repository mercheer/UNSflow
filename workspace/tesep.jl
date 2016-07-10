#First attempt to model separation for a cylinder


ntime = 10000
ndiv = 64
dx = (pi)/(n_div+1)
dt = 0.005
ttime = 1.5
cfl = 0.5
tottime = 0
uinf = 1.
c = 1.
re = 1000.

x = zeros(ndiv+2)

#Set sources
for i = 2:n_div+2
    x[i] = [i-1]*dx
end 

uet = zeros(ndiv+2)

for i = 1:n_div+2
    ue[i] = 2*sin(x[i])
    uex[i] = 2*cos(x[i])
end

E = zeros(ndiv+2)
B = zeros(ndiv+2)
del = zeros(ndiv+2)
F = zeros(ndiv+2)

#Initial
for i = 1:ndiv+2
    E[i] = 0.4142
    B[i] = 131.9*0.4142**3-167.32*0.4142**2+76.642*0.4142-11.068
    del[i] = sqrt(B(0:n_div+1)*dt)
    F[i] = 4.8274*0.4142**4-5.9816*0.4142**3+4.0274*0.4142**2+0.23247*0.4142+0.15174
end

ueprev = zeros[ndiv+2]
ueiterprev = zeros[ndiv+2]
delprev = zeros[ndiv+2]
Eprev = zeros[ndiv+2]

for itime = 1, 5
    
    ueprev = ue
    delprev = del
    Eprev = E

    iter_bl = 0
    while
        iter_bl = iter_bl+1
        del = del_prev
        E = E_prev
        ttime = 0.005
        
        #Solve IBL problem
        call calc_bl(ndiv,ue,uex,uet,del,E,ttime)
        
        for i = 1:n_div+2
            m[i] = ue[i]*del[i]
        end
        
        ueiterprev = ue
        for i = 2:ndiv+1
            cauchy[i] = 0.
            for j = 1:ndiv
                adder[j] = (m(j-1)-2*m(j)+m(j+1))/(dx*dx)
                xminus = x[j] - 0.5*dx
                xplus = x[j] + 0.5*dx
                cauchy[i] = cauchy[i] + adder[j]*(-(x[i]-xplus)*(log(abs(x[i]-xplus)))+(xminus-xplus)+(x(i)-xminus)*log(abs(x(i)-xminus)))
            end 
            
            cauchy[i] = (1./(pi*sqrt(re)))*cauchy[i]
            ue[i] = ueprev[i] + cauchy[i] 
        end 
                
        do i = 1, n_div+2
            uex[i] = (ue[i]-ue[i-1])/(x[i]-x(i_div-1))
        end
        uex[0] = 2*uex[1]-uex[2]

        res = sum(abs(ue-ueiterprev))
        if (res<eps) then
            exit
        end if
            
        end
    end 
    
end 
    

calc_bl(n_div,ue,uex,uet,del,E,ttime)

tottime = 0.
ntime = 2000
dx = pi/(dble(n_div+1))
do i_div = 0, n_div+1
    x(i_div) = (i_div)*dx
end do


    unkno(1,:) = del(:)
    unkno(2,:) = del(:)*(E(:)+1)

    do itime=1,ntime
        unkno(1:2,0)=2*unkno(1:2,1)-unkno(1:2,2)
        unkno(1:2,n_div+1)=2*unkno(1:2,n_div)-unkno(1:2,n_div-1)

        !Calculate derived quantities
        do i_div=0,n_div+1
            del(i_div) = unkno(1,i_div)
            E(i_div) = (unkno(2,i_div)/del(i_div))-1
            F(i_div) = 4.8274*E(i_div)**4-5.9816*E(i_div)**3+4.0274*E(i_div)**2+0.232\
            47*E(i_div)+0.15174
            if (E(i_div) < -0.0616) then
                B(i_div) = -225.86*E(i_div)**3-3016.6*E(i_div)**2-208.68*E(i_div)-17.9\
                15;
            elseif (E(i_div) > -0.0395) then
                B(i_div) = 131.9*E(i_div)**3-167.32*E(i_div)**2+76.642*E(i_div)-11.068\
                ;
            else
                B(i_div) = 0.5*(-225.86*E(i_div)**3-3016.6*E(i_div)**2-208.68*E(i_div)\
                                -&
                                &17.915+131.9*E(i_div)**3-167.32*E(i_div)**2+76.642*E(i_div)-11.0\
                                68);
            end if
                if (E(i_div) < -0.0582) then
                    S(i_div) = 451.55*E(i_div)**3+2010*E(i_div)**2+138.96*E(i_div)+11.296;
                elseif (E(i_div) > -0.042) then
                    S(i_div) = -96.739*E(i_div)**3+117.74*E(i_div)**2-46.432*E(i_div)+6.80\
                    74;
                else
                    S(i_div) = 0.5*(451.55*E(i_div)**3+2010*E(i_div)**2+138.96*E(i_div)+11\
                                    .296&
                                    &-96.739*E(i_div)**3+117.74*E(i_div)**2-46.432*E(i_div)+6.8074);
                end if
                    dfde(i_div) = 4*4.8274*E(i_div)**3-3*5.9816*E(i_div)**2+2*4.0274*E(i_div)\
                    +0.23247
                end do



                    !Compute eigenvalues
                    do i_div = 0, n_div+1
                        a_q = 1
                        b_q = -ue(i_div)*(dfde(i_div)-1)
                        c_q = ue(i_div)*ue(i_div)*(E(i_div)*dfde(i_div)-F(i_div))
                        lamb1(i_div) = (-b_q+sqrt(b_q*b_q-4*a_q*c_q))/(2*a_q)
                        lamb2(i_div) = (-b_q-sqrt(b_q*b_q-4*a_q*c_q))/(2*a_q)

                        !Always have lamb1 > lamb2
                        if (lamb2(i_div)>lamb1(i_div)) then
                            temp  = lamb2(i_div)
                            lamb2(i_div) = lamb1(i_div)
                            lamb1(i_div) = temp
                        end if
                        end do


                            !Compute timestep
                            dt=10000
                            do i_div = 1,n_div
                                dti=cfl*(dx/(abs(lamb1(i_div))+abs(lamb2(i_div))))
                                if (dti<dt) then
                                    dt=dti
                                end if
                                end do
                                    
                
