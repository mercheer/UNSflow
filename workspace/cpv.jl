# Evaluates the CPV integral of a user supplied function phi(x) over (a,b) at the point x0 in (a,b) by a n-point Double Epxponential (DE) scheme.

#The integal defintion is integral with respect to x over (a,b ) of  phi(x) /(x-x0).

# Double Exponential transformation between x & u is given by
#x= (b+a)/2 + (b-a)/2 tanh(pi/2 sinh(u))
#x0= (b+a)/2 + (b-a)/2 tanh(pi/2 sinh(u0))
#x, x0 lie in (a,b)
# u, u0 lie in (-inf,inf)
# Credits go to Euler, Mclaurin, C.Schwartz, Iri, Moriguti, Takasawa, Takahasi and Mori
#Loosely converted to Julia code by Kiran Ramesh 28 Jun 2016

function cpv(phi, a :: Real, b :: Real, x0 :: Real, n :: Int64)
    # Given x0, invert to get u0
    arg = (2*x0 - (b + a))/(b - a)
    u0 = asinh((2/pi)*atanh(arg))
    zmax = 3.5
    h = 2*zmax/n
    # Sampling points are chosen such that u0 lies between two successive quadrature points; Integral multiple of "h" is sliced off from (-zmax,zmax) and trapezoidal summation is done over over that interval, (umin,umax1)
    
    (umax1,hu,umin,hl) = inter(zmax,u0,h)
    z = umin - h
    sum=0
    while z <= umax1 - h
        (x,dvdu) = de(a,b,z)
        arg = (phi(x))/(x - x0)
        der = dvdu*2*(x - a)*(b - x)/(b - a)
	sum = sum + arg*der
        z=z+h
      end
    sum = h*sum

    # Weight has to be (h/2) for z=umin & z=umax1
    z=umin
    (x,dvdu) = de(a,b,z)
    arg = (phi(x))/(x - x0)
    der = dvdu*2*(x - a)*(b - x)/(b - a)
    sum = sum + der*arg*h/2
  
    z = umax1
    (x,dvdu) = de(a, b, z)
    arg = (phi(x))/(x - x0)
    der = dvdu*2*(x - a)*(b - x)/(b - a)
    sum = sum + der*arg*h/2
    
    # Add the contribution from the right fractional step by taking
    # z as the midpoint of  -zmax & umin;
    z = (umin - zmax)/2
    (x,dvdu) = de(a,b,z)
    arg = (phi(x))/(x - x0)
    der = dvdu*2*(x - a)*(b - x)/(b - a)
    sum = sum + hl*der*arg
    # Add the contribution from the right fractional step by taking
    # z as the midpoint of  zmax & umax1;
    z = (zmax + umax1)/2
    (x,dvdu) = de(a,b,z)
    arg = (phi(x))/(x - x0)
    der = dvdu*2*(x - a)*(b - x)/(b - a)
    sum = sum + hu*der*arg
    return sum
end

function de(a :: Real,b :: Real,u :: Real)
    #Given the interval (a,b) and u, the new variable  we get the old variable "x" & the transformation derivative x= (a+b)/2 +(b-a)/2 tanh(pi/2 sinh(u)); der= dx/du= (dx/dv)*(dv/du)
    # Note that the transformation derivative is arrived at as follows.
    # x= (b+a)/2 + (b-a)/2 tanh(v); v= (pi/2) sinh(u)
    # dx/du= dx/dv * dv/du; dx/dv= 2(b-x)(x-a)/(b-a); dv/du = (pi/2) cosh(u)
    # The term dvdu that is transferred here is just dv/du. It is multiplied by 
    # dx/dv= 2(b-x)(x-a)/(b-a) in the CPV routine 
    # 6-Feb-07
    # Credits go to Euler, Mclaurin, C.Schwartz, Iri, Moriguti, Takasawa, Takahasi and Mori
    #Loosely converted to Julia code by Kiran Ramesh 28 Jun 2016
    
    y = exp(u)
    yi = exp(-u)
    x1 = 0.5*(y + yi)
    x2 = 0.5*(y - yi)
    x3 = exp(pi*0.5*x2)
    x4 = exp(-pi*0.5*x2)
    x = (b*x3+a*x4)/(x3+x4)
    dvdu = 0.5*pi*x1
    return x, dvdu
end

function inter(umax :: Real,u0 :: Real,h :: Real)
    #umax =upper de integration limit; -umax=lower de integration limit
    # h=step size for DE; x0=demap(u0); x0=where you need the cpv integral in the old interval
    #Between (u0+ h/2) and umax find out the integer(j) multiples of "h"
    #The remainder is the fractional step size "hu"
    #Integrate in u from (u0+h/2) to umax1;
#Add the contribution for the fractional width hu
    # This is for the Left side of u0; For the Right side of u0
    # Find the integer multiples of h between -umax and (u0+h/2);
    #The remainder is the fraction step hl;
    # Integrate from umin to umax1 & just add
    # the contributions from the fractional widths on either side
    
    z = umax - u0 - (h/2)
    hu = rem(z,h)
    umax1 = umax - hu
    z = u0 + (h/2) + umax
    hl = rem(z,h)
    umin = -umax + hl
    return umax1, hu, umin, hl
end
