function [vr1,vr2,vr3] = func_rodrot(ang,e1,e2,e3,v1,v2,v3)

vc = cross([e1,e2,e3],[v1,v2,v3]);
dot_p = dot([e1,e2,e3],[v1,v2,v3]);

cr1=vc(1);
cr2=vc(2);
cr3=vc(3);

vr1 = v1*cos(ang)+cr1*sin(ang) + e1*dot_p*(1-cos(ang));
vr2 = v2*cos(ang)+cr2*sin(ang) + e2*dot_p*(1-cos(ang)); 
vr3 = v3*cos(ang)+cr3*sin(ang) + e3*dot_p*(1-cos(ang)); 

end