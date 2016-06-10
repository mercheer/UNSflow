function [u,v,w] = func_vortonstre(vorton1,vorton2,sig)

%% PAY ATTENTION TO SINGULAR OR NON-SINGULAR VORTONS

% sig = 0.01; % SINGULAR VORTON

%% VORTON AND POINT 

R_P = [vorton1.X;vorton1.Y;vorton1.Z];
R_Q = [vorton2.X;vorton2.Y;vorton2.Z];

ALFA_P = [vorton1.A_X;vorton1.A_Y;vorton1.A_Z];
ALFA_Q = [vorton2.A_X;vorton2.A_Y;vorton2.A_Z];


%% INDUCTION

DIST = norm(R_P- R_Q);

vec = 1/(4*pi)*(((DIST^2+5/2*sig^2)/(DIST^2+sig^2)^(5/2))*cross(ALFA_P,ALFA_Q)+3*...
    ((DIST^2+7/2*sig^2)/(DIST^2+sig^2)^(7/2))*(dot(ALFA_P,R_P-R_Q)*cross(R_P-R_Q,ALFA_Q)));


u = vec(1);
v = vec(2);
w = vec(3);

