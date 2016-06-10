function [u,v,w] = func_vortonind(X_CP,Y_CP,Z_CP,vorton,sig)

%% PAY ATTENTION TO SINGULAR OR NON-SINGULAR VORTONS

% sig = 0.01; % SINGULAR VORTON

%% VORTON AND POINT 

R_CP = [X_CP;Y_CP;Z_CP];
R_VO = [vorton.X;vorton.Y;vorton.Z];

ALFA_VO = [vorton.A_X;vorton.A_Y;vorton.A_Z];

%% INDUCTION

DIST = R_CP - R_VO;

K_sig = -1/(4*pi)*(norm(DIST)^2+5/2*sig^2)/(norm(DIST)^2+sig^2)^(5/2);

K_DIST = K_sig.*DIST;

INDUCTION = cross(K_DIST,ALFA_VO);


u = INDUCTION(1);
v = INDUCTION(2);
w = INDUCTION(3);


end
