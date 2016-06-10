function [u,v,w] = func_vortexl(X_CP,Y_CP,Z_CP,X_1,Y_1,Z_1,X_2,Y_2,Z_2,GAMMA)
% To calculate induced velocity of a vortex line element with strength per
% unit length = 1;The vortex line point in the direction (2-1);see Katz &
% Plotkin Low Speed Aerodynamics,VLM programme 13

R_CUT = 1*10^(-10);

%%%%%% CHANGED TO -7 FROM -5 IN 24/05 FOR SMALL WINGS
%%%%%% CHANGED TO 0.15% OF MINIMUM SIZE IN THE PANELS


% Calculation of R1 * R2

R1R2X = (Y_CP-Y_1)*(Z_CP-Z_2)-(Z_CP-Z_1)*(Y_CP-Y_2);
R1R2Y = -((X_CP-X_1)*(Z_CP-Z_2)-(Z_CP-Z_1)*(X_CP-X_2));
R1R2Z = (X_CP-X_1)*(Y_CP-Y_2)-(Y_CP-Y_1)*(X_CP-X_2);

% Calculation of (R1*R2)^2

SQUARE = (R1R2X)^2 + (R1R2Y)^2 + (R1R2Z)^2;

% Calculation of R0(R1/R(R1)-R2/R(R2))

R1 = sqrt((X_CP-X_1)^2+(Y_CP-Y_1)^2+(Z_CP-Z_1)^2);
R2 = sqrt((X_CP-X_2)^2+(Y_CP-Y_2)^2+(Z_CP-Z_2)^2);

if (R1<R_CUT) || (R2<R_CUT) || (SQUARE < R_CUT)
    u = 0;
    v = 0;
    w = 0;
else
    R0R1 = (X_2 - X_1)*(X_CP-X_1)+(Y_2 - Y_1)*(Y_CP-Y_1)+(Z_2 - Z_1)*(Z_CP-Z_1);
    R0R2 = (X_2 - X_1)*(X_CP-X_2)+(Y_2 - Y_1)*(Y_CP-Y_2)+(Z_2 - Z_1)*(Z_CP-Z_2);
    COEF = GAMMA/(4*pi*SQUARE)*(R0R1/R1-R0R2/R2);
    
    u = R1R2X *COEF;
    v = R1R2Y *COEF;
    w = R1R2Z *COEF;
    
end

