%% For each time step the code loads the flight path(given U,V,W ,p,q,r) and initial geometry and converts it into (X(X0,Y0,Z0),Y(X0,Y0,Z0),W(X0,Y0,Z0))

function [X,Y,Z] = func_kinematics_pos(t_now,X_0,Y_0,Z_0,Quat)

addpath([pwd,'/data']);

load([pwd,'/data/data_kinematics.mat']);

            
% the transformation matrix converts the vector in body reference system
% (DISTANCE OF POINT P FROM O at time 0) in a vector
% in the IRF at time t due to the attitude of BRF --> IRF
            
r = transpose(quat2dcm(Quat))*[X_0-X_D;Y_0-Y_D;Z_0-Z_D];

x = r(1);
y = r(2);
z = r(3);

%% Origin motion

% Integration for ORIGIN motion

vTime2 = linspace(0,t_now,100);
Uint = U(vTime2);
Vint = V(vTime2);
Wint = W(vTime2);

X_pst = trapz(vTime2,Uint);
Y_pst = trapz(vTime2,Vint);
Z_pst = trapz(vTime2,Wint);

% To add terms that come from distance from axes of rotations velocities p,q,r

X = x + X_pst;
Y = y + Y_pst;
Z = z + Z_pst;


end