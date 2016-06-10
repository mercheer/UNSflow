%% For each time step the code loads the flight path(given U,V,W ,p,q,r) and converts it into kinematic velocity of undisturbed fluid in BRF

% NEED THE POSITION OF THE POINT AT t_now IN THE BRF (x,y,z) AND LOAD THE
% KINEMATICS OF THE WING (U,V,W,p,q,r)

function [u,v,w] = func_kinematics_vel(t_now,x,y,z,Quat)

addpath([pwd,'/data']);

load([pwd,'/data/data_kinematics.mat']);

% Conversion to body-fixed frame of the speed of the wing
            
dvdt = quat2dcm(Quat)*[-U(t_now);-V(t_now);-W(t_now)];

u3 = dvdt(1);
v3 = dvdt(2);
w3 = dvdt(3);


% To add terms that come from distance from axes of rotations velocities p,q,r

u = u3 - q(t_now)*(z-Z_D) + r(t_now)*(y-Y_D);
v = v3 - r(t_now)*(x-X_D) + p(t_now)*(z-Z_D);
w = w3 - p(t_now)*(y-Y_D) + q(t_now)*(x-X_D); 


end