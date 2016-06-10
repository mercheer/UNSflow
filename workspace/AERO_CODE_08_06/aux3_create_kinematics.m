%% Auxiliary script for kinematics definition

clc;clear all;close all;

%% ADDPATH

addpath([pwd,'/data']);

t_fin = 8*pi;

%Initial Attitude

theta0 = convang(0,'deg','rad'); 
phi0   = convang(0,'deg','rad'); 
psi0   = convang(0,'deg','rad');

%Anonymous Function for time histories of ang. velocity

% Pitch Maneuver (try to copy the paper)
%
% t_pitch = 1;
% alfa_max = deg2rad(45);
% omega_0 = alfa_max/t_pitch;
%   
% 
% a = 11;
%   t1=1;
%    t2=t1+(t_pitch);
%    maxtime=t2+1;
% 
%    t_fin = maxtime;
% 
% alpha = @(t) ...
% ((omega_0/(2*a))*log(cosh(a*(t-t1))/cosh(a*(t-t2))))+(alfa_max/2);
% 
% K = omega_0/2;

% a = 0.925;
% 
% K = pi/8;
%   amp= 45;
% 
% a = pi*pi*K/(2*(amp*pi/180)*(1-a));
%   t1=1;
%    t2=t1+((amp*pi/180)/(2*K));
%    maxtime=t2+1;
% 
%    t_fin = maxtime;
% 
%    
% q =@(t) ...
%     ((K)*((tanh(a*(t-t1))-tanh(a*(t-t2)))));


%% U0 V0 W0 declaration in m/s


U0 = -1.0;
V0 = 0.0;
W0 = tan(deg2rad(20))*U0;

Q_INF = norm(U0);

RED_FREQ = 0.50;
AMP = 0.1;

W_MAX = AMP * RED_FREQ * 2 * Q_INF;



U=@(t) ...
    interp1( ...
    [0, t_fin],...
    [U0,        U0],...
    t,'linear'...
    );

V=@(t) ...
    interp1( ...
    [0, t_fin/10, t_fin*2/10, t_fin*3/10, t_fin*4/10, t_fin*5/10, t_fin*6/10,   t_fin*7/10, t_fin*8/10,t_fin*9/10,      t_fin],...
    [V0,         V0,       0.2*V0,          0.3*V0,         0.4*V0,         0.8*V0,        V0,           V0,       V0,        V0,         V0],...
    t,'pchip'...
    );

W=@(t)  ...
   interp1( ...
   [0, t_fin],...
   [W0,   W0],...
   t,'linear'...
   );
% 
%  W=@(t)  ...
%      W_MAX*cos(RED_FREQ*2*Q_INF*t);

% OMEGA_BEND = 0.8889*2*pi;
% 
% Z_BEND = @(t,x) ...
%     (-0.0032*x^4+1.4*10^(-6)*x^3+0.1*x^2-1.8*10^(-5)*x-0.2)*cos(OMEGA_BEND*t);
% 
% W_BEND = @(t,x) ...
%     -OMEGA_BEND*(-0.0032*x^4+1.4*10^(-6)*x^3+0.1*x^2-1.8*10^(-5)*x-0.2)*sin(OMEGA_BEND*t);

%% RED_FREQ AND AMP ARE VALID ONLY IF CHORD == 1


AMP_PHI = deg2rad(0);
AMP_TWIST = deg2rad(0);

% p q r declaration

p_max = convangvel(0.0,'deg/s','rad/s'); %Converts units of measurements of angular velocity
q_max = convangvel(0.0,'deg/s','rad/s');
r_max = convangvel(0.0,'deg/s','rad/s');



% p=@(t)  ...
%       -AMP_PHI*RED_FREQ*2*sin(RED_FREQ*2*Q_INF*t);    
% 
% p=@(t) ...
%     interp1(...
%     [0, t_fin/10, t_fin*2/10, t_fin*3/10, t_fin*4/10, t_fin*5/10, t_fin*6/10,   t_fin*7/10, t_fin*8/10,t_fin*9/10,      t_fin],...
%     [p_max,p_max/2,0,0,0,0,0,0,0,0,0],...
%     t,'linear'...
%     );
% q=@(t) ...
%     interp1(...
%     [0,  0.125,      0.375,   0.5, 0.625, 0.875,    1],...
%     [1.59, 0.0,        0.0, -1.59,   0.0,   0.0, 1.59],...
%     t,'phcip'...
%     );
% r=@(t) ...
%     interp1(...
%     [0,    0.05, 0.45, 0.5,  0.55, 0.95,   1],...
%     [0.0, -0.83,-0.83, 0.0, +0.83,+0.83, 0.0],...
%     t,'phcip'...
%     );
  
 p = @(t) ...
     0;
 q = @(t) ...
     0;
 r = @(t) ...
     0;

%% ROTATION AXES DISTANCES

X_D = 0.0;
Y_D = 0.0;
Z_D = 0.0;

%% WING-TIP TWIST TIME-LAW

% twist = @(t)...
%     -AMP_TWIST*sin(RED_FREQ*2*Q_INF*t);

 twist = @(t) ...
     0;

 
%%
save([pwd,'/data/data_kinematics.mat']);