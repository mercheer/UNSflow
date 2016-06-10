function Quat = func_upd_att(t_now,Quat_0)

addpath([pwd,'/data']);

load([pwd,'/data/data_kinematics.mat']);


%Kinematic Equations

dQuatdt=@(t,Q) ...
    0.5*[   0, -p(t), -q(t), -r(t);
         p(t),     0,  r(t), -q(t);
         q(t), -r(t),     0,  p(t);
         r(t),  q(t), -p(t),     0] * Q;
     
%Solution of quaternion component evolution equations

options = odeset( ...
    'RelTol', 1e-3, ...
    'AbsTol', 1e-3*ones(1,4) ...
    );

if t_now == 0 
    t_now = 0.00000000000001;
end

% Runge-Kutta Integration
[vTime, vQuat] = ode45(dQuatdt, [0 t_now], Quat_0, options);


% Time interpolation function for known quaternion history
Quat = @(t) ...
    [interp1(vTime,vQuat(:,1),t), ...
        interp1(vTime,vQuat(:,2),t), ...
            interp1(vTime,vQuat(:,3),t), ...
                interp1(vTime,vQuat(:,4),t)];


Quat = Quat(t_now);


end