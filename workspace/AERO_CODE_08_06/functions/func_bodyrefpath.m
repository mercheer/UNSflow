function [s_x,s_y,s_z] = func_bodyrefpath(t,dt,x,y,z,Q_inf)

addpath([pwd,'/functions']);


       [u1_tot1,v1_tot1,w1_tot1] = func_kinematics(t,...
           x,y,...
           z);

       u1_tot = u1_tot1;
       v1_tot = v1_tot1;
       w1_tot = w1_tot1;
       
       u1 = u1_tot-Q_inf(1);
       v1 = v1_tot-Q_inf(2);
       w1 = w1_tot-Q_inf(3);
       s_x = u1*dt;
       s_y = v1*dt;
       s_z = w1*dt;
end