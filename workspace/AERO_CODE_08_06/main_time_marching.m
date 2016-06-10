 %% Time marching solver of unsteady case

clearvars -except i_SS alfa_SS C_L_steady_state C_Di_steady_state C_L_steady_state_conv C_L_central_panel TIME_FIN;

global  NW_wake_vortex_rings LE_NW_wake_vortex_rings FW_wake_vortons LE_FW_wake_vortons...
       t_fin n_rows n_cols ...
       rho eta s chord

tot_time = tic;

%% Add paths

addpath([pwd,'/geometry']);
addpath([pwd,'/functions']);
addpath([pwd,'/data']);

%% load AIM,panels,airfoil

load('data_vortex_lattice.mat');
load('data_airfoil.mat');
load('data_refs.mat');
load('data_kinematics.mat');



rho = 0.09; % density kg/m^3


%% Flags (0 INACTIVE,1 ACTIVE)

flag_RU        = 0; % Flag to activate free wake evolution function,time consuming.

flag_MBD_IO    = 0; % Flag to activate MBDyn communication,do not activate if you are not using MBDyn.

flag_VORTONS   = 1; % Flag to transform wake vortex rings panels into vortons.

flag_LEVS      = 0; % Flag to activate wake panels shedding from the LEADING EDGE.

flag_FRSTR     = 0; % Flag to change FREE-STREAM reference velocity used in force calculation function,do not activate for simple cases(still issues).

flag_WAKEDIST  = 0;  % Flag to change the distance of the trailing segment of the first row of wake panels, LENGTH WAKE PANEL = cw * FREE-STREAM VELOCITY * TIME STEP.
cw = 0.50;

flag_BENDMOT   = 0; % Flag to activate bending motion of the wing,do not activate.

flag_FRAME     = 0; % Flag to activate frame recording.

%% Time-step calculation


if flag_FRSTR == 1

t1 = linspace(0,t_fin,100);
Q_inf_med = max(r(t1))*b_ref/2;



else

   t1 = zeros(1,100);
   Q_inf1 = cell(1,100);
   Q_inf_med = 0;

   for j = 1:101
     if j == 1
      t1(j) = 0;
      else
      if j == 101
      t1(j) = t_fin;
      else
      t1(j) = t_fin/101*(j+1);
      end
     end

    [u_inf1,v_inf1,w_inf1] = func_kinematics_vel(t1(j),0,0,0,angle2quat(0,0,0));
    Q_inf1{j} = [u_inf1,v_inf1,w_inf1];
    Q_inf_med = Q_inf_med+norm(Q_inf1{j});

   end

Q_inf_med = Q_inf_med/101;


end

if flag_WAKEDIST == 1

    N_t = 160;
    dt = t_fin/N_t;

else
    k   = 1/n_rows;
    dt  = k*c_ref/Q_inf_med;
    N_t = round(t_fin/dt);
end

sig = 1.3*Q_inf_med*dt/2;

disp(['Time step is : ' num2str(dt)]);
disp(['Sigma is : ' num2str(sig)]);

%% Display simulation conditions

fprintf('\n');
disp(['Number of total time steps is : ' num2str(N_t)]);
disp(['Number of rows is : ' num2str(n_rows)]);
disp(['Number of columns is : ' num2str(n_cols*2)]);
fprintf('\n');


%% Pre-allocation and storage


C_L_v          = zeros(1,N_t);
C_Di_v         = zeros(1,N_t);
C_Di_v_st      = zeros(1,N_t);
C_Di_v_unst    = zeros(1,N_t);


timerec_SOLVER = zeros(1,N_t);
step_time      = zeros(1,N_t);
est_time       = zeros(1,N_t);
time           = zeros(1,N_t);
adim_time      = zeros(1,N_t);

C_L_mat  = cell(1,N_t);
C_Di_mat = cell(1,N_t);
F_mat    = cell(1,N_t);
L_mat    = cell(1,N_t);

gamma_ad_C     = cell(1,N_t);
spanload       = cell(1,N_t);
adim_spanload  = cell(1,N_t);
KATZ_VERIFY_C  = cell(1,N_t);
gamma_mat_CELL = cell(1,N_t);

NW_wake_vortex_rings    = cell(N_t-1,n_cols);
gamma_w_mat             = zeros(N_t-1,n_cols);
FW_wake_vortons         = cell(N_t-1,n_cols);
LE_NW_wake_vortex_rings = cell(N_t-1,n_cols);
LE_gamma_w_mat          = zeros(N_t-1,n_cols);
LE_FW_wake_vortons      = cell(N_t-1,n_cols);


%% Save BFR vortex lattice positions at istant 0

vortices_mat_0    = vortices_mat;
vortices_mat_CELL = cell(1,N_t);
Quat_CELL         = cell(1,N_t);
Quat_CELL{1}      = angle2quat (psi0, theta0, phi0);  %Converts Euler angles into Quaternion


%% Time loop with wake shedding

for i_w = 1:N_t


%% Wait for the output.dat from MBDyn

if flag_MBD_IO == 1
waiting_time = tic;
 if i_w > 2

  pause(1)
  while (exist('INI_OUTPUT.dat','file')) == 0
  pause(1)
  end
  normal()
 end
waiting_time = toc(waiting_time);

else
end

%% Time

t = t_fin*(i_w)/N_t;

%% Velocity of the origin of BODY FIXED REF. SYSTEM in INERTIAL REFERENCE SYSTEM

if flag_FRSTR == 1
    Q_inf               = [Q_inf_med,0,0];
else
    [u_inf,v_inf,w_inf] = func_kinematics_vel(t,0,0,0,angle2quat(0,0,0));
    Q_inf               = [u_inf,v_inf,w_inf];
end

if flag_WAKEDIST == 1
adim_time(i_w) = (Q_inf_med*t)/c_ref;
else
adim_time(i_w) = k * i_w;
end

time(i_w)      = t;

disp(['Time step number : ' num2str(i_w) ' non-dimensional time is : ' num2str(adim_time(i_w)) ]);
disp(['Dimensional time is : ' num2str(time(i_w)) ]);


%% Updating wing position in inertial reference frame

if i_w == 1

time1_UPDPOS = tic;

    Quat_CELL{i_w}         = func_upd_att(t,Quat_CELL{1});
    vortices_mat           = func_upd_pos(t,vortices_mat_0,Quat_CELL{i_w},flag_BENDMOT);
    vortices_mat_CELL{i_w} = vortices_mat;

time1_UPDPOS = toc(time1_UPDPOS)

else

time1_UPDPOS = tic;
    Quat_CELL{i_w} = func_upd_att(t,Quat_CELL{1});
    vortices_mat = func_upd_pos(t,vortices_mat_0,Quat_CELL{i_w},flag_BENDMOT);
    vortices_mat_CELL{i_w} = vortices_mat;

time1_UPDPOS = toc(time1_UPDPOS)

end

%% Compute AIM ( once )
time_AIM = tic;
if i_w == 1
    [AIM,CHORDWISE_AIM,U_ind_mat,V_ind_mat,W_ind_mat] = func_AIM(vortices_mat_0);
else
    [AIM,CHORDWISE_AIM,U_ind_mat,V_ind_mat,W_ind_mat] = func_AIM(vortices_mat);
end
time_AIM = toc(time_AIM)

%% Update Wing Panels from the output of MBDyn


U_REL = zeros(n_rows,n_cols);
V_REL = zeros(n_rows,n_cols);
W_REL = zeros(n_rows,n_cols);


if flag_MBD_IO == 1
if i_w > 2
time_MBD_IO = tic;

   [U_REL,V_REL,W_REL] = func_upd_INFO();
   delete('OUTPUT.dat');
   delete('INI_OUTPUT.dat');
   display('Delete finished');
   vortices_mat_CELL{i_w} = vortices_mat;

time_MBD_IO = toc(time_MBD_IO)

end
else
end



%% Shed wake panels


if i_w ~= 1

 time2_SHEDWAKE = tic;

     if flag_WAKEDIST == 1
     func_shedwake_CW(i_w,vortices_mat_CELL,cw);
     else
     func_shedwake(i_w,vortices_mat_CELL);
     end

     % Assign Gamma from previous time step to wake panels
     for j = 1:n_cols
         gamma_w_mat(i_w-1,j) = gamma_mat_CELL{i_w-1}(n_rows,j);
         NW_wake_vortex_rings{i_w-1,j}.GAMMA = gamma_w_mat(i_w-1,j);
     end

 time2_SHEDWAKE = toc(time2_SHEDWAKE)

else
end

%% Shed LE wake panels

if flag_LEVS == 1
if i_w ~= 1


 time2_1_SHEDWAKE = tic;

     func_LE_shedwake(i_w,vortices_mat_CELL);

     % Assign Gamma from previous time step to wake panels
     for j = 1:n_cols
         LE_gamma_w_mat(i_w-1,j) = -gamma_mat_CELL{i_w-1}(1,j);
         LE_NW_wake_vortex_rings{i_w-1,j}.GAMMA = LE_gamma_w_mat(i_w-1,j);
     end

 time2_1_SHEDWAKE = toc(time2_1_SHEDWAKE)

else
end
end

%% Solve linear system for each time step with changing boundary conditions

time3_SOLVER = tic;
       [gamma_mat_CELL{i_w},W_panels_TV,u_mot_CP,v_mot_CP,w_mot_CP,U_wake,V_wake,W_wake,U_panels_TOT,V_panels_TOT,W_panels_TOT] = func_proc_solver (t,i_w,vortices_mat_CELL{i_w},vortices_mat_0,Quat_CELL{i_w},...
           AIM,CHORDWISE_AIM,U_ind_mat,V_ind_mat,W_ind_mat,U_REL,V_REL,W_REL,flag_LEVS,flag_VORTONS,flag_BENDMOT,sig);
time3_SOLVER = toc(time3_SOLVER)

%% Compute forces and store them into cells and arrays

if i_w ~= 1

time4_FORCES = tic;
  [C_L_v(i_w),C_Di_v(i_w),C_Di_v_st(i_w),C_Di_v_unst(i_w),C_L_mat{i_w},C_Di_mat{i_w},F_mat{i_w},L_mat{i_w},spanload{i_w},adim_spanload{i_w}] = func_forces_KATZPLOTKIN(dt,vortices_mat,gamma_mat_CELL{i_w},gamma_mat_CELL{i_w-1},Q_inf,...
      U_wake,V_wake,W_wake,W_panels_TV,U_panels_TOT,V_panels_TOT,W_panels_TOT,u_mot_CP,v_mot_CP,w_mot_CP);
time4_FORCES = toc(time4_FORCES)

else
time4_FORCES = tic;
  [C_L_v(i_w),C_Di_v(i_w),C_Di_v_st(i_w),C_Di_v_unst(i_w),C_L_mat{i_w},C_Di_mat{i_w},F_mat{i_w},L_mat{i_w},spanload{i_w},adim_spanload{i_w}] = func_forces_KATZPLOTKIN(dt,vortices_mat,gamma_mat_CELL{i_w},zeros(n_rows,n_cols),Q_inf,...
      U_wake,V_wake,W_wake,W_panels_TV,U_panels_TOT,V_panels_TOT,W_panels_TOT,u_mot_CP,v_mot_CP,w_mot_CP);
time4_FORCES = toc(time4_FORCES)

end


%% Create OUTPUT.dat for MBDyn analysis

if flag_MBD_IO == 1

time_INPUT = tic;

 func_INPUT(F_mat{i_w},vortices_mat_0);

time_INPUT = toc(time_INPUT)
else
end

%% Wake to Vortons conversion TE and LE

if flag_VORTONS == 1

  if i_w > 2

   time_RING2VORT = tic;

   func_TE_rings2vortons(i_w,NW_wake_vortex_rings);

   time_RING2VORT = toc(time_RING2VORT);

  end

 if flag_LEVS == 1

  if i_w > 2

   time_RING2VORT2 = tic;

   func_LE_rings2vortons(i_w,LE_NW_wake_vortex_rings);

   time_RING2VORT2 = toc(time_RING2VORT2);

  end
 else
 end

else
end




%% Calculate wake evolution


if flag_VORTONS == 1;


if i_w > 1

if  flag_RU == 1
 if  flag_LEVS == 1

  if i_w > 3

  time5_ROLLUP = tic;
    func_wakerollup_VORTONS_TELE_PARALL(dt,i_w,vortices_mat_CELL{i_w},gamma_mat_CELL{i_w},NW_wake_vortex_rings,LE_NW_wake_vortex_rings,FW_wake_vortons,LE_FW_wake_vortons,sig);
  time5_ROLLUP = toc(time5_ROLLUP)

  else

  time5_ROLLUP = tic;
    func_wakerollup_1_TELE(dt,i_w,vortices_mat_CELL{i_w},gamma_mat_CELL{i_w},NW_wake_vortex_rings,LE_NW_wake_vortex_rings);
  time5_ROLLUP = toc(time5_ROLLUP)

  end

 else
  if i_w > 3

  time5_ROLLUP = tic;
    func_wakerollup_VORTONS_PARALL(dt,i_w,vortices_mat_CELL{i_w},gamma_mat_CELL{i_w},NW_wake_vortex_rings,FW_wake_vortons,sig);
  time5_ROLLUP = toc(time5_ROLLUP)

  else

  time5_ROLLUP = tic;
    func_wakerollup_1(dt,i_w,vortices_mat_CELL{i_w},gamma_mat_CELL{i_w},NW_wake_vortex_rings);
  time5_ROLLUP = toc(time5_ROLLUP)


  end

 end


 else
 end
else

end

else

if i_w ~= 1

if  flag_RU == 1
 if  flag_LEVS == 1

  time5_ROLLUP = tic;
    func_wakerollup_MODIFIED(dt,i_w,vortices_mat_CELL{i_w},gamma_mat_CELL{i_w},NW_wake_vortex_rings,LE_NW_wake_vortex_rings);
  time5_ROLLUP = toc(time5_ROLLUP)
 else

  time5_ROLLUP = tic;
    func_wakerollup(dt,i_w,vortices_mat_CELL{i_w},gamma_mat_CELL{i_w},NW_wake_vortex_rings);
  time5_ROLLUP = toc(time5_ROLLUP)


 end


 else
 end
else

end
end


%% Time recording

timerec_SOLVER(i_w) = time3_SOLVER;
step_time(i_w)       = time1_UPDPOS + time3_SOLVER + time4_FORCES + time_AIM ;

if flag_MBD_IO == 1
est_time(i_w) = (time1_UPDPOS + time3_SOLVER + time_AIM + time_MBD_IO +waiting_time) * (N_t-i_w);
else
est_time(i_w) = (time1_UPDPOS + time3_SOLVER + time_AIM) * (N_t-i_w);
end


if flag_RU == 1
    if i_w > 1
est_time(i_w) = (time1_UPDPOS + time3_SOLVER + time_AIM + time5_ROLLUP) * (N_t-i_w);
    else
est_time(i_w) = (time1_UPDPOS + time3_SOLVER + time_AIM) * (N_t-i_w);
    end
else
est_time(i_w) = (time1_UPDPOS + time3_SOLVER + time_AIM) * (N_t-i_w);
end


%% Display current time-step calculations

fprintf('\n');
disp(['C_L is: ' num2str(C_L_v(i_w))]);
fprintf('\n');
disp(['U is: ' num2str(U(t))]);
disp(['V is: ' num2str(V(t))]);
disp(['W is: ' num2str(W(t))]);
disp(['p is: ' num2str(p(t))]);
disp(['q is: ' num2str(q(t))]);
disp(['r is: ' num2str(r(t))]);

if flag_MBD_IO == 1
fprintf('\n');
disp(['Time(in seconds) of current step was about : ' num2str(step_time(i_w)+time_MBD_IO+waiting_time)]);
disp(['Estimated minutes to end of computation is : ' num2str(est_time(i_w)/60)]);
fprintf('\n\n');
else
fprintf('\n');
disp(['Time(in seconds) of current step was about : ' num2str(step_time(i_w))]);
disp(['Estimated minutes to end of computation is : ' num2str(est_time(i_w)/60)]);
fprintf('\n\n');
end

%% Create Frame

if flag_FRAME == 1

time_FRAME = tic;
FRAME(i_w) = func_frame(vortices_mat,C_L_mat{i_w},i_w,NW_wake_vortex_rings,LE_NW_wake_vortex_rings,FW_wake_vortons,LE_FW_wake_vortons);
time_FRAME = toc(time_FRAME)

else
end

%% Saving current workspace to a .mat file

save('SIM_TEMP.mat');

end

%% Saving final workspace
save('SIM_FIN.mat');
tot_time = toc(tot_time);

disp(['Total elapsed time is : ' num2str(tot_time)]);
