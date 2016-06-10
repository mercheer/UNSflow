function [gamma_mat,W_panels_TV,u_mot_CP_mat,v_mot_CP_mat,w_mot_CP_mat,U_wake,V_wake,W_wake,U_panels_TOT,V_panels_TOT,W_panels_TOT] = func_proc_solver (t,now_index,vortices_mat_CELL,vortices_mat_0,Quat,AIM,CHORDWISE_AIM,U_ind_mat,V_ind_mat,W_ind_mat,U_REL,V_REL,W_REL,flag_LEVS,flag_VORTONS,flag_BENDMOT,sig)

% TO SOLVE THE LINEAR SYSTEM WITH ARBITRARY RHS

global NW_wake_vortex_rings LE_NW_wake_vortex_rings FW_wake_vortons LE_FW_wake_vortons... 
       s chord n_rows n_cols 

%% Adding paths

addpath([pwd,'/geometry']);
addpath([pwd,'/functions']);
addpath([pwd,'/data']);


%% To load AIM,panels,airfoil

load([pwd,'/data/data_airfoil.mat']);
load([pwd,'/data/data_refs.mat']);
load([pwd,'/data/data_kinematics.mat']);



%% To create RHS in BRF

N =  n_rows * n_cols;
RHS = zeros(N,1);    % RHS is a column vector

wake_ind_IRF    = cell(1,n_rows*n_cols);
LE_wake_ind_IRF = cell(1,n_rows*n_cols);
wake_ind_BRF    = cell(1,n_rows*n_cols);    
u_mot_CP        = zeros(1,n_rows*n_cols);
v_mot_CP        = zeros(1,n_rows*n_cols);
w_mot_CP        = zeros(1,n_rows*n_cols);
u_mot_CP_mat        = zeros(n_rows,n_cols);
v_mot_CP_mat        = zeros(n_rows,n_cols);
w_mot_CP_mat        = zeros(n_rows,n_cols);

rows = n_rows;
cols = n_cols;
s_buff = s;
chord_buff = chord;

wake_vortex_rings_buff    = NW_wake_vortex_rings;
LE_wake_vortex_rings_buff = LE_NW_wake_vortex_rings;
FW_wake_vortons_buff      = FW_wake_vortons;
LE_FW_wake_vortons_buff   = LE_FW_wake_vortons;

xx = xx;
yyc = yyc;
b_ref = b_ref;
twist = twist;

if flag_BENDMOT == 1
W_BEND = W_BEND;
else
W_BEND = 0;
end

parfor k = 1:N
        
    i = ceil(k/cols);
    j = k - cols*(i-1);
    
    N_x      = vortices_mat_CELL{i,j}.N_x;
    N_y      = vortices_mat_CELL{i,j}.N_y;
    N_z      = vortices_mat_CELL{i,j}.N_z;
    X_CP_IRF = vortices_mat_CELL{i,j}.X_C;
    Y_CP_IRF = vortices_mat_CELL{i,j}.Y_C;
    Z_CP_IRF = vortices_mat_CELL{i,j}.Z_C;
    
    X_CP_BRF = vortices_mat_0{i,j}.X_C;
    Y_CP_BRF = vortices_mat_0{i,j}.Y_C;
    Z_CP_BRF = vortices_mat_0{i,j}.Z_C;
    
    col = j;
    
    % To add effect of camber-line,twist and dihedral combination
    
    dihedral_loc = atan2((vortices_mat_CELL{i,j}.Z_2-vortices_mat_CELL{i,j}.Z_1),...
                          (vortices_mat_CELL{i,j}.Y_2-vortices_mat_CELL{i,j}.Y_1));
   
    c_CP   = interp1(s_buff,chord_buff,col,'linear','extrap');
    X_LE   = (vortices_mat_CELL{1,j}.X_2+vortices_mat_CELL{1,j}.X_1)/2-c_CP/rows*0.25;
    x_CP   = (X_CP_BRF-X_LE)/c_CP;
    
    yyc_CP = -interp1(xx,yyc,x_CP,'linear','extrap');
     
    twist_CP = vortices_mat_CELL{i,j}.TWIST + interp1([0 b_ref/2],[0 twist(t)],Y_CP_BRF,'linear','extrap');
    rot_angle = yyc_CP + twist_CP;
    
    [N_x_rot,N_y_rot,N_z_rot] = func_rodrot(rot_angle,0,cos(dihedral_loc),sin(dihedral_loc)...
        ,N_x,N_y,N_z);
    
    % To add effect of kinematics of the wing in IRF (Kinematics is given
    % in the initial RF that is the BRF )
    
    [u_mot_CP(k),v_mot_CP(k),w_mot_CP(k)] = func_kinematics_vel(t,X_CP_BRF,Y_CP_BRF,Z_CP_BRF,[1,0,0,0]);
    
    % To add effect of structural deformations velocities
    
    u_mot_CP(k) = u_mot_CP(k) + U_REL(k);
    v_mot_CP(k) = v_mot_CP(k) + V_REL(k);
    
    if flag_BENDMOT == 1
    
    w_mot_CP(k) = w_mot_CP(k) + W_REL(k) + W_BEND(t,Y_CP_BRF);
    
    else
        
    w_mot_CP(k) = w_mot_CP(k) + W_REL(k);
    
    end

    
if now_index > 1    
    
%% Loop for NEAR WAKE induction on C_p( Right semi-span and left semi-span)
 if flag_VORTONS == 1
    
    l = now_index-1;
    u_w_tot = 0;
    v_w_tot = 0;
    w_w_tot = 0;
        
        for m = 1:cols 
        [u_w_add_R,v_w_add_R,w_w_add_R] = func_voring_comp(X_CP_IRF,Y_CP_IRF,Z_CP_IRF,wake_vortex_rings_buff{l,m},...
                    wake_vortex_rings_buff{l,m}.GAMMA,0);
        [u_w_add_L,v_w_add_L,w_w_add_L] = func_voring_comp(X_CP_IRF,-Y_CP_IRF,Z_CP_IRF,wake_vortex_rings_buff{l,m},...
                    wake_vortex_rings_buff{l,m}.GAMMA,0);
        u_w_tot = u_w_tot  + u_w_add_R +u_w_add_L;
        v_w_tot = v_w_tot + v_w_add_R -v_w_add_L;
        w_w_tot = w_w_tot + w_w_add_R +w_w_add_L;
        end
        
        if now_index > 2
        l = now_index-2;    
        for m = 1:cols 
        [u_w_add_R,v_w_add_R,w_w_add_R] = func_voring_comp(X_CP_IRF,Y_CP_IRF,Z_CP_IRF,wake_vortex_rings_buff{l,m},...
                    wake_vortex_rings_buff{l,m}.GAMMA,0);
        [u_w_add_L,v_w_add_L,w_w_add_L] = func_voring_comp(X_CP_IRF,-Y_CP_IRF,Z_CP_IRF,wake_vortex_rings_buff{l,m},...
                    wake_vortex_rings_buff{l,m}.GAMMA,0);
        u_w_tot = u_w_tot  + u_w_add_R +u_w_add_L;
        v_w_tot = v_w_tot + v_w_add_R -v_w_add_L;
        w_w_tot = w_w_tot + w_w_add_R +w_w_add_L;
        end
        else
        end
        
    
%% Loop for NEAR WAKE induction on C_p( Right semi-span and left semi-span)
    
    l = now_index-1;
    u_w_LE_tot = 0;
    v_w_LE_tot = 0;
    w_w_LE_tot = 0;
    
   if flag_LEVS == 1 

        for m = 1:cols
            
        [u_w_LE_add_R,v_w_LE_add_R,w_w_LE_add_R] = func_voring_comp(X_CP_IRF,Y_CP_IRF,Z_CP_IRF,LE_wake_vortex_rings_buff{l,m},...
                    LE_wake_vortex_rings_buff{l,m}.GAMMA,0);
        [u_w_LE_add_L,v_w_LE_add_L,w_w_LE_add_L] = func_voring_comp(X_CP_IRF,-Y_CP_IRF,Z_CP_IRF,LE_wake_vortex_rings_buff{l,m},...
                    LE_wake_vortex_rings_buff{l,m}.GAMMA,0);
        u_w_LE_tot = u_w_LE_tot  + u_w_LE_add_R +u_w_LE_add_L;
        v_w_LE_tot = v_w_LE_tot + v_w_LE_add_R -v_w_LE_add_L;
        w_w_LE_tot = w_w_LE_tot + w_w_LE_add_R +w_w_LE_add_L;
        
        end
        
        if now_index > 2
        l = now_index-2;    
        for m = 1:cols
            
        [u_w_LE_add_R,v_w_LE_add_R,w_w_LE_add_R] = func_voring_comp(X_CP_IRF,Y_CP_IRF,Z_CP_IRF,LE_wake_vortex_rings_buff{l,m},...
                    LE_wake_vortex_rings_buff{l,m}.GAMMA,0);
        [u_w_LE_add_L,v_w_LE_add_L,w_w_LE_add_L] = func_voring_comp(X_CP_IRF,-Y_CP_IRF,Z_CP_IRF,LE_wake_vortex_rings_buff{l,m},...
                    LE_wake_vortex_rings_buff{l,m}.GAMMA,0);
        u_w_LE_tot = u_w_LE_tot  + u_w_LE_add_R +u_w_LE_add_L;
        v_w_LE_tot = v_w_LE_tot + v_w_LE_add_R -v_w_LE_add_L;
        w_w_LE_tot = w_w_LE_tot + w_w_LE_add_R +w_w_LE_add_L;
        end
        else
        end
        
        
        
   end
   
   
   
   
   
%% Loop for FAR WAKE induction on C_P ( Right semi-span and left semi-span )
   
   
    l = 1;
    FW_u_w_tot = 0;
    FW_v_w_tot = 0;
    FW_w_w_tot = 0;
    
     while l < now_index-2
        for m = 1:cols 
        [u_w_add_R,v_w_add_R,w_w_add_R] = func_vortonind(X_CP_IRF,Y_CP_IRF,Z_CP_IRF,FW_wake_vortons_buff{l,m},sig);
        [u_w_add_L,v_w_add_L,w_w_add_L] = func_vortonind(X_CP_IRF,-Y_CP_IRF,Z_CP_IRF,FW_wake_vortons_buff{l,m},sig);
        FW_u_w_tot = FW_u_w_tot  + u_w_add_R +u_w_add_L;
        FW_v_w_tot = FW_v_w_tot + v_w_add_R -v_w_add_L;
        FW_w_w_tot = FW_w_w_tot + w_w_add_R +w_w_add_L;
        end
     l = l +1 ;
     end
   
%% Loop for FAR WAKE induction on C_p( Right semi-span and left semi-span)
  
  
    l = 1;
    FW_u_w_LE_tot = 0;
    FW_v_w_LE_tot = 0;
    FW_w_w_LE_tot = 0;
    
   if flag_LEVS == 1 
        
       while l < now_index-2
         for m = 1:cols 
         [u_w_LE_add_R,v_w_LE_add_R,w_w_LE_add_R] = func_vortonind(X_CP_IRF,Y_CP_IRF,Z_CP_IRF,LE_FW_wake_vortons_buff{l,m},sig);
         [u_w_LE_add_L,v_w_LE_add_L,w_w_LE_add_L] = func_vortonind(X_CP_IRF,-Y_CP_IRF,Z_CP_IRF,LE_FW_wake_vortons_buff{l,m},sig);
         FW_u_w_LE_tot = FW_u_w_LE_tot  + u_w_LE_add_R +u_w_LE_add_L;
         FW_v_w_LE_tot = FW_v_w_LE_tot + v_w_LE_add_R -v_w_LE_add_L;
         FW_w_w_LE_tot = FW_w_w_LE_tot + w_w_LE_add_R +w_w_LE_add_L;
         end
        l = l +1 ;
        end
        
   end
   
%% Loop for NEAR-WAKE - FAR-WAKE separation LINE
         
   u_w_NWFW_tot = 0;
   v_w_NWFW_tot = 0;
   w_w_NWFW_tot = 0;
   
   
   if now_index > 3
         
    for m = 1:cols
        X_1 = wake_vortex_rings_buff{now_index-2,m}.X_3;
        X_2 = wake_vortex_rings_buff{now_index-2,m}.X_4;
        Y_1 = wake_vortex_rings_buff{now_index-2,m}.Y_3;
        Y_2 = wake_vortex_rings_buff{now_index-2,m}.Y_4;        
        Z_1 = wake_vortex_rings_buff{now_index-2,m}.Z_3;
        Z_2 = wake_vortex_rings_buff{now_index-2,m}.Z_4;
        GAMMA = wake_vortex_rings_buff{now_index-3,m}.GAMMA;
        
         [u_w_NWFW_add_R,v_w_NWFW_add_R,w_w_NWFW_add_R] = func_vortexl(X_CP_IRF,Y_CP_IRF,Z_CP_IRF,X_1,Y_1,Z_1,X_2,Y_2,Z_2,GAMMA);
         [u_w_NWFW_add_L,v_w_NWFW_add_L,w_w_NWFW_add_L] = func_vortexl(X_CP_IRF,-Y_CP_IRF,Z_CP_IRF,X_1,Y_1,Z_1,X_2,Y_2,Z_2,GAMMA);
         u_w_NWFW_tot = u_w_NWFW_tot  + u_w_NWFW_add_R +u_w_NWFW_add_L;
         v_w_NWFW_tot = v_w_NWFW_tot + v_w_NWFW_add_R -v_w_NWFW_add_L;
         w_w_NWFW_tot = w_w_NWFW_tot + w_w_NWFW_add_R +w_w_NWFW_add_L;
    end
   else
   end
   
   LE_u_w_NWFW_tot = 0;
   LE_v_w_NWFW_tot = 0;
   LE_w_w_NWFW_tot = 0;

  if flag_LEVS == 1 
   if now_index > 3
         
     for m = 1:cols
        X_1 = LE_wake_vortex_rings_buff{now_index-2,m}.X_3;
        X_2 = LE_wake_vortex_rings_buff{now_index-2,m}.X_4;
        Y_1 = LE_wake_vortex_rings_buff{now_index-2,m}.Y_3;
        Y_2 = LE_wake_vortex_rings_buff{now_index-2,m}.Y_4;        
        Z_1 = LE_wake_vortex_rings_buff{now_index-2,m}.Z_3;
        Z_2 = LE_wake_vortex_rings_buff{now_index-2,m}.Z_4;
        GAMMA = LE_wake_vortex_rings_buff{now_index-3,m}.GAMMA;
        
         [LE_u_w_NWFW_add_R,LE_v_w_NWFW_add_R,LE_w_w_NWFW_add_R] = func_vortexl(X_CP_IRF,Y_CP_IRF,Z_CP_IRF,X_1,Y_1,Z_1,X_2,Y_2,Z_2,GAMMA);
         [LE_u_w_NWFW_add_L,LE_v_w_NWFW_add_L,LE_w_w_NWFW_add_L] = func_vortexl(X_CP_IRF,-Y_CP_IRF,Z_CP_IRF,X_1,Y_1,Z_1,X_2,Y_2,Z_2,GAMMA);
         LE_u_w_NWFW_tot = LE_u_w_NWFW_tot  + LE_u_w_NWFW_add_R +LE_u_w_NWFW_add_L;
         LE_v_w_NWFW_tot = LE_v_w_NWFW_tot + LE_v_w_NWFW_add_R -LE_v_w_NWFW_add_L;
         LE_w_w_NWFW_tot = LE_w_w_NWFW_tot + LE_w_w_NWFW_add_R +LE_w_w_NWFW_add_L;
     end
   else
   end
  else
  end
     
 else
     
    FW_u_w_tot = 0;
    FW_v_w_tot = 0;
    FW_w_w_tot = 0;
    FW_u_w_LE_tot = 0;
    FW_v_w_LE_tot = 0;
    FW_w_w_LE_tot = 0;
    u_w_NWFW_tot = 0;
    v_w_NWFW_tot = 0;
    w_w_NWFW_tot = 0;
    LE_u_w_NWFW_tot = 0;
    LE_v_w_NWFW_tot = 0;
    LE_w_w_NWFW_tot = 0; 
    
%% Loop for wake induction on C_p( Right semi-span and left semi-span)
    l = 1;
    u_w_tot = 0;
    v_w_tot = 0;
    w_w_tot = 0;

    
    while l < now_index
        for m = 1:cols 
        [u_w_add_R,v_w_add_R,w_w_add_R] = func_voring_comp(X_CP_IRF,Y_CP_IRF,Z_CP_IRF,wake_vortex_rings_buff{l,m},...
                    wake_vortex_rings_buff{l,m}.GAMMA,0);
        [u_w_add_L,v_w_add_L,w_w_add_L] = func_voring_comp(X_CP_IRF,-Y_CP_IRF,Z_CP_IRF,wake_vortex_rings_buff{l,m},...
                    wake_vortex_rings_buff{l,m}.GAMMA,0);
        u_w_tot = u_w_tot  + u_w_add_R +u_w_add_L;
        v_w_tot = v_w_tot + v_w_add_R -v_w_add_L;
        w_w_tot = w_w_tot + w_w_add_R +w_w_add_L;
        end
        l = l + 1;
    end
    
%% Loop for LEADING EDGE wake induction on C_p( Right semi-span and left semi-span)

    l = 1;
    u_w_LE_tot = 0;
    v_w_LE_tot = 0;
    w_w_LE_tot = 0;
    
   if flag_LEVS == 1 
   while l < now_index
        for m = 1:cols 
        [u_w_LE_add_R,v_w_LE_add_R,w_w_LE_add_R] = func_voring_comp(X_CP_IRF,Y_CP_IRF,Z_CP_IRF,LE_wake_vortex_rings_buff{l,m},...
                    LE_wake_vortex_rings_buff{l,m}.GAMMA,0);
        [u_w_LE_add_L,v_w_LE_add_L,w_w_LE_add_L] = func_voring_comp(X_CP_IRF,-Y_CP_IRF,Z_CP_IRF,LE_wake_vortex_rings_buff{l,m},...
                    LE_wake_vortex_rings_buff{l,m}.GAMMA,0);
        u_w_LE_tot = u_w_LE_tot  + u_w_LE_add_R +u_w_LE_add_L;
        v_w_LE_tot = v_w_LE_tot + v_w_LE_add_R -v_w_LE_add_L;
        w_w_LE_tot = w_w_LE_tot + w_w_LE_add_R +w_w_LE_add_L;
        end
        l = l + 1;
   end
   end
   
   
   
   
 end
 
%% Wake inductions in IRF
    
     wake_ind_IRF{k} = [u_w_tot+FW_u_w_tot+u_w_NWFW_tot;v_w_tot+FW_v_w_tot+v_w_NWFW_tot;w_w_tot+FW_w_w_tot+w_w_NWFW_tot];
     
     LE_wake_ind_IRF{k} = [u_w_LE_tot+FW_u_w_LE_tot+LE_u_w_NWFW_tot;v_w_LE_tot+FW_v_w_LE_tot+LE_v_w_NWFW_tot;w_w_LE_tot+FW_w_w_LE_tot+LE_w_w_NWFW_tot];

%% Wake inductions in BRF 
    
     wake_ind_BRF{k} = quat2dcm(Quat)*[wake_ind_IRF{k}(1);wake_ind_IRF{k}(2);wake_ind_IRF{k}(3)]; 
 
else
   
    wake_ind_IRF{k}    = [0;0;0];
    LE_wake_ind_IRF{k} = [0;0;0];
    wake_ind_BRF{k}    = [0;0;0];
    
end

    
%% RHS computation as Q_motion + Q_wake in IRF
    
    RHS(k) = -((+wake_ind_IRF{k}(1)+LE_wake_ind_IRF{k}(1)+u_mot_CP(k))*N_x_rot +...
        (+wake_ind_IRF{k}(2)+LE_wake_ind_IRF{k}(2)+v_mot_CP(k))*N_y_rot + ...
        (+wake_ind_IRF{k}(3)+LE_wake_ind_IRF{k}(3)+w_mot_CP(k))*N_z_rot);
    
    
end

%% To solve the system

gamma = AIM\RHS;


%% To compute wing self-inductions( and chordwise vortices inductions)

U_ind = U_ind_mat * gamma;
V_ind = V_ind_mat * gamma;
W_ind = W_ind_mat * gamma;
W_ind_TV = CHORDWISE_AIM*gamma;

%% To re-order into a matrix gamma and calculate induction in the CP

gamma_mat = zeros(n_rows,n_cols);
W_panels_TV = zeros(n_rows,n_cols); % Trailing Vortices induction to compute induced drag
U_wake = zeros(n_rows,n_cols);
V_wake = zeros(n_rows,n_cols);
W_wake = zeros(n_rows,n_cols);
U_panels_TOT = zeros(n_rows,n_cols);
V_panels_TOT = zeros(n_rows,n_cols);
W_panels_TOT = zeros(n_rows,n_cols);
V_wing_SELFIND = cell(n_rows,n_cols);



for k = 1:rows*cols
    i = ceil(k/cols);
    j = k - cols*(i-1);
    u_mot_CP_mat(i,j) = u_mot_CP(k); 
    v_mot_CP_mat(i,j) = v_mot_CP(k); 
    w_mot_CP_mat(i,j) = w_mot_CP(k); 
end

k = 1;

for i = 1:rows

    for j = 1:cols        
        
gamma_mat(i,j) = gamma(j+cols*(i-1));

V_wing_SELFIND{i,j} =  [U_ind(j+cols*(i-1));...
                        V_ind(j+cols*(i-1));...
                        W_ind(j+cols*(i-1))];
                    
W_panels_TV(i,j) = W_ind_TV(j+cols*(i-1));
                    
% VELOCITIES ON CP IN IRF (THESE VELOCITIES ARE ONLY DUE TO WAKE AND WING MOTION )    
% Notice that wake induction should be computed on bound vortices instead
% of on control points
   

U_wake(i,j) =  wake_ind_IRF{k}(1) + LE_wake_ind_IRF{k}(1);
V_wake(i,j) =  wake_ind_IRF{k}(2) + LE_wake_ind_IRF{k}(2);
W_wake(i,j) =  wake_ind_IRF{k}(3) + LE_wake_ind_IRF{k}(3);

U_panels_TOT(i,j) = u_mot_CP_mat(i,j) + wake_ind_IRF{k}(1) + LE_wake_ind_IRF{k}(1) + V_wing_SELFIND{i,j}(1);
V_panels_TOT(i,j) = v_mot_CP_mat(i,j) + wake_ind_IRF{k}(2) + LE_wake_ind_IRF{k}(2) + V_wing_SELFIND{i,j}(2);
W_panels_TOT(i,j) = w_mot_CP_mat(i,j) + wake_ind_IRF{k}(3) + LE_wake_ind_IRF{k}(3) + V_wing_SELFIND{i,j}(3);

k = k + 1;
    end
    
end


end