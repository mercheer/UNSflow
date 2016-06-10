function [X,Y,Z,U_SL,V_SL,W_SL,sx,sy,sz] = func_domain_induction(vortices_mat_CELL,gamma_mat,Q_inf,N_t)
%func_domain_induction compute domain points and velocities
%   Detailed explanation goes here

global n_rows n_cols wake_vortex_rings chord

%% Adding path

addpath([pwd,'/functions']);


% PAY ATTENTION THIS IS VALID ONLY FOR FORWARD MOTION IN -X DIRECTION


ang = atan2(abs(vortices_mat_CELL{1,1}.Z_3-vortices_mat_CELL{1,1}.Z_1),abs(vortices_mat_CELL{1,1}.X_3-vortices_mat_CELL{1,1}.X_1));



LENGTH_Z = 3;



wake_vortex_rings_buff = wake_vortex_rings;

U_SL = zeros(n_cols,n_rows,LENGTH_Z);
V_SL = zeros(n_cols,n_rows,LENGTH_Z);
W_SL = zeros(n_cols,n_rows,LENGTH_Z);




for j = 1:n_cols
    
    X_START(j) = vortices_mat_CELL{1,j}.X_C;  
    Y_P(j) = vortices_mat_CELL{1,j}.Y_C;

    for i = 1:n_rows
        
        
        for k = 1:LENGTH_Z
        
            
        X_buff = vortices_mat_CELL{i,j}.X_C;
        Y_buff = vortices_mat_CELL{i,j}.Y_C;
        
        switch k 
            case 1    
                Z_buff = vortices_mat_CELL{i,j}.Z_C + chord(1)*0.05;
            case 2
                Z_buff = vortices_mat_CELL{i,j}.Z_C;
            case 3
                Z_buff = vortices_mat_CELL{i,j}.Z_C - chord(1)*0.05;

        end

 % Loop for wake induction on POINTS( Right semi-span and left semi-span)

    u_w_tot = 0;
    v_w_tot = 0;
    w_w_tot = 0;
    u_wing_tot = 0;
    v_wing_tot = 0;
    w_wing_tot = 0;
            
   for l = 1:N_t-1 
       for m = 1:n_cols 
        [u_w_add_R,v_w_add_R,w_w_add_R] = func_voring_comp(X_buff,Y_buff,Z_buff,wake_vortex_rings_buff{l,m},...
                    wake_vortex_rings_buff{l,m}.GAMMA,0);
        [u_w_add_L,v_w_add_L,w_w_add_L] = func_voring_comp(X_buff,-Y_buff,Z_buff,wake_vortex_rings_buff{l,m},...
                    wake_vortex_rings_buff{l,m}.GAMMA,0);
                
        u_w_tot = u_w_tot  + u_w_add_R +u_w_add_L;
        v_w_tot = v_w_tot + v_w_add_R -v_w_add_L;
        w_w_tot = w_w_tot + w_w_add_R +w_w_add_L;
        end
    end
    
   for a = 1:n_rows
       for b = 1:n_cols
        [u_wing_add_R,v_wing_add_R,w_wing_add_R] = func_voring_comp(X_buff,Y_buff,Z_buff,vortices_mat_CELL{a,b},...
                    gamma_mat(a,b),0);
        [u_wing_add_L,v_wing_add_L,w_wing_add_L] = func_voring_comp(X_buff,-Y_buff,Z_buff,vortices_mat_CELL{a,b},...
                    gamma_mat(a,b),0);
                
        u_wing_tot = u_wing_tot  + u_wing_add_R +u_wing_add_L;
        v_wing_tot = v_wing_tot + v_wing_add_R -v_wing_add_L;
        w_wing_tot = w_wing_tot + w_wing_add_R +w_wing_add_L;    
           
           
       end
   end
   
   X(j,i,k) = X_buff;
   Y(j,i,k) = Y_buff;
   Z(j,i,k) = Z_buff;
   
   Z_P(k) = Z_buff;
   
   U_SL(j,i,k) = u_w_tot +u_wing_tot + Q_inf(1);
   V_SL(j,i,k) = v_w_tot +v_wing_tot + Q_inf(2);
   W_SL(j,i,k) = w_w_tot +w_wing_tot + Q_inf(3);
        end
    end
end



sx = [X_START;X_START;X_START];
sy = [Y_P;Y_P;Y_P];
sz = [ones(1,n_cols).*Z_P(1);ones(1,n_cols).*Z_P(2);ones(1,n_cols).*Z_P(3)];

end

