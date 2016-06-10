function [AIM,CHORDWISE_AIM,U_ind_mat,V_ind_mat,W_ind_mat] = func_AIM(vortices_mat)

%% TO FORM AERODYNAMIC INFLUENCE MATRIX FOR THE VORTEX LATTICE

global n_rows n_cols 

%% Adding paths

addpath([pwd,'/functions']);

%% To import vortex lattice

N = n_rows*n_cols;
vortex_rings = cell(1,N);

for i = 1:n_rows
    for j = 1:n_cols
vortex_rings(1,j+n_cols*(i-1)) = vortices_mat(i,j);
    end
end


%% To form the AIM

AIM = zeros(N,N);
CHORDWISE_AIM = zeros(N,N);
U_ind_mat= zeros(N,N);
V_ind_mat= zeros(N,N);
W_ind_mat= zeros(N,N);


GAMMA = 1.0;
ONOFF = [0,1,2];

% loop on Collocation Points and Vortex Rings


parfor i = 1 : N
    
    X_CP = vortex_rings{1,i}.X_C;
    Y_CP = vortex_rings{1,i}.Y_C;
    Z_CP = vortex_rings{1,i}.Z_C;
    N_x  = vortex_rings{1,i}.N_x;
    N_y  = vortex_rings{1,i}.N_y;
    N_z  = vortex_rings{1,i}.N_z;
    
%     % Middle point of bound vortices sides (DE-ACTIVATED)
%     
%     X_BV = (vortex_rings{1,i}.X_1+vortex_rings{1,i}.X_2)/2;
%     Y_BV = (vortex_rings{1,i}.Y_1+vortex_rings{1,i}.Y_2)/2;
%     Z_BV = (vortex_rings{1,i}.Z_1+vortex_rings{1,i}.Z_2)/2;
        
    X_BV = X_CP;
    Y_BV = Y_CP;
    Z_BV = Z_CP;
    
    
    for j = 1 : N
   
if j >= (n_rows-1)*n_cols+1
    
      ind = func_voring(X_CP,Y_CP,Z_CP,vortex_rings{1,j},GAMMA,ONOFF(1));
      ind_image = func_voring(X_CP,-Y_CP,Z_CP,vortex_rings{1,j},GAMMA,ONOFF(1));
      c_ind = func_voring(X_CP,Y_CP,Z_CP,vortex_rings{1,j},GAMMA,ONOFF(3));
      c_ind_image = func_voring(X_CP,-Y_CP,Z_CP,vortex_rings{1,j},GAMMA,ONOFF(3));
      
   AIM(i,j)           = (ind(1)+ind_image(1))*N_x + (ind(2)-ind_image(2))*N_y + (ind(3)+ind_image(3))*N_z;

       ind_BV = func_voring(X_BV,Y_BV,Z_BV,vortex_rings{1,j},GAMMA,ONOFF(1));
       ind_image_BV = func_voring(X_BV,-Y_BV,Z_BV,vortex_rings{1,j},GAMMA,ONOFF(1));


      
   U_ind_mat(i,j)              = (ind_BV(1)+ind_image_BV(1));
   V_ind_mat(i,j)              = (ind_BV(2)-ind_image_BV(2));
   W_ind_mat(i,j)              = (ind_BV(3)+ind_image_BV(3));
   CHORDWISE_AIM(i,j) = (c_ind(1)+c_ind_image(1))*N_x + ...
       (c_ind(2)-c_ind_image(2))*N_y +(c_ind(3)+c_ind_image(3))*N_z;

    
else
    
      ind = func_voring(X_CP,Y_CP,Z_CP,vortex_rings{1,j},GAMMA,ONOFF(1));
      ind_image = func_voring(X_CP,-Y_CP,Z_CP,vortex_rings{1,j},GAMMA,ONOFF(1));
      c_ind = func_voring(X_CP,Y_CP,Z_CP,vortex_rings{1,j},GAMMA,ONOFF(2));
      c_ind_image = func_voring(X_CP,-Y_CP,Z_CP,vortex_rings{1,j},GAMMA,ONOFF(2));
      
   AIM(i,j)           = (ind(1)+ind_image(1))*N_x + (ind(2)-ind_image(2))*N_y + (ind(3)+ind_image(3))*N_z;

       ind_BV = func_voring(X_BV,Y_BV,Z_BV,vortex_rings{1,j},GAMMA,ONOFF(1));
       ind_image_BV = func_voring(X_BV,-Y_BV,Z_BV,vortex_rings{1,j},GAMMA,ONOFF(1));


      
   U_ind_mat(i,j)              = (ind_BV(1)+ind_image_BV(1));
   V_ind_mat(i,j)              = (ind_BV(2)-ind_image_BV(2));
   W_ind_mat(i,j)              = (ind_BV(3)+ind_image_BV(3));
   CHORDWISE_AIM(i,j) = (c_ind(1)+c_ind_image(1))*N_x + ...
       (c_ind(2)-c_ind_image(2))*N_y +(c_ind(3)+c_ind_image(3))*N_z;
 
    
end
      
   
    end
end

end