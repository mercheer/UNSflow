function func_wakerollup_1(dt,now_index,vortices_CELL,gamma_mat,wake_vortex_ring_now)

global NW_wake_vortex_rings 

%% ADDPATH

addpath([pwd,'/geometry']);
addpath([pwd,'/functions']);
addpath([pwd,'/data']);

%% LOAD

load([pwd,'/data/data_vortex_lattice.mat']);

%% LOOPS

if now_index == 2
    
  i = now_index-1;
else
  
  i = now_index-2;
end

while i < now_index
  if abs(now_index-i)<30
    if i  == 1
     for j = 1:n_cols
        if j == 1
        X_1 = NW_wake_vortex_rings{i,j}.X_1;
        X_2 = NW_wake_vortex_rings{i,j}.X_2;
        X_3 = NW_wake_vortex_rings{i,j}.X_3;
        X_4 = NW_wake_vortex_rings{i,j}.X_4;
        Y_1 = NW_wake_vortex_rings{i,j}.Y_1;
        Y_2 = NW_wake_vortex_rings{i,j}.Y_2;
        Y_3 = NW_wake_vortex_rings{i,j}.Y_3;
        Y_4 = NW_wake_vortex_rings{i,j}.Y_4;
        Z_1 = NW_wake_vortex_rings{i,j}.Z_1;
        Z_2 = NW_wake_vortex_rings{i,j}.Z_2;
        Z_3 = NW_wake_vortex_rings{i,j}.Z_3;
        Z_4 = NW_wake_vortex_rings{i,j}.Z_4;
         
        % Loop for inducted velocity of bound vortex rings
u1 = 0;
u2 = 0;
u3 = 0;
u4 = 0;
v1 = 0;
v2 = 0;
v3 = 0;
v4 = 0;
w1 = 0;
w2 = 0;
w3 = 0;
w4 = 0;
u1w = 0;
u2w = 0;
u3w = 0;
u4w = 0;
v1w = 0;
v2w = 0;
v3w = 0;
v4w = 0;
w1w = 0;
w2w = 0;
w3w = 0;
w4w = 0;
u1w_LE = 0;
u2w_LE = 0;
u3w_LE = 0;
u4w_LE = 0;
v1w_LE = 0;
v2w_LE = 0;
v3w_LE = 0;
v4w_LE = 0;
w1w_LE = 0;
w2w_LE = 0;
w3w_LE = 0;
w4w_LE = 0;

        for l = 1:n_rows
           for m = 1:n_cols
               
               [u1_ind_add_R,v1_ind_add_R,w1_ind_add_R] = func_voring_comp(X_1,Y_1,Z_1,vortices_CELL{l,m},...
                   gamma_mat(l,m),0);
               [u2_ind_add_R,v2_ind_add_R,w2_ind_add_R] = func_voring_comp(X_2,Y_2,Z_2,vortices_CELL{l,m},...
                   gamma_mat(l,m),0);
               [u3_ind_add_R,v3_ind_add_R,w3_ind_add_R] = func_voring_comp(X_3,Y_3,Z_3,vortices_CELL{l,m},...
                   gamma_mat(l,m),0);
               [u4_ind_add_R,v4_ind_add_R,w4_ind_add_R] = func_voring_comp(X_4,Y_4,Z_4,vortices_CELL{l,m},...
                   gamma_mat(l,m),0);
               [u1_ind_add_L,v1_ind_add_L,w1_ind_add_L] = func_voring_comp(X_1,-Y_1,Z_1,vortices_CELL{l,m},...
                   gamma_mat(l,m),0);
               [u2_ind_add_L,v2_ind_add_L,w2_ind_add_L] = func_voring_comp(X_2,-Y_2,Z_2,vortices_CELL{l,m},...
                   gamma_mat(l,m),0);
               [u3_ind_add_L,v3_ind_add_L,w3_ind_add_L] = func_voring_comp(X_3,-Y_3,Z_3,vortices_CELL{l,m},...
                   gamma_mat(l,m),0);
               [u4_ind_add_L,v4_ind_add_L,w4_ind_add_L] = func_voring_comp(X_4,-Y_4,Z_4,vortices_CELL{l,m},...
                   gamma_mat(l,m),0);
               
            u1 = u1 + u1_ind_add_R+ u1_ind_add_L;
            u2 = u2 + u2_ind_add_R+ u2_ind_add_L;
            u3 = u3 + u3_ind_add_R+ u3_ind_add_L;
            u4 = u4 + u4_ind_add_R+ u4_ind_add_L;
            v1 = v1 + v1_ind_add_R- v1_ind_add_L;
            v2 = v2 + v2_ind_add_R- v2_ind_add_L;
            v3 = v3 + v3_ind_add_R- v3_ind_add_L;
            v4 = v4 + v4_ind_add_R- v4_ind_add_L;
            w1 = w1 + w1_ind_add_R+ w1_ind_add_L;
            w2 = w2 + w2_ind_add_R+ w2_ind_add_L;
            w3 = w3 + w3_ind_add_R+ w3_ind_add_L;
            w4 = w4 + w4_ind_add_R+ w4_ind_add_L;
           end
        end
        
        % Loop for inducted velocity of wake vortex rings
        r = 1;
        
        while r < now_index
           for s = 1:n_cols
               
               [u1w_ind_add_R,v1w_ind_add_R,w1w_ind_add_R] = func_voring_comp(X_1,Y_1,Z_1,wake_vortex_ring_now{r,s},...
                   wake_vortex_ring_now{r,s}.GAMMA,0);
               [u2w_ind_add_R,v2w_ind_add_R,w2w_ind_add_R] = func_voring_comp(X_2,Y_2,Z_2,wake_vortex_ring_now{r,s},...
                   wake_vortex_ring_now{r,s}.GAMMA,0);
               [u3w_ind_add_R,v3w_ind_add_R,w3w_ind_add_R] = func_voring_comp(X_3,Y_3,Z_3,wake_vortex_ring_now{r,s},...
                   wake_vortex_ring_now{r,s}.GAMMA,0);
               [u4w_ind_add_R,v4w_ind_add_R,w4w_ind_add_R] = func_voring_comp(X_4,Y_4,Z_4,wake_vortex_ring_now{r,s},...
                   wake_vortex_ring_now{r,s}.GAMMA,0);
               [u1w_ind_add_L,v1w_ind_add_L,w1w_ind_add_L] = func_voring_comp(X_1,-Y_1,Z_1,wake_vortex_ring_now{r,s},...
                   wake_vortex_ring_now{r,s}.GAMMA,0);
               [u2w_ind_add_L,v2w_ind_add_L,w2w_ind_add_L] = func_voring_comp(X_2,-Y_2,Z_2,wake_vortex_ring_now{r,s},...
                   wake_vortex_ring_now{r,s}.GAMMA,0);
               [u3w_ind_add_L,v3w_ind_add_L,w3w_ind_add_L] = func_voring_comp(X_3,-Y_3,Z_3,wake_vortex_ring_now{r,s},...
                   wake_vortex_ring_now{r,s}.GAMMA,0);
               [u4w_ind_add_L,v4w_ind_add_L,w4w_ind_add_L] = func_voring_comp(X_4,-Y_4,Z_4,wake_vortex_ring_now{r,s},...
                   wake_vortex_ring_now{r,s}.GAMMA,0);
            u1w = u1w + u1w_ind_add_R+ u1w_ind_add_L;
            u2w = u2w + u2w_ind_add_R+ u2w_ind_add_L;
            u3w = u3w + u3w_ind_add_R+ u3w_ind_add_L;
            u4w = u4w + u4w_ind_add_R+ u4w_ind_add_L;
            v1w = v1w + v1w_ind_add_R- v1w_ind_add_L;
            v2w = v2w + v2w_ind_add_R- v2w_ind_add_L;
            v3w = v3w + v3w_ind_add_R- v3w_ind_add_L;
            v4w = v4w + v4w_ind_add_R- v4w_ind_add_L;
            w1w = w1w + w1w_ind_add_R+ w1w_ind_add_L;
            w2w = w2w + w2w_ind_add_R+ w2w_ind_add_L;
            w3w = w3w + w3w_ind_add_R+ w3w_ind_add_L;
            w4w = w4w + w4w_ind_add_R+ w4w_ind_add_L;
           end
           r = r + 1;
        end
        
        
        
        
        dX_1 = (u1+u1w+u1w_LE)*dt;
        dX_2 = (u2+u2w+u2w_LE)*dt;
        dX_3 = (u3+u3w+u3w_LE)*dt;
        dX_4 = (u4+u4w+u4w_LE)*dt;
        dY_1 = (v1+v1w+v1w_LE)*dt;
        dY_2 = (v2+v2w+v2w_LE)*dt;
        dY_3 = (v3+v3w+v3w_LE)*dt;
        dY_4 = (v4+v4w+v4w_LE)*dt;
        dZ_1 = (w1+w1w+w1w_LE)*dt;
        dZ_2 = (w2+w2w+w2w_LE)*dt;
        dZ_3 = (w3+w3w+w3w_LE)*dt;
        dZ_4 = (w4+w4w+w4w_LE)*dt;
        
      
        
        
        
        NW_wake_vortex_rings{i,j}.X_1 = wake_vortex_ring_now{i,j}.X_1 + dX_1 ;
        NW_wake_vortex_rings{i,j}.X_2 = wake_vortex_ring_now{i,j}.X_2 + dX_2 ;
        NW_wake_vortex_rings{i,j}.X_3 = wake_vortex_ring_now{i,j}.X_3 + dX_3 ;
        NW_wake_vortex_rings{i,j}.X_4 = wake_vortex_ring_now{i,j}.X_4 + dX_4 ;
        NW_wake_vortex_rings{i,j}.Y_1 = wake_vortex_ring_now{i,j}.Y_1 + dY_1 ;
        NW_wake_vortex_rings{i,j}.Y_2 = wake_vortex_ring_now{i,j}.Y_2 + dY_2 ;
        NW_wake_vortex_rings{i,j}.Y_3 = wake_vortex_ring_now{i,j}.Y_3 + dY_3 ;
        NW_wake_vortex_rings{i,j}.Y_4 = wake_vortex_ring_now{i,j}.Y_4 + dY_4 ;
        NW_wake_vortex_rings{i,j}.Z_1 = wake_vortex_ring_now{i,j}.Z_1 + dZ_1 ;
        NW_wake_vortex_rings{i,j}.Z_2 = wake_vortex_ring_now{i,j}.Z_2 + dZ_2 ;
        NW_wake_vortex_rings{i,j}.Z_3 = wake_vortex_ring_now{i,j}.Z_3 + dZ_3 ;
        NW_wake_vortex_rings{i,j}.Z_4 = wake_vortex_ring_now{i,j}.Z_4 + dZ_4 ;
        
        
        else
            
        X_2 = NW_wake_vortex_rings{i,j}.X_2;
        X_4 = NW_wake_vortex_rings{i,j}.X_4;
        Y_2 = NW_wake_vortex_rings{i,j}.Y_2;
        Y_4 = NW_wake_vortex_rings{i,j}.Y_4;
        Z_2 = NW_wake_vortex_rings{i,j}.Z_2;
        Z_4 = NW_wake_vortex_rings{i,j}.Z_4;
         
        % Loop for inducted velocity of bound vortex rings
u2 = 0;
u4 = 0;
v2 = 0;
v4 = 0;
w2 = 0;
w4 = 0;
u2w = 0;
u4w = 0;
v2w = 0;
v4w = 0;
w2w = 0;
w4w = 0;
u2w_LE = 0;
u4w_LE = 0;
v2w_LE = 0;
v4w_LE = 0;
w2w_LE = 0;
w4w_LE = 0;
        for l = 1:n_rows
           for m = 1:n_cols
               
               [u2_ind_add_R,v2_ind_add_R,w2_ind_add_R] = func_voring_comp(X_2,Y_2,Z_2,vortices_CELL{l,m},...
                   gamma_mat(l,m),0);
               [u4_ind_add_R,v4_ind_add_R,w4_ind_add_R] = func_voring_comp(X_4,Y_4,Z_4,vortices_CELL{l,m},...
                   gamma_mat(l,m),0);
               [u2_ind_add_L,v2_ind_add_L,w2_ind_add_L] = func_voring_comp(X_2,-Y_2,Z_2,vortices_CELL{l,m},...
                   gamma_mat(l,m),0);
               [u4_ind_add_L,v4_ind_add_L,w4_ind_add_L] = func_voring_comp(X_4,-Y_4,Z_4,vortices_CELL{l,m},...
                   gamma_mat(l,m),0);
               
            u2 = u2 + u2_ind_add_R+ u2_ind_add_L;
            u4 = u4 + u4_ind_add_R+ u4_ind_add_L;
            v2 = v2 + v2_ind_add_R- v2_ind_add_L;
            v4 = v4 + v4_ind_add_R- v4_ind_add_L;
            w2 = w2 + w2_ind_add_R+ w2_ind_add_L;
            w4 = w4 + w4_ind_add_R+ w4_ind_add_L;
           end
        end
        
        % Loop for inducted velocity of wake vortex rings
        r = 1;
        
        while r < now_index
           for s = 1:n_cols
               
               [u2w_ind_add_R,v2w_ind_add_R,w2w_ind_add_R] = func_voring_comp(X_2,Y_2,Z_2,wake_vortex_ring_now{r,s},...
                   wake_vortex_ring_now{r,s}.GAMMA,0);
               [u4w_ind_add_R,v4w_ind_add_R,w4w_ind_add_R] = func_voring_comp(X_4,Y_4,Z_4,wake_vortex_ring_now{r,s},...
                   wake_vortex_ring_now{r,s}.GAMMA,0);
               [u2w_ind_add_L,v2w_ind_add_L,w2w_ind_add_L] = func_voring_comp(X_2,-Y_2,Z_2,wake_vortex_ring_now{r,s},...
                   wake_vortex_ring_now{r,s}.GAMMA,0);
               [u4w_ind_add_L,v4w_ind_add_L,w4w_ind_add_L] = func_voring_comp(X_4,-Y_4,Z_4,wake_vortex_ring_now{r,s},...
                   wake_vortex_ring_now{r,s}.GAMMA,0);
            
            u2w = u2w + u2w_ind_add_R+ u2w_ind_add_L;
            u4w = u4w + u4w_ind_add_R+ u4w_ind_add_L;
            v2w = v2w + v2w_ind_add_R- v2w_ind_add_L;
            v4w = v4w + v4w_ind_add_R- v4w_ind_add_L;
            w2w = w2w + w2w_ind_add_R+ w2w_ind_add_L;
            w4w = w4w + w4w_ind_add_R+ w4w_ind_add_L;
           end
           r = r + 1;
        end
        
        
        
        dX_2 = (u2+u2w+u2w_LE)*dt;
        dX_4 = (u4+u4w+u4w_LE)*dt;
        dY_2 = (v2+v2w+v2w_LE)*dt;
        dY_4 = (v4+v4w+v4w_LE)*dt;
        dZ_2 = (w2+w2w+w2w_LE)*dt;
        dZ_4 = (w4+w4w+w4w_LE)*dt;
        
      
        
        
        
        NW_wake_vortex_rings{i,j}.X_1 = NW_wake_vortex_rings{i,j-1}.X_2 ;
        NW_wake_vortex_rings{i,j}.X_2 = wake_vortex_ring_now{i,j}.X_2 + dX_2 ;
        NW_wake_vortex_rings{i,j}.X_3 = NW_wake_vortex_rings{i,j-1}.X_4 ;
        NW_wake_vortex_rings{i,j}.X_4 = wake_vortex_ring_now{i,j}.X_4 + dX_4 ;
        NW_wake_vortex_rings{i,j}.Y_1 = NW_wake_vortex_rings{i,j-1}.Y_2 ;
        NW_wake_vortex_rings{i,j}.Y_2 = wake_vortex_ring_now{i,j}.Y_2 + dY_2 ;
        NW_wake_vortex_rings{i,j}.Y_3 = NW_wake_vortex_rings{i,j-1}.Y_4 ;
        NW_wake_vortex_rings{i,j}.Y_4 = wake_vortex_ring_now{i,j}.Y_4 + dY_4 ;
        NW_wake_vortex_rings{i,j}.Z_1 = NW_wake_vortex_rings{i,j-1}.Z_2 ;
        NW_wake_vortex_rings{i,j}.Z_2 = wake_vortex_ring_now{i,j}.Z_2 + dZ_2 ;
        NW_wake_vortex_rings{i,j}.Z_3 = NW_wake_vortex_rings{i,j-1}.Z_4 ;
        NW_wake_vortex_rings{i,j}.Z_4 = wake_vortex_ring_now{i,j}.Z_4 + dZ_4 ;    
            
        end
    
        
    end
    else
     for j = 1:n_cols
        if j == 1
        X_1 = NW_wake_vortex_rings{i,j}.X_1;
        X_2 = NW_wake_vortex_rings{i,j}.X_2;
        Y_1 = NW_wake_vortex_rings{i,j}.Y_1;
        Y_2 = NW_wake_vortex_rings{i,j}.Y_2;
        Z_1 = NW_wake_vortex_rings{i,j}.Z_1;
        Z_2 = NW_wake_vortex_rings{i,j}.Z_2;
    
         
        % Loop for inducted velocity of bound vortex rings
u1 = 0;
u2 = 0;
v1 = 0;
v2 = 0;
w1 = 0;
w2 = 0;
u1w = 0;
u2w = 0;
v1w = 0;
v2w = 0;
w1w = 0;
w2w = 0;
u1w_LE = 0;
u2w_LE = 0;
v1w_LE = 0;
v2w_LE = 0;
w1w_LE = 0;
w2w_LE = 0;
        for l = 1:n_rows
           for m = 1:n_cols
               
               [u1_ind_add_R,v1_ind_add_R,w1_ind_add_R] = func_voring_comp(X_1,Y_1,Z_1,vortices_CELL{l,m},...
                   gamma_mat(l,m),0);
               [u2_ind_add_R,v2_ind_add_R,w2_ind_add_R] = func_voring_comp(X_2,Y_2,Z_2,vortices_CELL{l,m},...
                   gamma_mat(l,m),0);
               [u1_ind_add_L,v1_ind_add_L,w1_ind_add_L] = func_voring_comp(X_1,-Y_1,Z_1,vortices_CELL{l,m},...
                   gamma_mat(l,m),0);
               [u2_ind_add_L,v2_ind_add_L,w2_ind_add_L] = func_voring_comp(X_2,-Y_2,Z_2,vortices_CELL{l,m},...
                   gamma_mat(l,m),0);

               
            u1 = u1 + u1_ind_add_R+ u1_ind_add_L;
            u2 = u2 + u2_ind_add_R+ u2_ind_add_L;
            v1 = v1 + v1_ind_add_R- v1_ind_add_L;
            v2 = v2 + v2_ind_add_R- v2_ind_add_L;
            w1 = w1 + w1_ind_add_R+ w1_ind_add_L;
            w2 = w2 + w2_ind_add_R+ w2_ind_add_L;
           end
        end
        
        % Loop for inducted velocity of wake vortex rings
        r = 1;
        
        while r < now_index
           for s = 1:n_cols
               
               [u1w_ind_add_R,v1w_ind_add_R,w1w_ind_add_R] = func_voring_comp(X_1,Y_1,Z_1,wake_vortex_ring_now{r,s},...
                   wake_vortex_ring_now{r,s}.GAMMA,0);
               [u2w_ind_add_R,v2w_ind_add_R,w2w_ind_add_R] = func_voring_comp(X_2,Y_2,Z_2,wake_vortex_ring_now{r,s},...
                   wake_vortex_ring_now{r,s}.GAMMA,0);
               [u1w_ind_add_L,v1w_ind_add_L,w1w_ind_add_L] = func_voring_comp(X_1,-Y_1,Z_1,wake_vortex_ring_now{r,s},...
                   wake_vortex_ring_now{r,s}.GAMMA,0);
               [u2w_ind_add_L,v2w_ind_add_L,w2w_ind_add_L] = func_voring_comp(X_2,-Y_2,Z_2,wake_vortex_ring_now{r,s},...
                   wake_vortex_ring_now{r,s}.GAMMA,0);

            u1w = u1w + u1w_ind_add_R+ u1w_ind_add_L;
            u2w = u2w + u2w_ind_add_R+ u2w_ind_add_L;
            v1w = v1w + v1w_ind_add_R- v1w_ind_add_L;
            v2w = v2w + v2w_ind_add_R- v2w_ind_add_L;
            w1w = w1w + w1w_ind_add_R+ w1w_ind_add_L;
            w2w = w2w + w2w_ind_add_R+ w2w_ind_add_L;

           end
           r = r + 1;
        end
        
        
        
        
        dX_1 = (u1+u1w+u1w_LE)*dt;
        dX_2 = (u2+u2w+u2w_LE)*dt;
        dY_1 = (v1+v1w+v1w_LE)*dt;
        dY_2 = (v2+v2w+v2w_LE)*dt;
        dZ_1 = (w1+w1w+w1w_LE)*dt;
        dZ_2 = (w2+w2w+w2w_LE)*dt;

        
        
        NW_wake_vortex_rings{i,j}.X_1 = wake_vortex_ring_now{i,j}.X_1 + dX_1 ;
        NW_wake_vortex_rings{i,j}.X_2 = wake_vortex_ring_now{i,j}.X_2 + dX_2 ;
        NW_wake_vortex_rings{i,j}.X_3 = NW_wake_vortex_rings{i-1,j}.X_1 ;
        NW_wake_vortex_rings{i,j}.X_4 = NW_wake_vortex_rings{i-1,j}.X_2 ;
        NW_wake_vortex_rings{i,j}.Y_1 = wake_vortex_ring_now{i,j}.Y_1 + dY_1 ;
        NW_wake_vortex_rings{i,j}.Y_2 = wake_vortex_ring_now{i,j}.Y_2 + dY_2 ;
        NW_wake_vortex_rings{i,j}.Y_3 = NW_wake_vortex_rings{i-1,j}.Y_1 ;
        NW_wake_vortex_rings{i,j}.Y_4 = NW_wake_vortex_rings{i-1,j}.Y_2 ;
        NW_wake_vortex_rings{i,j}.Z_1 = wake_vortex_ring_now{i,j}.Z_1 + dZ_1 ;
        NW_wake_vortex_rings{i,j}.Z_2 = wake_vortex_ring_now{i,j}.Z_2 + dZ_2 ;
        NW_wake_vortex_rings{i,j}.Z_3 = NW_wake_vortex_rings{i-1,j}.Z_1 ;
        NW_wake_vortex_rings{i,j}.Z_4 = NW_wake_vortex_rings{i-1,j}.Z_2 ;
        
        
        else
            
        X_2 = NW_wake_vortex_rings{i,j}.X_2;
        Y_2 = NW_wake_vortex_rings{i,j}.Y_2;
        Z_2 = NW_wake_vortex_rings{i,j}.Z_2;
         
        % Loop for inducted velocity of bound vortex rings
u2 = 0;
v2 = 0;
w2 = 0;
u2w = 0;
v2w = 0;
w2w = 0;
u2w_LE = 0;
v2w_LE = 0;
w2w_LE = 0;
        for l = 1:n_rows
           for m = 1:n_cols
               
               [u2_ind_add_R,v2_ind_add_R,w2_ind_add_R] = func_voring_comp(X_2,Y_2,Z_2,vortices_CELL{l,m},...
                   gamma_mat(l,m),0);
               [u2_ind_add_L,v2_ind_add_L,w2_ind_add_L] = func_voring_comp(X_2,-Y_2,Z_2,vortices_CELL{l,m},...
                   gamma_mat(l,m),0);

               
            u2 = u2 + u2_ind_add_R+ u2_ind_add_L;
            v2 = v2 + v2_ind_add_R- v2_ind_add_L;
            w2 = w2 + w2_ind_add_R+ w2_ind_add_L;
            
           end
        end
        
        % Loop for inducted velocity of wake vortex rings
        r = 1;
        
        while r < now_index
           for s = 1:n_cols
               
               [u2w_ind_add_R,v2w_ind_add_R,w2w_ind_add_R] = func_voring_comp(X_2,Y_2,Z_2,wake_vortex_ring_now{r,s},...
                   wake_vortex_ring_now{r,s}.GAMMA,0);
               [u2w_ind_add_L,v2w_ind_add_L,w2w_ind_add_L] = func_voring_comp(X_2,-Y_2,Z_2,wake_vortex_ring_now{r,s},...
                   wake_vortex_ring_now{r,s}.GAMMA,0);

            
            u2w = u2w + u2w_ind_add_R+ u2w_ind_add_L;
            v2w = v2w + v2w_ind_add_R- v2w_ind_add_L;
            w2w = w2w + w2w_ind_add_R+ w2w_ind_add_L;
           end
           r = r + 1;
        end
        
        
        
        dX_2 = (u2+u2w+u2w_LE)*dt;
        dY_2 = (v2+v2w+v2w_LE)*dt;
        dZ_2 = (w2+w2w+w2w_LE)*dt;
        
      
        
        
        
        NW_wake_vortex_rings{i,j}.X_1 = NW_wake_vortex_rings{i,j-1}.X_2 ;
        NW_wake_vortex_rings{i,j}.X_2 = wake_vortex_ring_now{i,j}.X_2 + dX_2 ;
        NW_wake_vortex_rings{i,j}.X_3 = NW_wake_vortex_rings{i,j-1}.X_4 ;
        NW_wake_vortex_rings{i,j}.X_4 = NW_wake_vortex_rings{i-1,j}.X_2 ;
        NW_wake_vortex_rings{i,j}.Y_1 = NW_wake_vortex_rings{i,j-1}.Y_2 ;
        NW_wake_vortex_rings{i,j}.Y_2 = wake_vortex_ring_now{i,j}.Y_2 + dY_2 ;
        NW_wake_vortex_rings{i,j}.Y_3 = NW_wake_vortex_rings{i,j-1}.Y_4 ;
        NW_wake_vortex_rings{i,j}.Y_4 = NW_wake_vortex_rings{i-1,j}.Y_2 ;
        NW_wake_vortex_rings{i,j}.Z_1 = NW_wake_vortex_rings{i,j-1}.Z_2 ;
        NW_wake_vortex_rings{i,j}.Z_2 = wake_vortex_ring_now{i,j}.Z_2 + dZ_2 ;
        NW_wake_vortex_rings{i,j}.Z_3 = NW_wake_vortex_rings{i,j-1}.Z_4 ;
        NW_wake_vortex_rings{i,j}.Z_4 = NW_wake_vortex_rings{i-1,j}.Z_2 ;    
            
        end
    
        
    end    
    end
  else
  end
i = i + 1;   
end

end