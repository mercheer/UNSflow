function func_LE_shedwake(i,vortices_mat_CELL)
% Add a row of vortex rings to the wake

global LE_NW_wake_vortex_rings n_cols n_rows

%% ADDPATH

addpath([pwd,'/functions']);
addpath([pwd,'/data']);

%% To load vortex lattice mesh

LE_wake_vortex_ring = struct('X_1',[],'Y_1',[],'Z_1',[],...
                     'X_2',[],'Y_2',[],'Z_2',[],...
                     'X_3',[],'Y_3',[],'Z_3',[],...
                     'X_4',[],'Y_4',[],'Z_4',[],...
                     'X_C',[],'Y_C',[],'Z_C',[],...
                     'GAMMA',[],...
                     'N_x',[],'N_y',[],'N_z',[],...
                     'Z_SIGN',[]);


%% Shedding procedure
                 
if i == 2             
    

for j = 1 : n_cols
       
       LE_NW_wake_vortex_rings{i-1,j} = LE_wake_vortex_ring;
       
       LE_NW_wake_vortex_rings{i-1,j}.X_1 =  vortices_mat_CELL{i}{1,j}.X_1 ;
       LE_NW_wake_vortex_rings{i-1,j}.X_2 =  vortices_mat_CELL{i}{1,j}.X_2 ;
       LE_NW_wake_vortex_rings{i-1,j}.Y_1 =  vortices_mat_CELL{i}{1,j}.Y_1 ;
       LE_NW_wake_vortex_rings{i-1,j}.Y_2 =  vortices_mat_CELL{i}{1,j}.Y_2 ;
       LE_NW_wake_vortex_rings{i-1,j}.Z_1 =  vortices_mat_CELL{i}{1,j}.Z_1 ;
       LE_NW_wake_vortex_rings{i-1,j}.Z_2 =  vortices_mat_CELL{i}{1,j}.Z_2 ;
       
       LE_NW_wake_vortex_rings{i-1,j}.X_3 = vortices_mat_CELL{i-1}{1,j}.X_1 ;
       LE_NW_wake_vortex_rings{i-1,j}.X_4 = vortices_mat_CELL{i-1}{1,j}.X_2 ;
       LE_NW_wake_vortex_rings{i-1,j}.Y_3 = vortices_mat_CELL{i-1}{1,j}.Y_1 ;
       LE_NW_wake_vortex_rings{i-1,j}.Y_4 = vortices_mat_CELL{i-1}{1,j}.Y_2 ;
       LE_NW_wake_vortex_rings{i-1,j}.Z_3 = vortices_mat_CELL{i-1}{1,j}.Z_1 ;
       LE_NW_wake_vortex_rings{i-1,j}.Z_4 = vortices_mat_CELL{i-1}{1,j}.Z_2 ;
       LE_NW_wake_vortex_rings{i-1,j}.Z_SIGN = 0;


end
    
else    

for j = 1 : n_cols
       
       LE_NW_wake_vortex_rings{i-1,j} = LE_wake_vortex_ring;
       
       LE_NW_wake_vortex_rings{i-1,j}.X_1 =  vortices_mat_CELL{i}{1,j}.X_1 ;
       LE_NW_wake_vortex_rings{i-1,j}.X_2 =  vortices_mat_CELL{i}{1,j}.X_2 ;
       LE_NW_wake_vortex_rings{i-1,j}.Y_1 =  vortices_mat_CELL{i}{1,j}.Y_1 ;
       LE_NW_wake_vortex_rings{i-1,j}.Y_2 =  vortices_mat_CELL{i}{1,j}.Y_2 ;
       LE_NW_wake_vortex_rings{i-1,j}.Z_1 =  vortices_mat_CELL{i}{1,j}.Z_1 ;
       LE_NW_wake_vortex_rings{i-1,j}.Z_2 =  vortices_mat_CELL{i}{1,j}.Z_2 ;
       
       LE_NW_wake_vortex_rings{i-1,j}.X_3 =  LE_NW_wake_vortex_rings{i-2,j}.X_1 ;
       LE_NW_wake_vortex_rings{i-1,j}.X_4 =  LE_NW_wake_vortex_rings{i-2,j}.X_2 ;
       LE_NW_wake_vortex_rings{i-1,j}.Y_3 =  LE_NW_wake_vortex_rings{i-2,j}.Y_1 ;
       LE_NW_wake_vortex_rings{i-1,j}.Y_4 =  LE_NW_wake_vortex_rings{i-2,j}.Y_2 ;
       LE_NW_wake_vortex_rings{i-1,j}.Z_3 =  LE_NW_wake_vortex_rings{i-2,j}.Z_1 ;
       LE_NW_wake_vortex_rings{i-1,j}.Z_4 =  LE_NW_wake_vortex_rings{i-2,j}.Z_2 ;
       LE_NW_wake_vortex_rings{i-1,j}.Z_SIGN = 0;


end    

end

