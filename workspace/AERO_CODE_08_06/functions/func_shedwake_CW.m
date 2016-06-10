function func_shedwake_CW(i,vortices_mat_CELL,cw)
% Add a row of vortex rings to the wake

global NW_wake_vortex_rings n_cols n_rows

%% ADDPATH

addpath([pwd,'/functions']);
addpath([pwd,'/data']);

%% To load vortex lattice mesh

wake_vortex_ring = struct('X_1',[],'Y_1',[],'Z_1',[],...
                     'X_2',[],'Y_2',[],'Z_2',[],...
                     'X_3',[],'Y_3',[],'Z_3',[],...
                     'X_4',[],'Y_4',[],'Z_4',[],...
                     'X_C',[],'Y_C',[],'Z_C',[],...
                     'GAMMA',[],...
                     'N_x',[],'N_y',[],'N_z',[],...
                     'Z_SIGN',[]);


%% Shedding procedure
        
    

for j = 1 : n_cols
       
       NW_wake_vortex_rings{i-1,j} = wake_vortex_ring;
       
       NW_wake_vortex_rings{i-1,j}.X_1 =  vortices_mat_CELL{i}{n_rows,j}.X_3 ;
       NW_wake_vortex_rings{i-1,j}.X_2 =  vortices_mat_CELL{i}{n_rows,j}.X_4 ;
       NW_wake_vortex_rings{i-1,j}.Y_1 =  vortices_mat_CELL{i}{n_rows,j}.Y_3 ;
       NW_wake_vortex_rings{i-1,j}.Y_2 =  vortices_mat_CELL{i}{n_rows,j}.Y_4 ;
       NW_wake_vortex_rings{i-1,j}.Z_1 =  vortices_mat_CELL{i}{n_rows,j}.Z_3 ;
       NW_wake_vortex_rings{i-1,j}.Z_2 =  vortices_mat_CELL{i}{n_rows,j}.Z_4 ;
       
       NW_wake_vortex_rings{i-1,j}.X_3 = (vortices_mat_CELL{i}{n_rows,j}.X_3- cw*(vortices_mat_CELL{i}{n_rows,j}.X_3 - vortices_mat_CELL{i-1}{n_rows,j}.X_3)) ;
       NW_wake_vortex_rings{i-1,j}.X_4 = (vortices_mat_CELL{i}{n_rows,j}.X_4- cw*(vortices_mat_CELL{i}{n_rows,j}.X_4 - vortices_mat_CELL{i-1}{n_rows,j}.X_4)) ;
       NW_wake_vortex_rings{i-1,j}.Y_3 = (vortices_mat_CELL{i}{n_rows,j}.Y_3- cw*(vortices_mat_CELL{i}{n_rows,j}.Y_3 - vortices_mat_CELL{i-1}{n_rows,j}.Y_3)) ;
       NW_wake_vortex_rings{i-1,j}.Y_4 = (vortices_mat_CELL{i}{n_rows,j}.Y_4- cw*(vortices_mat_CELL{i}{n_rows,j}.Y_4 - vortices_mat_CELL{i-1}{n_rows,j}.Y_4)) ;
       NW_wake_vortex_rings{i-1,j}.Z_3 = (vortices_mat_CELL{i}{n_rows,j}.Z_3- cw*(vortices_mat_CELL{i}{n_rows,j}.Z_3 - vortices_mat_CELL{i-1}{n_rows,j}.Z_3)) ;
       NW_wake_vortex_rings{i-1,j}.Z_4 = (vortices_mat_CELL{i}{n_rows,j}.Z_4- cw*(vortices_mat_CELL{i}{n_rows,j}.Z_4 - vortices_mat_CELL{i-1}{n_rows,j}.Z_4)) ;
       NW_wake_vortex_rings{i-1,j}.Z_SIGN = 0;


end
    
if i > 2

for j = 1 : n_cols
       
       
       NW_wake_vortex_rings{i-2,j}.X_1 =  NW_wake_vortex_rings{i-1,j}.X_3 ;
       NW_wake_vortex_rings{i-2,j}.X_2 =  NW_wake_vortex_rings{i-1,j}.X_4 ;
       NW_wake_vortex_rings{i-2,j}.Y_1 =  NW_wake_vortex_rings{i-1,j}.Y_3 ;
       NW_wake_vortex_rings{i-2,j}.Y_2 =  NW_wake_vortex_rings{i-1,j}.Y_4 ;
       NW_wake_vortex_rings{i-2,j}.Z_1 =  NW_wake_vortex_rings{i-1,j}.Z_3 ;
       NW_wake_vortex_rings{i-2,j}.Z_2 =  NW_wake_vortex_rings{i-1,j}.Z_4 ;
       


end
else
end

end

