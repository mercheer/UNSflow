function  vortices_mat = func_upd_pos(t,vortices_mat_0,Quat,flag_BENDMOT)


global n_rows n_cols

%% Adding paths

addpath([pwd,'/geometry']);
addpath([pwd,'/functions']);
addpath([pwd,'/data']);


%% To load panels

load([pwd,'/data/data_refs.mat']);
load([pwd,'/data/data_kinematics.mat']);

%% To update panels position in Inertial reference system

rows = n_rows;
cols = n_cols;

vortices_mat_buff = vortices_mat_0;
vortices_mat_buff_cycle = cell(1,rows*cols);
vortices_mat_0_cycle    = cell(1,rows*cols);

k = 1;

for i = 1:rows
    for j = 1:cols
        vortices_mat_buff_cycle{k} = vortices_mat_buff{i,j};
        vortices_mat_0_cycle{k}    = vortices_mat_0{i,j};
    k = k + 1;
    end
end

N = rows*cols;

if flag_BENDMOT == 1
    
Z_BEND = Z_BEND;

else
    
Z_BEND = 0;

end


parfor cyc = 1:N
    
%            N_x       = vortices_mat_0_cycle{cyc}.N_x;
%            N_y       = vortices_mat_0_cycle{cyc}.N_y;
%            N_z       = vortices_mat_0_cycle{cyc}.N_z;
%            
%            N_v = transpose(quat2dcm(Quat))*[N_x;N_y;N_z];
%            
%            vortices_mat_buff_cycle{cyc}.N_x = N_v(1);
%            vortices_mat_buff_cycle{cyc}.N_y = N_v(2);
%            vortices_mat_buff_cycle{cyc}.N_z = N_v(3);
        
        [vortices_mat_buff_cycle{cyc}.X_1,vortices_mat_buff_cycle{cyc}.Y_1,vortices_mat_buff_cycle{cyc}.Z_1] = ...
            func_kinematics_pos(t,vortices_mat_0_cycle{cyc}.X_1,vortices_mat_0_cycle{cyc}.Y_1,...
            vortices_mat_0_cycle{cyc}.Z_1,Quat);
        [vortices_mat_buff_cycle{cyc}.X_2,vortices_mat_buff_cycle{cyc}.Y_2,vortices_mat_buff_cycle{cyc}.Z_2] = ...
            func_kinematics_pos(t,vortices_mat_0_cycle{cyc}.X_2,vortices_mat_0_cycle{cyc}.Y_2,...
            vortices_mat_0_cycle{cyc}.Z_2,Quat);
        [vortices_mat_buff_cycle{cyc}.X_3,vortices_mat_buff_cycle{cyc}.Y_3,vortices_mat_buff_cycle{cyc}.Z_3] = ...
            func_kinematics_pos(t,vortices_mat_0_cycle{cyc}.X_3,vortices_mat_0_cycle{cyc}.Y_3,...
            vortices_mat_0_cycle{cyc}.Z_3,Quat);
        [vortices_mat_buff_cycle{cyc}.X_4,vortices_mat_buff_cycle{cyc}.Y_4,vortices_mat_buff_cycle{cyc}.Z_4] = ...
            func_kinematics_pos(t,vortices_mat_0_cycle{cyc}.X_4,vortices_mat_0_cycle{cyc}.Y_4,...
            vortices_mat_0_cycle{cyc}.Z_4,Quat);
        [vortices_mat_buff_cycle{cyc}.X_C,vortices_mat_buff_cycle{cyc}.Y_C,vortices_mat_buff_cycle{cyc}.Z_C] = ...
            func_kinematics_pos(t,vortices_mat_0_cycle{cyc}.X_C,vortices_mat_0_cycle{cyc}.Y_C,...
            vortices_mat_0_cycle{cyc}.Z_C,Quat);
        
        %% BENDING MOTION
        
        if flag_BENDMOT == 1
        
        vortices_mat_buff_cycle{cyc}.Z_1 = vortices_mat_buff_cycle{cyc}.Z_1 + Z_BEND(t,vortices_mat_0_cycle{cyc}.Y_1);
        vortices_mat_buff_cycle{cyc}.Z_2 = vortices_mat_buff_cycle{cyc}.Z_2 + Z_BEND(t,vortices_mat_0_cycle{cyc}.Y_2);
        vortices_mat_buff_cycle{cyc}.Z_3 = vortices_mat_buff_cycle{cyc}.Z_3 + Z_BEND(t,vortices_mat_0_cycle{cyc}.Y_3);
        vortices_mat_buff_cycle{cyc}.Z_4 = vortices_mat_buff_cycle{cyc}.Z_4 + Z_BEND(t,vortices_mat_0_cycle{cyc}.Y_4);
        vortices_mat_buff_cycle{cyc}.Z_C = vortices_mat_buff_cycle{cyc}.Z_C + Z_BEND(t,vortices_mat_0_cycle{cyc}.Y_C);
        
        else
        end
        %% UPDATE NORMALS USING CROSS-PRODUCT OF DIAGONAL VECTORS
        
        A1 = vortices_mat_buff_cycle{cyc}.X_3 - vortices_mat_buff_cycle{cyc}.X_2;
        A2 = vortices_mat_buff_cycle{cyc}.Y_3 - vortices_mat_buff_cycle{cyc}.Y_2;
        A3 = vortices_mat_buff_cycle{cyc}.Z_3 - vortices_mat_buff_cycle{cyc}.Z_2;
        B1 = vortices_mat_buff_cycle{cyc}.X_4 - vortices_mat_buff_cycle{cyc}.X_1;
        B2 = vortices_mat_buff_cycle{cyc}.Y_4 - vortices_mat_buff_cycle{cyc}.Y_1;
        B3 = vortices_mat_buff_cycle{cyc}.Z_4 - vortices_mat_buff_cycle{cyc}.Z_1;
        XX = A2*B3 - A3*B2;
        YY = B1*A3 - A1*B3;
        ZZ = A1*B2 - A2*B1;
        A = sqrt(XX^2+YY^2+ZZ^2);
        
        vortices_mat_buff_cycle{cyc}.N_x = XX/A;
        vortices_mat_buff_cycle{cyc}.N_y = YY/A;
        vortices_mat_buff_cycle{cyc}.N_z = ZZ/A;

end

k = 1;

for i = 1:rows
    for j = 1:cols
        vortices_mat_buff{i,j} = vortices_mat_buff_cycle{k};
    k = k + 1;
    end
end



vortices_mat = vortices_mat_buff;

end
