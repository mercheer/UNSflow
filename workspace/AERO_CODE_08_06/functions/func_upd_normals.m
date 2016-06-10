function  func_upd_INFO()


global n_rows n_cols vortices_mat



%% Adding paths

addpath([pwd,'/geometry']);
addpath([pwd,'/functions']);
addpath([pwd,'/data']);


%% Load new OUTPUT.dat

[colheaders,data]=func_import_OUTPUT('OUTPUT.dat');

for i = 1:size(colheaders,2)
    assignin('caller' , genvarname(colheaders{i}),data(:,i));
end

%% To update panels position in Inertial reference system

rows = n_rows;
cols = n_cols;

k = 1;

for j = 1:cols
    for i = 1:rows
	
           vortices_mat{i,j}.N_x = -data(k,4);
           vortices_mat{i,j}.N_y = data(k,3);
           vortices_mat{i,j}.N_z = data(k,5);
           
    k = k + 1;
        
    end
end

end
