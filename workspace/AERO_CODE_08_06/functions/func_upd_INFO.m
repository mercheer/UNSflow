function  [UO,VO,WO] = func_upd_INFO()


global n_rows n_cols vortices_mat



%% Adding paths

addpath([pwd,'/geometry']);
addpath([pwd,'/functions']);
addpath([pwd,'/data']);


%% Load new OUTPUT.dat

[colheaders1,data1]=func_import_OUTPUT('OUTPUT.dat');
[colheaders2,data2]=func_import_OUTPUT('OUTPUTVEL.dat');


for i = 1:size(colheaders1,2)
    assignin('caller' , genvarname(colheaders1{i}),data1(:,i));
end

for i = 1:size(colheaders2,2)
    assignin('caller' , genvarname(colheaders2{i}),data2(:,i));
end


%% To update panels position in Inertial reference system

rows = n_rows;
cols = n_cols;

U0 = zeros(rows,cols);
V0 = zeros(rows,cols);
W0 = zeros(rows,cols);


k = 1;

for j = 1:cols
    for i = 1:rows
	
           vortices_mat{i,j}.N_x = -data1(k,4);
           vortices_mat{i,j}.N_y = data1(k,3);
           vortices_mat{i,j}.N_z = data1(k,5);
           UO(i,j) = -data2(k,4);
           VO(i,j) =  data2(k,3);
           WO(i,j) =  data2(k,5);
    k = k + 1;
        
    end
end

end
