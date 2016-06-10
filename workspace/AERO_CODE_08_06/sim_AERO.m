F_mat = cell(n_rows, n_cols);
for i = 1:1:n_rows
    for j = 1:1:n_cols
        F_mat{i,j} = [1,1,1];
    end
end
pause(1);
save('F_mat.mat', 'F_mat');