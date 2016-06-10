function  FRAME = func_frame(vortices_mat,C_L_mat,i_w,NW_wake_vortex_rings,LE_NW_wake_vortex_rings,FW_wake_vortons,LE_FW_wake_vortons)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

global n_rows n_cols

for i = 1:n_rows
    for j = 1:n_cols
vortex_rings(1,j+n_cols*(i-1)) = vortices_mat(i,j);
    end
end

for i = 1:length(vortex_rings)
    X_1(i) = vortex_rings{i}.X_1;
    X_2(i) = vortex_rings{i}.X_2;
    X_3(i) = vortex_rings{i}.X_3;
    X_4(i) = vortex_rings{i}.X_4;
    Y_1(i) = vortex_rings{i}.Y_1;
    Y_2(i) = vortex_rings{i}.Y_2;
    Y_3(i) = vortex_rings{i}.Y_3;
    Y_4(i) = vortex_rings{i}.Y_4;
    Z_1(i) = vortex_rings{i}.Z_1;
    Z_2(i) = vortex_rings{i}.Z_2;
    Z_3(i) = vortex_rings{i}.Z_3;
    Z_4(i) = vortex_rings{i}.Z_4;
end

C_L_vector = zeros(1,n_rows*n_cols);

for i = 1:n_rows

    for j = 1:n_cols
    
 C_L_vector(j+n_cols*(i-1))=C_L_mat(i,j);

    end
    
end


h = figure(i_w);

if i_w > 1
    
index = i_w;
while index > 1    
for j = 1:n_cols
NW_wake_vortex_rings_plot(j) = NW_wake_vortex_rings(index-1,j);
end
for i = 1:length(NW_wake_vortex_rings_plot)
    X_1_W(i) = NW_wake_vortex_rings_plot{i}.X_1;
    X_2_W(i) = NW_wake_vortex_rings_plot{i}.X_2;
    X_3_W(i) = NW_wake_vortex_rings_plot{i}.X_3;
    X_4_W(i) = NW_wake_vortex_rings_plot{i}.X_4;
    Y_1_W(i) = NW_wake_vortex_rings_plot{i}.Y_1;
    Y_2_W(i) = NW_wake_vortex_rings_plot{i}.Y_2;
    Y_3_W(i) = NW_wake_vortex_rings_plot{i}.Y_3;
    Y_4_W(i) = NW_wake_vortex_rings_plot{i}.Y_4;
    Z_1_W(i) = NW_wake_vortex_rings_plot{i}.Z_1;
    Z_2_W(i) = NW_wake_vortex_rings_plot{i}.Z_2;
    Z_3_W(i) = NW_wake_vortex_rings_plot{i}.Z_3;
    Z_4_W(i) = NW_wake_vortex_rings_plot{i}.Z_4;
end
for j = 1:n_cols
LE_NW_wake_vortex_rings_plot(j) = LE_NW_wake_vortex_rings(index-1,j);
end
for i = 1:length(NW_wake_vortex_rings_plot)
    LE_X_1_W(i) = LE_NW_wake_vortex_rings_plot{i}.X_1;
    LE_X_2_W(i) = LE_NW_wake_vortex_rings_plot{i}.X_2;
    LE_X_3_W(i) = LE_NW_wake_vortex_rings_plot{i}.X_3;
    LE_X_4_W(i) = LE_NW_wake_vortex_rings_plot{i}.X_4;
    LE_Y_1_W(i) = LE_NW_wake_vortex_rings_plot{i}.Y_1;
    LE_Y_2_W(i) = LE_NW_wake_vortex_rings_plot{i}.Y_2;
    LE_Y_3_W(i) = LE_NW_wake_vortex_rings_plot{i}.Y_3;
    LE_Y_4_W(i) = LE_NW_wake_vortex_rings_plot{i}.Y_4;
    LE_Z_1_W(i) = LE_NW_wake_vortex_rings_plot{i}.Z_1;
    LE_Z_2_W(i) = LE_NW_wake_vortex_rings_plot{i}.Z_2;
    LE_Z_3_W(i) = LE_NW_wake_vortex_rings_plot{i}.Z_3;
    LE_Z_4_W(i) = LE_NW_wake_vortex_rings_plot{i}.Z_4;
end

patch([LE_X_1_W;LE_X_2_W;LE_X_4_W;LE_X_3_W],[LE_Y_1_W;LE_Y_2_W;LE_Y_4_W;LE_Y_3_W],[LE_Z_1_W;LE_Z_2_W;LE_Z_4_W;LE_Z_3_W],[1,1,1]);hold on;
patch([LE_X_1_W;LE_X_2_W;LE_X_4_W;LE_X_3_W],[-LE_Y_1_W;-LE_Y_2_W;-LE_Y_4_W;-LE_Y_3_W],[LE_Z_1_W;LE_Z_2_W;LE_Z_4_W;LE_Z_3_W],[1,1,1]);hold on;
patch([X_1_W;X_2_W;X_4_W;X_3_W],[Y_1_W;Y_2_W;Y_4_W;Y_3_W],[Z_1_W;Z_2_W;Z_4_W;Z_3_W],[1,1,1]);hold on;
patch([X_1_W;X_2_W;X_4_W;X_3_W],[-Y_1_W;-Y_2_W;-Y_4_W;-Y_3_W],[Z_1_W;Z_2_W;Z_4_W;Z_3_W],[1,1,1]);hold on;
xlim([-50 2]);ylim([-5 5]);zlim([-2 2]);
colormap hot
map = colormap;
view(3);
axis equal;

index = index - 1;
end

patch([X_1;X_2;X_4;X_3],[Y_1;Y_2;Y_4;Y_3],[Z_1;Z_2;Z_4;Z_3],C_L_vector);hold on;
patch([X_1;X_2;X_4;X_3],[-Y_1;-Y_2;-Y_4;-Y_3],[Z_1;Z_2;Z_4;Z_3],C_L_vector);hold on;
xlim([-50 2]);ylim([-5 5]);zlim([-2 2]);
colormap hot
map = colormap;
view(3);
axis equal;


else
    
patch([X_1;X_2;X_4;X_3],[Y_1;Y_2;Y_4;Y_3],[Z_1;Z_2;Z_4;Z_3],C_L_vector);hold on;
patch([X_1;X_2;X_4;X_3],[-Y_1;-Y_2;-Y_4;-Y_3],[Z_1;Z_2;Z_4;Z_3],C_L_vector);hold on;
xlim([-50 2]);ylim([-5 5]);zlim([-2 2]);
colormap hot
map = colormap;
view(3);
axis equal;


end

FRAME = getframe;
close(h);



end

