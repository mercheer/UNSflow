%% FORM VORTEX RING LATTICES UPON GIVEN PATCHES

clearvars -except i_SS alfa_SS C_L_steady_state C_Di_steady_state C_L_steady_state_conv C_L_central_panel TIME_FIN;


%% Adding paths

addpath([pwd,'/functions']);
addpath([pwd,'/data']);

%% Import wing patches

load([pwd,'/data/data_patches.mat']);
load([pwd,'/data/data_refs.mat']);

% A patch contains these informations:
% 1 - Leftmost X-coordinate of L.E
% 2 - Leftmost Y-coordinate of L.E
% 3 - Leftmost Z-coordinate of L.E
% 4 - Rightmost X-coordinate of L.E
% 5 - Rightmost Y-coordinate of L.E
% 6 - Rightmost Z-coordinate of L.E
% 7- Leftmost chord-length
% 8 - Rightmost chord length
% 9 - Twist angle
%
% Patches are stored from the leftmost tip to the rightmost tip

N_p = length(patches);

%% To choose numbers of spanwise and chordwise elements for each patch and Wake Position

% N_span = zeros(1,length(patches));
% 
% for i = 1:length(patches)
%     N_span(i) = 1;
% end

N_span = 24;        % spanwise divisions,need the number of divisions first

M = 4;               % chordwise divisions


%% Vortex lattice and rings pre-allocation

vortex_ring = struct('X_1',[],'Y_1',[],'Z_1',[],...
                     'X_2',[],'Y_2',[],'Z_2',[],...
                     'X_3',[],'Y_3',[],'Z_3',[],...
                     'X_4',[],'Y_4',[],'Z_4',[],...
                     'X_C',[],'Y_C',[],'Z_C',[],...
                     'S',[],...
                     'N_x',[],'N_y',[],'N_z',[], ...
                     'Z_SIGN',[], ...
                     'TWIST',[]);
                 
vortex_lattice = cell(1,N_p);


%% Vortex lattice creation 

for i=1:N_p
    
patch_index = i;
N =N_span(i) ;
tic;
vortex_lattice{1,i} = func_panel(M,N,patch_index,patches,vortex_ring);
toc;

end

lattice_index = length(vortex_lattice);
N_vortices = 0; 

for k = 1:lattice_index

    [m,n] = size(vortex_lattice{1,k});
    N_vortices = N_vortices + m*n;

end

%% Vortices vector creation

vortex_rings = cell(1,N_vortices);
[n_rows,null] = size(vortex_lattice{1,1});
n_cols = N_vortices/m;
ratio = n_rows/n_cols;


for t = 1 : n_rows

col_buff = 0;    
    
 for i = 1 : lattice_index
    
  vortex_lattice_scan = vortex_lattice{1,i};
  [n_rows_i,col_i] = size(vortex_lattice_scan);

   for j = 1 : col_i
      vortex_rings{1,j+(col_buff)+n_cols*(t-1)} = vortex_lattice_scan{t,j};
   end
  
  col_buff = col_buff + col_i;

 end
end

vortices_mat = cell(n_rows,n_cols);

AR_panels_max = 0;
dist_min = 1000;

for i = 1:n_rows
    for j = 1:n_cols
    vortices_mat(i,j) = vortex_rings(j+n_cols*(i-1));
    
    b_buff = abs(vortices_mat{i,j}.Y_2-vortices_mat{i,j}.Y_1);
    c_buff = abs(vortices_mat{i,j}.X_3-vortices_mat{i,j}.X_1);
    
    AR_panels_buff = b_buff/c_buff;
    
    if b_buff < dist_min || c_buff < dist_min
       dist_min = min([b_buff,c_buff]);
    else
    end
    
    if AR_panels_buff  > AR_panels_max
       AR_panels_max = AR_panels_buff;
    end
        
    end
end

%% To create spanwise mesh and chord distribution

eta = zeros(1,n_cols);
s = zeros(1,n_cols+1);
chord_vec = zeros(1,n_cols+1);
staz = zeros(1,n_cols+1);

for j = 1:n_cols
    
    eta(j) = vortices_mat{1,j}.Y_C/(b_ref/2);
    s(j+1) = j/n_cols*b_tot/2;
    chord_vec(j) = abs(vortices_mat{n_rows,j}.X_3-vortices_mat{1,j}.X_1);
    staz(j) = vortices_mat{1,j}.Y_1;
end

s(1) = 0;
chord_vec(n_cols+1) = abs(vortices_mat{n_rows,n_cols}.X_4-vortices_mat{1,n_cols}.X_2);
staz(n_cols+1) = vortices_mat{1,n_cols}.Y_2;

chord = chord_vec;

%%

for i = 1:n_rows
    for j = 1:n_cols
        XX(i,j) = (vortices_mat{i,j}.X_C);
        YY(i,j) = (vortices_mat{i,j}.Y_C);
        ZZ(i,j) = (vortices_mat{i,j}.Z_C);
        
        NX(i,j) = (vortices_mat{i,j}.N_x);
        NY(i,j) = (vortices_mat{i,j}.N_y);
        NZ(i,j) = (vortices_mat{i,j}.N_z);
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

patch([X_1;X_2;X_4;X_3],[Y_1;Y_2;Y_4;Y_3],[Z_1;Z_2;Z_4;Z_3],[1,1,1],'Marker','o','markersize',3);hold on;
patch([X_1;X_2;X_4;X_3],[-Y_1;-Y_2;-Y_4;-Y_3],[Z_1;Z_2;Z_4;Z_3],[1,1,1],'Marker','o','markersize',3);hold on;
quiver3(XX,YY,ZZ,NX,NY,NZ,'k','Marker','x','markersize',3);hold on;
quiver3(XX,-YY,ZZ,NX,-NY,NZ,'k','Marker','o','markersize',3);hold on;

for index = 1:length(patches)
plot3([patches{1,index}.X_L,patches{1,index}.X_R,patches{1,index}.X_R+patches{1,index}.c_R,...
    patches{1,index}.X_L+patches{1,index}.c_L,patches{1,index}.X_L],...
      [patches{1,index}.Y_L,patches{1,index}.Y_R,patches{1,index}.Y_R,...
    patches{1,index}.Y_L,patches{1,index}.Y_L],...
    [patches{1,index}.Z_L,patches{1,index}.Z_R,patches{1,index}.Z_R,...
    patches{1,index}.Z_L,patches{1,index}.Z_L],...
    'r-');hold on
plot3([patches{1,index}.X_L,patches{1,index}.X_R,patches{1,index}.X_R+patches{1,index}.c_R,...
    patches{1,index}.X_L+patches{1,index}.c_L,patches{1,index}.X_L],...
      -[patches{1,index}.Y_L,patches{1,index}.Y_R,patches{1,index}.Y_R,...
    patches{1,index}.Y_L,patches{1,index}.Y_L],...
    [patches{1,index}.Z_L,patches{1,index}.Z_R,patches{1,index}.Z_R,...
    patches{1,index}.Z_L,patches{1,index}.Z_L],...
    'r-');hold on
end
xlabel('X');ylabel('Y');zlabel('Z');
axis equal;
grid on;
view(3);




%%

save([pwd,'/data/data_vortex_lattice.mat'],'vortex_rings','vortices_mat',...
    'n_rows','n_cols',...
    'eta','chord','s','dist_min');