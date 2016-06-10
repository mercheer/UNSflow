%% READ INPUT DATA FOR GIVEN WING GEOMETRY
% AND FORM PATCHES

close all;clc;clear all;

%% Adding paths

addpath([pwd,'/geometry']);
addpath([pwd,'/functions']);
addpath([pwd,'/data']);

%% Load geometry data

filename = [pwd,'/geometry/wing_geometry.dat'];
[colheaders,data] =func_import_wing_geometry(filename);


%% Create new variables in the base workspace from import.
for i = 1:size(colheaders, 2)
    assignin('base' , genvarname(colheaders{i}),data(:,i));
end

N_PATCHES = length(data(:,1));

patches = cell(1,N_PATCHES);

p_ix = 1;
S_ref = 0;
c_ref = 0;
b_ref_buff = zeros(1,N_PATCHES+1);
b_tot = 0;

 while p_ix <= N_PATCHES


     
     
%% Patches info


% CORNER POINTS

ix = 1;
X_LE_L = 0;
Y_LE_L = 0;
Z_LE_L = 0;


while ix < p_ix

    X_LE_L = X_LE_L + SPAN_LENGTH(ix)*sin(deg2rad(SWEEP_ANGLE(ix)));
    Y_LE_L = Y_LE_L + SPAN_LENGTH(ix)*cos(deg2rad(DIHEDRAL_ANGLE(ix)));
    Z_LE_L = Z_LE_L + SPAN_LENGTH(ix)*sin(deg2rad(DIHEDRAL_ANGLE(ix)));
    
    ix = ix + 1;
end

ix = 1;
X_LE_R = 0;
Y_LE_R = 0;
Z_LE_R = 0;

while ix <= p_ix

    X_LE_R = X_LE_R + SPAN_LENGTH(ix)*sin(deg2rad(SWEEP_ANGLE(ix)));
    Y_LE_R = Y_LE_R + SPAN_LENGTH(ix)*cos(deg2rad(DIHEDRAL_ANGLE(ix)));
    Z_LE_R = Z_LE_R + SPAN_LENGTH(ix)*sin(deg2rad(DIHEDRAL_ANGLE(ix)));
    
    ix = ix + 1;
end

% Collect patch infos in a struct(CORNER POINTS,CHORD AND TWIST)

patch_p_ix = struct('X_L',X_LE_L,'Y_L',Y_LE_L,'Z_L',Z_LE_L,...
    'X_R',X_LE_R,'Y_R',Y_LE_R,'Z_R',Z_LE_R,'c_L',CHORD_IN(p_ix),'c_R',CHORD_OUT(p_ix),...
    'twist_L',deg2rad(TWIST_ANGLE_IN(p_ix)),'twist_R',deg2rad(TWIST_ANGLE_OUT(p_ix)));

patches{1,p_ix} = patch_p_ix;

S_ref = S_ref + (CHORD_OUT(p_ix)+CHORD_IN(p_ix))/2 * SPAN_LENGTH(p_ix);
b_ref_buff(p_ix+1) = b_ref_buff(p_ix) + SPAN_LENGTH(p_ix)*cos(deg2rad(SWEEP_ANGLE(p_ix)))*cos(deg2rad(DIHEDRAL_ANGLE(p_ix))); 
c_ref = c_ref + (CHORD_OUT(p_ix)+CHORD_IN(p_ix))/2;
b_tot = b_tot + abs(SPAN_LENGTH(p_ix)*cos(deg2rad(SWEEP_ANGLE(p_ix)))*cos(deg2rad(DIHEDRAL_ANGLE(p_ix))));

p_ix = p_ix + 1;


 end


%% Reference parameters

S_ref = S_ref*2;
b_ref = max(b_ref_buff)*2;
c_ref = c_ref/N_PATCHES;
b_tot = b_tot*2;


%% Plot the wing planform (correct planform check)

% X-Y PLANFORM

figure()

for i = 1 : N_PATCHES
   
    
plot([patches{1,i}.X_L,patches{1,i}.X_R,...
    patches{1,i}.X_R+patches{1,i}.c_R,patches{1,i}.X_L+patches{1,i}.c_L,patches{1,i}.X_L]...
    ,[patches{1,i}.Y_L,patches{1,i}.Y_R,patches{1,i}.Y_R,patches{1,i}.Y_L,patches{1,i}.Y_L]...
    ,'-k');hold on;
axis equal;

plot([patches{1,i}.X_L,patches{1,i}.X_R,...
    patches{1,i}.X_R+patches{1,i}.c_R,patches{1,i}.X_L+patches{1,i}.c_L,patches{1,i}.X_L]...
    ,-[patches{1,i}.Y_L,patches{1,i}.Y_R,patches{1,i}.Y_R,patches{1,i}.Y_L,patches{1,i}.Y_L]...
    ,'-k');hold on;
axis equal;

xlabel('X'),ylabel('Y');

end

figure()

for i = 1 : N_PATCHES
   
    
plot3([patches{1,i}.X_L,patches{1,i}.X_R,...
    patches{1,i}.X_R+patches{1,i}.c_R,patches{1,i}.X_L+patches{1,i}.c_L,patches{1,i}.X_L]...
    ,[patches{1,i}.Y_L,patches{1,i}.Y_R,patches{1,i}.Y_R,patches{1,i}.Y_L,patches{1,i}.Y_L]...
    ,[patches{1,i}.Z_L,patches{1,i}.Z_R,patches{1,i}.Z_R,patches{1,i}.Z_L,patches{1,i}.Z_L],'-k');hold on;
axis equal;

plot3([patches{1,i}.X_L,patches{1,i}.X_R,...
    patches{1,i}.X_R+patches{1,i}.c_R,patches{1,i}.X_L+patches{1,i}.c_L,patches{1,i}.X_L]...
    ,-[patches{1,i}.Y_L,patches{1,i}.Y_R,patches{1,i}.Y_R,patches{1,i}.Y_L,patches{1,i}.Y_L]...
    ,[patches{1,i}.Z_L,patches{1,i}.Z_R,patches{1,i}.Z_R,patches{1,i}.Z_L,patches{1,i}.Z_L],'-k');hold on;
axis equal;

xlabel('X'),ylabel('Y');zlabel('Z');

end



Y_kink = 0;

%% Save Patches into a .mat file

save([pwd,'/data/data_patches.mat'],'patches');
save([pwd,'/data/data_refs.mat'],'S_ref','b_ref','b_tot','c_ref');
