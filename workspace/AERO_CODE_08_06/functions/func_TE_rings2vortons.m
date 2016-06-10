function func_TE_rings2vortons(now_index,vortex_rings)
%rings2vortons This function converts NEAR WAKE TRAILING EDGE vortex rings into FAR WAKE
%vortons
%   Detailed explanation goes here

global  n_cols FW_wake_vortons


vorton = struct('X',[],'Y',[],'Z',[],...
                'A_X',[],'A_Y',[],'A_Z',[]...
                );


%Maany ways to do the conversion of vortex rings to vortons

%% Vortex rings to vortons conversion


for j = 1:n_cols

    FW_wake_vortons{now_index-2,j} = vorton;

    FW_wake_vortons{now_index-2,j}.X = (vortex_rings{now_index-2,j}.X_4+vortex_rings{now_index-2,j}.X_3)/2;
    FW_wake_vortons{now_index-2,j}.Y = (vortex_rings{now_index-2,j}.Y_4+vortex_rings{now_index-2,j}.Y_3)/2;
    FW_wake_vortons{now_index-2,j}.Z = (vortex_rings{now_index-2,j}.Z_4+vortex_rings{now_index-2,j}.Z_3)/2;

    t1_kj = [vortex_rings{now_index-2,j}.X_2-vortex_rings{now_index-2,j}.X_1;
          vortex_rings{now_index-2,j}.Y_2-vortex_rings{now_index-2,j}.Y_1;
          vortex_rings{now_index-2,j}.Z_2-vortex_rings{now_index-2,j}.Z_1];
    t2_kj = [vortex_rings{now_index-2,j}.X_4-vortex_rings{now_index-2,j}.X_2;
          vortex_rings{now_index-2,j}.Y_4-vortex_rings{now_index-2,j}.Y_2;
          vortex_rings{now_index-2,j}.Z_4-vortex_rings{now_index-2,j}.Z_2];
    t3_kj = [vortex_rings{now_index-2,j}.X_3-vortex_rings{now_index-2,j}.X_4;
          vortex_rings{now_index-2,j}.Y_3-vortex_rings{now_index-2,j}.Y_4;
          vortex_rings{now_index-2,j}.Z_3-vortex_rings{now_index-2,j}.Z_4];
    t4_kj = [vortex_rings{now_index-2,j}.X_1-vortex_rings{now_index-2,j}.X_3;
          vortex_rings{now_index-2,j}.Y_1-vortex_rings{now_index-2,j}.Y_3;
          vortex_rings{now_index-2,j}.Z_1-vortex_rings{now_index-2,j}.Z_3];



   if now_index == 3

    if j == 1

      A_vect = + (vortex_rings{now_index-2,j}.GAMMA).* t3_kj +...
                   0.5*(vortex_rings{now_index-2,j}.GAMMA - vortex_rings{now_index-2,j+1}.GAMMA).*t2_kj;
    else
      if j == n_cols

      A_vect = 0.5*(vortex_rings{now_index-2,j}.GAMMA - vortex_rings{now_index-2,j-1}.GAMMA).*t4_kj +(vortex_rings{now_index-2,j}.GAMMA).* t3_kj +...
                   0.5*(vortex_rings{now_index-2,j}.GAMMA).*t2_kj;
      else
      A_vect = 0.5*(vortex_rings{now_index-2,j}.GAMMA - vortex_rings{now_index-2,j-1}.GAMMA).*t4_kj +(vortex_rings{now_index-2,j}.GAMMA).* t3_kj +...
                   0.5*(vortex_rings{now_index-2,j}.GAMMA - vortex_rings{now_index-2,j+1}.GAMMA).*t2_kj;
      end
    end
   else
     if j == 1

      A_vect =   + (vortex_rings{now_index-2,j}.GAMMA- vortex_rings{now_index-3,j}.GAMMA).* t3_kj +...
                   0.5*(vortex_rings{now_index-2,j}.GAMMA - vortex_rings{now_index-2,j+1}.GAMMA).*t2_kj;
     else
       if j == n_cols

      A_vect = 0.5*(vortex_rings{now_index-2,j}.GAMMA - vortex_rings{now_index-2,j-1}.GAMMA).*t4_kj +(vortex_rings{now_index-2,j}.GAMMA - vortex_rings{now_index-3,j}.GAMMA).* t3_kj +...
                   0.5*(vortex_rings{now_index-2,j}.GAMMA).*t2_kj;
       else
      A_vect = 0.5*(vortex_rings{now_index-2,j}.GAMMA - vortex_rings{now_index-2,j-1}.GAMMA).*t4_kj +(vortex_rings{now_index-2,j}.GAMMA - vortex_rings{now_index-3,j}.GAMMA).* t3_kj +...
                   0.5*(vortex_rings{now_index-2,j}.GAMMA - vortex_rings{now_index-2,j+1}.GAMMA).*t2_kj;
       end
     end

   end

 %  0.5*(vortex_rings{now_index-2,j}.GAMMA).*t4_kj


   FW_wake_vortons{now_index-2,j}.A_X = A_vect(1);
   FW_wake_vortons{now_index-2,j}.A_Y = A_vect(2);
   FW_wake_vortons{now_index-2,j}.A_Z = A_vect(3);


end


end
