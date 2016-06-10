function [C_L,C_Di,C_Di_st,C_Di_unst,C_L_add_mat,C_Di_add_mat,F_mat,L_mat,spanload,adim_spanload] = func_forces_KATZPLOTKIN(dt,vortices_mat,gamma_mat,gamma_mat_p,Q_inf,U_wake,V_wake,W_wake,W_ind_TV,U,V,W,u_mot_CP,v_mot_CP,w_mot_CP)

global n_rows n_cols eta s chord


%% Adding paths

addpath([pwd,'/functions']);
addpath([pwd,'/data']);



%% To load geometry and gamma distribution

load([pwd,'/data/data_refs.mat']);

abs_Q_inf = sqrt(Q_inf(1)^2+Q_inf(2)^2+Q_inf(3)^2);
vers_Q_inf = Q_inf./abs_Q_inf;

for i = 1:n_rows
    for j = 1:n_cols
vortex_rings(1,j+n_cols*(i-1)) = vortices_mat(i,j);
    end
end


%% To compute Total Lift and Induced Drag Coefficients of the wing 

% First creating empty matrices of C_L and C_Di for each vorticity segment

delta_b = zeros(n_rows,n_cols);
delta_c = zeros(n_rows,n_cols);
tau_i   = cell(n_rows,n_cols);
tau_j   = cell(n_rows,n_cols);
V1_wake  = cell(n_rows,n_cols);
V_mot   = cell(n_rows,n_cols);
V_tot   = cell(n_rows,n_cols);
abs_V_mot = zeros(n_rows,n_cols);
S       = zeros(n_rows,n_cols);
N_x     = zeros(n_rows,n_cols);
N_y     = zeros(n_rows,n_cols);
N_z     = zeros(n_rows,n_cols);
N_v     = cell(n_rows,n_cols);
P       = cell(n_rows,n_cols);

C_Di_add_mat = zeros(n_rows,n_cols);
C_Di_st_add_mat = zeros(n_rows,n_cols);
C_Di_unst_add_mat = zeros(n_rows,n_cols);
C_L_add_mat = zeros(n_rows,n_cols);
L_mat = zeros(n_rows,n_cols);
C_F_mat = cell(n_rows,n_cols);
F_mat   = cell(n_rows,n_cols);
L_mat = zeros(n_rows,n_cols);


temp_CL = 0;
temp_CD = 0;
temp_CD_st = 0;
temp_CD_unst = 0;


% For each span-wise bound-vorticity segment Kutta-Joukovski theorem is
% applied considering the free-stream velocity acting on CPs

for i = 1:n_rows
    
    for j = 1:n_cols
           delta_b(i,j) = sqrt((vortices_mat{i,j}.X_2-vortices_mat{i,j}.X_1)^2+...
               (vortices_mat{i,j}.Y_2-vortices_mat{i,j}.Y_1)^2+...
               (vortices_mat{i,j}.Z_2-vortices_mat{i,j}.Z_1)^2);
           delta_c(i,j) = sqrt((vortices_mat{i,j}.X_3-vortices_mat{i,j}.X_1)^2+...
               (vortices_mat{i,j}.Y_3-vortices_mat{i,j}.Y_1)^2+...
               (vortices_mat{i,j}.Z_3-vortices_mat{i,j}.Z_1)^2);
           
           tau_i{i,j} = [(vortices_mat{i,j}.X_3-vortices_mat{i,j}.X_1);...
               (vortices_mat{i,j}.Y_3-vortices_mat{i,j}.Y_1);...
               (vortices_mat{i,j}.Z_3-vortices_mat{i,j}.Z_1)]./delta_c(i,j);
           tau_j{i,j} = [(vortices_mat{i,j}.X_2-vortices_mat{i,j}.X_1);...
               (vortices_mat{i,j}.Y_2-vortices_mat{i,j}.Y_1);...
               (vortices_mat{i,j}.Z_2-vortices_mat{i,j}.Z_1)]./delta_b(i,j);
           
           V1_wake{i,j}   = [U_wake(i,j);V_wake(i,j);W_wake(i,j)+W_ind_TV(i,j)];
           V_mot{i,j} = [u_mot_CP(i,j);v_mot_CP(i,j);w_mot_CP(i,j)];
           V_tot{i,j} = [U_wake(i,j)+u_mot_CP(i,j);V_wake(i,j)+v_mot_CP(i,j);W_wake(i,j)+w_mot_CP(i,j)];
           
           abs_V_mot(i,j) = norm(V_mot{i,j});
           norm_V_mot = V_mot{i,j}./abs_V_mot(i,j);
           
           S(i,j)    = vortices_mat{i,j}.S;
           N_x(i,j)       = vortices_mat{i,j}.N_x;
           N_y(i,j)       = vortices_mat{i,j}.N_y;
           N_z(i,j)       = vortices_mat{i,j}.N_z;
           P{i,j}   = eye(3,3) - norm_V_mot*transpose(norm_V_mot);
           
           N_v{i,j} = [N_x(i,j);N_y(i,j);N_z(i,j)];
           
        if i == 1
            if j == 1
               
            alfa = atan2(dot(V_mot{i,j},N_v{i,j}),dot(V_mot{i,j},tau_i{i,j}));    
            C_L_add_mat(i,j) = (dot(V_tot{i,j},tau_i{i,j})*gamma_mat(i,j)/delta_c(i,j)+... 
                         (gamma_mat(i,j)-gamma_mat_p(i,j))/dt)*S(i,j)*cos(alfa)/...
                      (1/2*S_ref*abs_Q_inf^2);
            C_Di_add_mat(i,j) = 2*(-(dot((V1_wake{i,j}),P{i,j}*N_v{i,j})*gamma_mat(i,j))*...
                 delta_b(i,j) + (gamma_mat(i,j)-gamma_mat_p(i,j))*S(i,j)*...
                 sin(alfa)/dt)/(S_ref*abs_Q_inf^2);
            C_Di_st_add_mat(i,j) = 2*(-(dot((V1_wake{i,j}),P{i,j}*N_v{i,j})*gamma_mat(i,j))*...
                 delta_b(i,j))/(S_ref*abs_Q_inf^2);
            C_Di_unst_add_mat(i,j) = 2*((gamma_mat(i,j)-gamma_mat_p(i,j))*S(i,j)*...
                 sin(alfa)/dt)/(S_ref*abs_Q_inf^2);
                                    
            else
            alfa = atan2(dot(V_mot{i,j},N_v{i,j}),dot(V_mot{i,j},tau_i{i,j}));    
            C_L_add_mat(i,j) = (dot(V_tot{i,j},tau_i{i,j})*gamma_mat(i,j)/delta_c(i,j)+... 
                         dot(V_tot{i,j},tau_j{i,j})*(gamma_mat(i,j)-gamma_mat(i,j-1))/delta_b(i,j)+... 
                         (gamma_mat(i,j)-gamma_mat_p(i,j))/dt)*S(i,j)*cos(alfa)/...
                      (1/2*S_ref*abs_Q_inf^2);
            C_Di_add_mat(i,j) = 2*(-(dot((V1_wake{i,j}),P{i,j}*N_v{i,j})*gamma_mat(i,j))*...
                 delta_b(i,j) + (gamma_mat(i,j)-gamma_mat_p(i,j))*S(i,j)*...
                 sin(alfa)/dt)/(S_ref*abs_Q_inf^2);
            C_Di_st_add_mat(i,j) = 2*(-(dot((V1_wake{i,j}),P{i,j}*N_v{i,j})*gamma_mat(i,j))*...
                 delta_b(i,j))/(S_ref*abs_Q_inf^2);
            C_Di_unst_add_mat(i,j) = 2*((gamma_mat(i,j)-gamma_mat_p(i,j))*S(i,j)*...
                 sin(alfa)/dt)/(S_ref*abs_Q_inf^2);          
            end    
                

        else
           if j == 1
               
            alfa = atan2(dot(V_mot{i,j},N_v{i,j}),dot(V_mot{i,j},tau_i{i,j}));    
            C_L_add_mat(i,j) = (dot(V_tot{i,j},tau_i{i,j})*(gamma_mat(i,j)-gamma_mat(i-1,j))/delta_c(i,j)+... 
                         (gamma_mat(i,j)-gamma_mat_p(i,j))/dt)*S(i,j)*cos(alfa)/...
                      (1/2*S_ref*abs_Q_inf^2);
            C_Di_add_mat(i,j) = 2*(-dot((V1_wake{i,j}),P{i,j}*N_v{i,j})*(gamma_mat(i,j)-gamma_mat(i-1,j))*...
                 delta_b(i,j) + (gamma_mat(i,j)-gamma_mat_p(i,j))*S(i,j)*...
                 sin(alfa)/dt)/(S_ref*abs_Q_inf^2);
            C_Di_st_add_mat(i,j) = 2*(-dot((V1_wake{i,j}),P{i,j}*N_v{i,j})*(gamma_mat(i,j)-gamma_mat(i-1,j))*...
                 delta_b(i,j))/(S_ref*abs_Q_inf^2);
            C_Di_unst_add_mat(i,j) = 2*((gamma_mat(i,j)-gamma_mat_p(i,j))*S(i,j)*...
                 sin(alfa)/dt)/(S_ref*abs_Q_inf^2);             
            
           else
            alfa = atan2(dot(V_mot{i,j},N_v{i,j}),dot(V_mot{i,j},tau_i{i,j}));    
            C_L_add_mat(i,j) = (dot(V_tot{i,j},tau_i{i,j})*(gamma_mat(i,j)-gamma_mat(i-1,j))/delta_c(i,j)+... 
                         dot(V_tot{i,j},tau_j{i,j})*(gamma_mat(i,j)-gamma_mat(i,j-1))/delta_b(i,j)+...  
                         (gamma_mat(i,j)-gamma_mat_p(i,j))/dt)*S(i,j)*cos(alfa)/...
                      (1/2*S_ref*abs_Q_inf^2);
            C_Di_add_mat(i,j) = 2*(-dot((V1_wake{i,j}),P{i,j}*N_v{i,j})*(gamma_mat(i,j)-gamma_mat(i-1,j))*...
                 delta_b(i,j) + (gamma_mat(i,j)-gamma_mat_p(i,j))*S(i,j)*...
                 sin(alfa)/dt)/(S_ref*abs_Q_inf^2);
            C_Di_st_add_mat(i,j) = 2*(-dot((V1_wake{i,j}),P{i,j}*N_v{i,j})*(gamma_mat(i,j)-gamma_mat(i-1,j))*...
                 delta_b(i,j))/(S_ref*abs_Q_inf^2);
            C_Di_unst_add_mat(i,j) = 2*((gamma_mat(i,j)-gamma_mat_p(i,j))*S(i,j)*...
                 sin(alfa)/dt)/(S_ref*abs_Q_inf^2); 
          
           end
                   
        end
        
        C_F_mat{i,j} = C_Di_add_mat(i,j).*norm_V_mot + C_L_add_mat(i,j).* (P{i,j}*N_v{i,j});
        F_mat{i,j}   = (C_Di_add_mat(i,j).*norm_V_mot + C_L_add_mat(i,j).* (P{i,j}*N_v{i,j})).*(1/2*S_ref*abs_Q_inf^2);
        L_mat(i,j)   = C_L_add_mat(i,j).* (1/2*S_ref*abs_Q_inf^2);
       
        
        temp_CL = temp_CL + C_L_add_mat(i,j);  % Adding terms to the total C_L coefficient
        temp_CD = temp_CD + C_Di_add_mat(i,j); % Adding terms to the total C_Di coefficient
        temp_CD_st = temp_CD_st + C_Di_st_add_mat(i,j);
        temp_CD_unst = temp_CD_unst + C_Di_unst_add_mat(i,j);

    end
    

end

% To take into account of left-half of the wing

C_Di = 2*temp_CD;
C_Di_st = 2*temp_CD_st;
C_Di_unst = 2*temp_CD_unst;
C_L = 2*temp_CL;



gamma_ad = zeros(1,n_cols);
spanload = zeros(1,n_cols);
chord_s  = zeros(1,n_cols);

for j = 1:n_cols
    
    chord_s(j) = interp1(s,chord,j,'linear','extrap');
    
    
    for i = 1:n_rows
        
        if i == 1
            
        gamma_ad(j) = gamma_ad(j)+gamma_mat(i,j)/(b_ref*abs_Q_inf);
        
        else
        
        gamma_ad(j) = gamma_ad(j)+(gamma_mat(i,j)-gamma_mat(i-1,j))/(b_ref*abs_Q_inf);
    
            
        end
        
        spanload(j) = gamma_ad(j)*2*b_ref;
    end
    

end

adim_spanload = spanload./chord_s;
KATZ_VERIFY = adim_spanload/C_L;


 %% To Plot wing with forces

% XX = zeros(n_rows,n_cols);
% YY = zeros(n_rows,n_cols);
% ZZ = zeros(n_rows,n_cols);
% CFX = zeros(n_rows,n_cols);
% CFY = zeros(n_rows,n_cols);
% CFZ = zeros(n_rows,n_cols);
% 
% for i = 1:n_rows
%     for j = 1:n_cols
%         XX(i,j) = (vortices_mat{i,j}.X_C);
%         YY(i,j) = (vortices_mat{i,j}.Y_C);
%         ZZ(i,j) = (vortices_mat{i,j}.Z_C);
%         CFX(i,j)  = C_F_mat{i,j}(1);
%         CFY(i,j)  = C_F_mat{i,j}(2);
%         CFZ(i,j)  = C_F_mat{i,j}(3);
%     end
% end
% 
% for i = 1:length(vortex_rings)
%     X_1(i) = vortex_rings{i}.X_1;
%     X_2(i) = vortex_rings{i}.X_2;
%     X_3(i) = vortex_rings{i}.X_3;
%     X_4(i) = vortex_rings{i}.X_4;
%     Y_1(i) = vortex_rings{i}.Y_1;
%     Y_2(i) = vortex_rings{i}.Y_2;
%     Y_3(i) = vortex_rings{i}.Y_3;
%     Y_4(i) = vortex_rings{i}.Y_4;
%     Z_1(i) = vortex_rings{i}.Z_1;
%     Z_2(i) = vortex_rings{i}.Z_2;
%     Z_3(i) = vortex_rings{i}.Z_3;
%     Z_4(i) = vortex_rings{i}.Z_4;
% end
% 
% C_L_vector = zeros(1,n_rows*n_cols);
% 
% for i = 1:n_rows
% 
%     for j = 1:n_cols
%     
%  C_L_vector(j+n_cols*(i-1))=C_L_add_mat(i,j);
% 
%     end
%     
% end
% 
% drawnow
% figure()
% title(['C_L = ' num2str(2*temp_CL) ]);
% patch([X_1;X_2;X_4;X_3],[Y_1;Y_2;Y_4;Y_3],[Z_1;Z_2;Z_4;Z_3],C_L_vector);hold on;
% patch([X_1;X_2;X_4;X_3],[-Y_1;-Y_2;-Y_4;-Y_3],[Z_1;Z_2;Z_4;Z_3],C_L_vector);hold on;
% quiver3(XX,YY,ZZ,CFX,CFY,CFZ,'k');hold on;
% quiver3(XX,-YY,ZZ,CFX,-CFY,CFZ,'k');hold on;
% 
% xlabel('X');ylabel('Y');zlabel('Z');
% axis equal;
% h = colorbar;
% set(get(h,'title'),'string','C_{L}');
% colormap('hot'); 


% quiver3(XX,YY,ZZ,U,V,W,'k');hold on;
% quiver3(XX,-YY,ZZ,U,-V,W,'k');hold on;

% figure
% plot([-fliplr(s(1:end-1)),s(1:end-1)],[fliplr(spanload),spanload],'.-k');
% xlabel('\eta');ylabel('c(\eta) C_{l}(\eta)');
% set(gca,'Xlim',[-s(end)-0.05 s(end)+0.05]);
% % set(gca,'Ylim',[0 max(spanload)+max(spanload)*0.05]);
% 
% figure
% plot([-fliplr(s(1:end-1)),s(1:end-1)],[fliplr(adim_spanload),adim_spanload],'.-k');
% xlabel('\eta');ylabel('C_{l}(\eta)');
% set(gca,'Xlim',[-s(end)-0.05 s(end)+0.05]);
% % set(gca,'Ylim',[0 max(adim_spanload)+max(adim_spanload)*0.05]);


end