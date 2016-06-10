function [C_L,C_Di,C_Di_st,C_Di_unst,D_p,spanload,adim_spanload] = func_forces_plot(dt,gamma_mat,gamma_mat_p,Q_inf,U_ind,V_ind,W_ind,W_ind_TV,U,V,W)

global vortices_mat n_rows n_cols eta s chord


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
V_v     = cell(n_rows,n_cols);
S       = zeros(n_rows,n_cols);
N_x     = zeros(n_rows,n_cols);
N_y     = zeros(n_rows,n_cols);
N_z     = zeros(n_rows,n_cols);
N_v     = cell(n_rows,n_cols);

C_Di_add_mat = zeros(n_rows,n_cols);
C_Di_st_add_mat = zeros(n_rows,n_cols);
C_Di_unst_add_mat = zeros(n_rows,n_cols);
C_L_add_mat = zeros(n_rows,n_cols);
D_p = zeros(n_rows,n_cols);
C_F = cell(n_rows,n_cols);
abs_C_F = zeros(n_rows,n_cols);
vers_C_F = cell(n_rows,n_cols);

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
           V_v{i,j}  = [U_ind(i,j);V_ind(i,j);W_ind(i,j)];
           S(i,j)    = vortices_mat{i,j}.S;
           N_x(i,j)       = vortices_mat{i,j}.N_x;
           N_y(i,j)       = vortices_mat{i,j}.N_y;
           N_z(i,j)       = vortices_mat{i,j}.N_z;
           
           
           N_v{i,j} = [N_x(i,j);N_y(i,j);N_z(i,j)];
           
        if i == 1
            if j == 1
            D_p(i,j) = dot(V_v{i,j},tau_i{i,j})*gamma_mat(i,j)/delta_c(i,j)+... 
                         (gamma_mat(i,j)-gamma_mat_p(i,j))/dt;
                     
            C_F{i,j}  = [( D_p(i,j)* S(i,j) * N_v{i,j}(1));...
                          ( D_p(i,j)* S(i,j) * N_v{i,j}(2));...
                          ( D_p(i,j)* S(i,j) * N_v{i,j}(3))]./...
                      (1/2*S_ref*abs_Q_inf^2);
               
                              
            abs_C_F(i,j) = sqrt(C_F{i,j}(1)^2+C_F{i,j}(2)^2+C_F{i,j}(3)^2);
            vers_C_F{i,j} = C_F{i,j}./abs_C_F(i,j);
            
            alfa = atan2(W_ind(i,j),U_ind(i,j));
            ang = acos(dot(vers_C_F{i,j},vers_Q_inf));          
           
            if D_p(i,j)  >= 0 
            C_L_add_mat(i,j)  = abs_C_F(i,j)*sin(ang);
            else
            C_L_add_mat(i,j)  = - abs_C_F(i,j)*sin(ang); 
            end
            C_Di_add_mat(i,j) = 2*(((W_ind_TV(i,j)+V_v{i,j}(3)*cos(alfa)-V_v{i,j}(3)*sin(alfa))*gamma_mat(i,j))*...
                 delta_b(i,j) + (gamma_mat(i,j)-gamma_mat_p(i,j))*S(i,j)*...
                 sin(alfa)/dt)/(S_ref*abs_Q_inf^2);
            C_Di_st_add_mat(i,j) = 2*(((W_ind_TV(i,j)+V_v{i,j}(3)*cos(alfa)-V_v{i,j}(3)*sin(alfa))*gamma_mat(i,j))*...
                 delta_b(i,j))/(S_ref*abs_Q_inf^2);
            C_Di_unst_add_mat(i,j) = 2*((gamma_mat(i,j)-gamma_mat_p(i,j))*S(i,j)*...
                 sin(alfa)/dt)/(S_ref*abs_Q_inf^2);          
            
            else
                
            D_p(i,j) = dot(V_v{i,j},tau_i{i,j})*gamma_mat(i,j)/delta_c(i,j)+... 
                         dot(V_v{i,j},tau_j{i,j})*(gamma_mat(i,j)-gamma_mat(i,j-1))/delta_b(i,j)+...
                         (gamma_mat(i,j)-gamma_mat_p(i,j))/dt;
            C_F{i,j}  = [( D_p(i,j)* S(i,j) * N_v{i,j}(1));...
                          ( D_p(i,j)* S(i,j) * N_v{i,j}(2));...
                          ( D_p(i,j)* S(i,j) * N_v{i,j}(3))]./...
                      (1/2*S_ref*abs_Q_inf^2);
                 
            abs_C_F(i,j) = sqrt(C_F{i,j}(1)^2+C_F{i,j}(2)^2+C_F{i,j}(3)^2);
            vers_C_F{i,j} = C_F{i,j}./abs_C_F(i,j);
            
            alfa = atan2(W_ind(i,j),U_ind(i,j));
            ang = acos(dot(vers_C_F{i,j},vers_Q_inf));          
                      
            if D_p(i,j)  >= 0 
            C_L_add_mat(i,j)  = abs_C_F(i,j)*sin(ang);
            else
            C_L_add_mat(i,j)  = - abs_C_F(i,j)*sin(ang); 
            end
            C_Di_add_mat(i,j) = 2*(((W_ind_TV(i,j)+V_v{i,j}(3)*cos(alfa)-V_v{i,j}(3)*sin(alfa))*gamma_mat(i,j))*...
                 delta_b(i,j) + (gamma_mat(i,j)-gamma_mat_p(i,j))*S(i,j)*...
                 sin(alfa)/dt)/(S_ref*abs_Q_inf^2);
            C_Di_st_add_mat(i,j) = 2*(((W_ind_TV(i,j)+V_v{i,j}(3)*cos(alfa)-V_v{i,j}(3)*sin(alfa))*gamma_mat(i,j))*...
                 delta_b(i,j))/(S_ref*abs_Q_inf^2);
            C_Di_unst_add_mat(i,j) = 2*((gamma_mat(i,j)-gamma_mat_p(i,j))*S(i,j)*...
                 sin(alfa)/dt)/(S_ref*abs_Q_inf^2);       
            end
        else
           if j == 1
            D_p(i,j) = dot(V_v{i,j},tau_i{i,j})*(gamma_mat(i,j)-gamma_mat(i-1,j))/delta_c(i,j)+... 
                         (gamma_mat(i,j)-gamma_mat_p(i,j))/dt;
            C_F{i,j}  = [( D_p(i,j)* S(i,j) * N_v{i,j}(1));...
                          ( D_p(i,j)* S(i,j) * N_v{i,j}(2));...
                          ( D_p(i,j)* S(i,j) * N_v{i,j}(3))]./...
                      (1/2*S_ref*abs_Q_inf^2);
            
            abs_C_F(i,j) = sqrt(C_F{i,j}(1)^2+C_F{i,j}(2)^2+C_F{i,j}(3)^2);
            vers_C_F{i,j} = C_F{i,j}./abs_C_F(i,j);
            
            alfa = atan2(W_ind(i,j),U_ind(i,j));
            ang = acos(dot(vers_C_F{i,j},vers_Q_inf));          
                      
            if D_p(i,j)  >= 0 
            C_L_add_mat(i,j)  = abs_C_F(i,j)*sin(ang);
            else
            C_L_add_mat(i,j)  = - abs_C_F(i,j)*sin(ang); 
            end
            C_Di_add_mat(i,j) = 2*(((W_ind_TV(i,j)+V_v{i,j}(3)*cos(alfa)-V_v{i,j}(3)*sin(alfa))*(gamma_mat(i,j)-gamma_mat(i-1,j)))*...
                 delta_b(i,j) + (gamma_mat(i,j)-gamma_mat_p(i,j))*S(i,j)*...
                 sin(alfa)/dt)/(S_ref*abs_Q_inf^2);         
            C_Di_st_add_mat(i,j) = 2*(((W_ind_TV(i,j)+V_v{i,j}(3)*cos(alfa)-V_v{i,j}(3)*sin(alfa))*(gamma_mat(i,j)-gamma_mat(i-1,j)))*...
                 delta_b(i,j))/(S_ref*abs_Q_inf^2);
            C_Di_unst_add_mat(i,j) = 2*((gamma_mat(i,j)-gamma_mat_p(i,j))*S(i,j)*...
                 sin(alfa)/dt)/(S_ref*abs_Q_inf^2);                
            else
                
            D_p(i,j) = dot(V_v{i,j},tau_i{i,j})*(gamma_mat(i,j)-gamma_mat(i-1,j))/delta_c(i,j)+... 
                         dot(V_v{i,j},tau_j{i,j})*(gamma_mat(i,j)-gamma_mat(i,j-1))/delta_b(i,j)+...
                         (gamma_mat(i,j)-gamma_mat_p(i,j))/dt;
            C_F{i,j}  = [( D_p(i,j)* S(i,j) * N_v{i,j}(1));...
                          ( D_p(i,j)* S(i,j) * N_v{i,j}(2));...
                          ( D_p(i,j)* S(i,j) * N_v{i,j}(3))]./...
                      (1/2*S_ref*abs_Q_inf^2);
            

            abs_C_F(i,j) = sqrt(C_F{i,j}(1)^2+C_F{i,j}(2)^2+C_F{i,j}(3)^2);
            vers_C_F{i,j} = C_F{i,j}./abs_C_F(i,j);
            
            alfa = atan2(W_ind(i,j),U_ind(i,j));
            ang = acos(dot(vers_C_F{i,j},vers_Q_inf));          
                      
            if D_p(i,j)  >= 0 
            C_L_add_mat(i,j)  = abs_C_F(i,j)*sin(ang);
            else
            C_L_add_mat(i,j)  = - abs_C_F(i,j)*sin(ang); 
            end
            C_Di_add_mat(i,j) = 2*(((W_ind_TV(i,j)+V_v{i,j}(3)*cos(alfa)-V_v{i,j}(3)*sin(alfa))*(gamma_mat(i,j)-gamma_mat(i-1,j)))*...
                 delta_b(i,j) + (gamma_mat(i,j)-gamma_mat_p(i,j))*S(i,j)*...
                 sin(alfa)/dt)/(S_ref*abs_Q_inf^2);
            C_Di_st_add_mat(i,j) = 2*(((W_ind_TV(i,j)+V_v{i,j}(3)*cos(alfa)-V_v{i,j}(3)*sin(alfa))*(gamma_mat(i,j)-gamma_mat(i-1,j)))*...
                 delta_b(i,j))/(S_ref*abs_Q_inf^2);
            C_Di_unst_add_mat(i,j) = 2*((gamma_mat(i,j)-gamma_mat_p(i,j))*S(i,j)*...
                 sin(alfa)/dt)/(S_ref*abs_Q_inf^2);  
           end
                   
        end
        
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

XX = zeros(n_rows,n_cols);
YY = zeros(n_rows,n_cols);
ZZ = zeros(n_rows,n_cols);


for i = 1:n_rows
    for j = 1:n_cols
        XX(i,j) = (vortices_mat{i,j}.X_C);
        YY(i,j) = (vortices_mat{i,j}.Y_C);
        ZZ(i,j) = (vortices_mat{i,j}.Z_C);
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

C_p_vector = zeros(1,n_rows*n_cols);

for i = 1:n_rows

    for j = 1:n_cols
    
 C_p_vector(j+n_cols*(i-1))=abs_C_F(i,j);

    end
    
end

% figure()
% title(['C_L = ' num2str(2*temp_CL) ]);
% quiver3(XX,YY,ZZ,U,V,W,'k');hold on;
% q1 = quiver3(XX,YY,ZZ,U_ind,V_ind,W_ind,'k');hold on;
% patch([X_1;X_2;X_4;X_3],[Y_1;Y_2;Y_4;Y_3],[Z_1;Z_2;Z_4;Z_3],C_p_vector);
% q2 = quiver3(XX,-YY,ZZ,U_ind,-V_ind,W_ind,'k');hold on;
% patch([X_1;X_2;X_4;X_3],[-Y_1;-Y_2;-Y_4;-Y_3],[Z_1;Z_2;Z_4;Z_3],C_p_vector);
% set(q1,'Color',[.8 .8 .8]);
% set(q2,'Color',[.8 .8 .8]);
% axis equal;
% colorbar;
% colormap('hot'); 
% 
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