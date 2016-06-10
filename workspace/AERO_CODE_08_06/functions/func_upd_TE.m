function vortices_mat_NOW = func_upd_TE(K,dt,Q_inf,vortices_mat_TE_NOW)

global n_rows n_cols 


vortices_mat_NOW = vortices_mat_TE_NOW;


%% To update TE segment position in Inertial reference system
i_last = n_rows;

    for j = 1:n_cols
        vortices_mat_NOW{i_last,j}.X_3 = vortices_mat_TE_NOW{i_last,j}.X_3 + K* Q_inf(1)*dt;
        vortices_mat_NOW{i_last,j}.X_4 = vortices_mat_TE_NOW{i_last,j}.X_4 + K* Q_inf(1)*dt;
        vortices_mat_NOW{i_last,j}.X_C = vortices_mat_TE_NOW{i_last,j}.X_C ;
        vortices_mat_NOW{i_last,j}.Y_3 = vortices_mat_TE_NOW{i_last,j}.Y_3 + K* Q_inf(2)*dt;
        vortices_mat_NOW{i_last,j}.Y_4 = vortices_mat_TE_NOW{i_last,j}.Y_4 + K* Q_inf(2)*dt;
        vortices_mat_NOW{i_last,j}.Y_C = vortices_mat_TE_NOW{i_last,j}.Y_C ;
        vortices_mat_NOW{i_last,j}.Z_3 = vortices_mat_TE_NOW{i_last,j}.Z_3 + K* Q_inf(3)*dt;
        vortices_mat_NOW{i_last,j}.Z_4 = vortices_mat_TE_NOW{i_last,j}.Z_4 + K* Q_inf(3)*dt;
        vortices_mat_NOW{i_last,j}.Z_C = vortices_mat_TE_NOW{i_last,j}.Z_C ;
        % Update panel area
        
        C1 = vortices_mat_NOW{i_last,j}.X_3 - vortices_mat_NOW{i_last,j}.X_1;
        C2 = vortices_mat_NOW{i_last,j}.Y_3 - vortices_mat_NOW{i_last,j}.Y_1;
        C3 = vortices_mat_NOW{i_last,j}.Z_3 - vortices_mat_NOW{i_last,j}.Z_1;
        D1 = vortices_mat_NOW{i_last,j}.X_2 - vortices_mat_NOW{i_last,j}.X_1;
        D2 = vortices_mat_NOW{i_last,j}.Y_2 - vortices_mat_NOW{i_last,j}.Y_1;
        D3 = vortices_mat_NOW{i_last,j}.Z_2 - vortices_mat_NOW{i_last,j}.Z_1;
        E1 = C2*D3 - C3*D2;
        E2 = D1*C3 - C1*D3;
        E3 = C1*D2 - C2*D1;
        
        vortices_mat_NOW{i_last,j}.S = sqrt(E1^2+E2^2+E3^2);
                                           
                                           
    end
