function func_wakerollup_VORTONS_TELE_PARALL(dt,now_index,vortices_CELL,gamma_mat,wake_vortex_ring_now,LE_wake_vortex_ring_now,FW_wake_vortons_now,LE_FW_wake_vortons_now,sig)

global NW_wake_vortex_rings  LE_NW_wake_vortex_rings FW_wake_vortons LE_FW_wake_vortons n_cols n_rows

%% ADDPATH

addpath([pwd,'/geometry']);
addpath([pwd,'/functions']);
addpath([pwd,'/data']);

%% LOAD

load([pwd,'/data/data_vortex_lattice.mat']);

%% LOOPS VORTONS


i = now_index-1;

while i > now_index-3
     for j = 1:n_cols
        X_1 = NW_wake_vortex_rings{i,j}.X_1;
        X_2 = NW_wake_vortex_rings{i,j}.X_2;
        X_3 = NW_wake_vortex_rings{i,j}.X_3;
        X_4 = NW_wake_vortex_rings{i,j}.X_4;
        Y_1 = NW_wake_vortex_rings{i,j}.Y_1;
        Y_2 = NW_wake_vortex_rings{i,j}.Y_2;
        Y_3 = NW_wake_vortex_rings{i,j}.Y_3;
        Y_4 = NW_wake_vortex_rings{i,j}.Y_4;
        Z_1 = NW_wake_vortex_rings{i,j}.Z_1;
        Z_2 = NW_wake_vortex_rings{i,j}.Z_2;
        Z_3 = NW_wake_vortex_rings{i,j}.Z_3;
        Z_4 = NW_wake_vortex_rings{i,j}.Z_4;

        % Loop for inducted velocity of bound vortex rings
u1 = 0;
u2 = 0;
u3 = 0;
u4 = 0;
v1 = 0;
v2 = 0;
v3 = 0;
v4 = 0;
w1 = 0;
w2 = 0;
w3 = 0;
w4 = 0;
u1w = 0;
u2w = 0;
u3w = 0;
u4w = 0;
v1w = 0;
v2w = 0;
v3w = 0;
v4w = 0;
w1w = 0;
w2w = 0;
w3w = 0;
w4w = 0;
u1w_v = 0;
u2w_v = 0;
u3w_v = 0;
u4w_v = 0;
v1w_v = 0;
v2w_v = 0;
v3w_v = 0;
v4w_v = 0;
w1w_v = 0;
w2w_v = 0;
w3w_v = 0;
w4w_v = 0;
LE_u1w = 0;
LE_u2w = 0;
LE_u3w = 0;
LE_u4w = 0;
LE_v1w = 0;
LE_v2w = 0;
LE_v3w = 0;
LE_v4w = 0;
LE_w1w = 0;
LE_w2w = 0;
LE_w3w = 0;
LE_w4w = 0;
LE_u1w_v = 0;
LE_u2w_v = 0;
LE_u3w_v = 0;
LE_u4w_v = 0;
LE_v1w_v = 0;
LE_v2w_v = 0;
LE_v3w_v = 0;
LE_v4w_v = 0;
LE_w1w_v = 0;
LE_w2w_v = 0;
LE_w3w_v = 0;
LE_w4w_v = 0;
u1_w_NWFW_tot = 0;
v1_w_NWFW_tot = 0;
w1_w_NWFW_tot = 0;
u2_w_NWFW_tot = 0;
v2_w_NWFW_tot = 0;
w2_w_NWFW_tot = 0;
u3_w_NWFW_tot = 0;
v3_w_NWFW_tot = 0;
w3_w_NWFW_tot = 0;
u4_w_NWFW_tot = 0;
v4_w_NWFW_tot = 0;
w4_w_NWFW_tot = 0;
LE_u1_w_NWFW_tot = 0;
LE_v1_w_NWFW_tot = 0;
LE_w1_w_NWFW_tot = 0;
LE_u2_w_NWFW_tot = 0;
LE_v2_w_NWFW_tot = 0;
LE_w2_w_NWFW_tot = 0;
LE_u3_w_NWFW_tot = 0;
LE_v3_w_NWFW_tot = 0;
LE_w3_w_NWFW_tot = 0;
LE_u4_w_NWFW_tot = 0;
LE_v4_w_NWFW_tot = 0;
LE_w4_w_NWFW_tot = 0;

%% Loop for wing induction


        for l = 1:n_rows
           for m = 1:n_cols

               [u1_ind_add_R,v1_ind_add_R,w1_ind_add_R] = func_voring_comp(X_1,Y_1,Z_1,vortices_CELL{l,m},...
                   gamma_mat(l,m),0);
               [u2_ind_add_R,v2_ind_add_R,w2_ind_add_R] = func_voring_comp(X_2,Y_2,Z_2,vortices_CELL{l,m},...
                   gamma_mat(l,m),0);
               [u3_ind_add_R,v3_ind_add_R,w3_ind_add_R] = func_voring_comp(X_3,Y_3,Z_3,vortices_CELL{l,m},...
                   gamma_mat(l,m),0);
               [u4_ind_add_R,v4_ind_add_R,w4_ind_add_R] = func_voring_comp(X_4,Y_4,Z_4,vortices_CELL{l,m},...
                   gamma_mat(l,m),0);
               [u1_ind_add_L,v1_ind_add_L,w1_ind_add_L] = func_voring_comp(X_1,-Y_1,Z_1,vortices_CELL{l,m},...
                   gamma_mat(l,m),0);
               [u2_ind_add_L,v2_ind_add_L,w2_ind_add_L] = func_voring_comp(X_2,-Y_2,Z_2,vortices_CELL{l,m},...
                   gamma_mat(l,m),0);
               [u3_ind_add_L,v3_ind_add_L,w3_ind_add_L] = func_voring_comp(X_3,-Y_3,Z_3,vortices_CELL{l,m},...
                   gamma_mat(l,m),0);
               [u4_ind_add_L,v4_ind_add_L,w4_ind_add_L] = func_voring_comp(X_4,-Y_4,Z_4,vortices_CELL{l,m},...
                   gamma_mat(l,m),0);

            u1 = u1 + u1_ind_add_R+ u1_ind_add_L;
            u2 = u2 + u2_ind_add_R+ u2_ind_add_L;
            u3 = u3 + u3_ind_add_R+ u3_ind_add_L;
            u4 = u4 + u4_ind_add_R+ u4_ind_add_L;
            v1 = v1 + v1_ind_add_R- v1_ind_add_L;
            v2 = v2 + v2_ind_add_R- v2_ind_add_L;
            v3 = v3 + v3_ind_add_R- v3_ind_add_L;
            v4 = v4 + v4_ind_add_R- v4_ind_add_L;
            w1 = w1 + w1_ind_add_R+ w1_ind_add_L;
            w2 = w2 + w2_ind_add_R+ w2_ind_add_L;
            w3 = w3 + w3_ind_add_R+ w3_ind_add_L;
            w4 = w4 + w4_ind_add_R+ w4_ind_add_L;
           end
        end

%% Loop for NEAR WAKE wake vortex rings

        r = now_index-1;

        while r > now_index-3
           for s = 1:n_cols

               [u1w_ind_add_R,v1w_ind_add_R,w1w_ind_add_R] = func_voring_comp(X_1,Y_1,Z_1,wake_vortex_ring_now{r,s},...
                   wake_vortex_ring_now{r,s}.GAMMA,0);
               [u2w_ind_add_R,v2w_ind_add_R,w2w_ind_add_R] = func_voring_comp(X_2,Y_2,Z_2,wake_vortex_ring_now{r,s},...
                   wake_vortex_ring_now{r,s}.GAMMA,0);
               [u3w_ind_add_R,v3w_ind_add_R,w3w_ind_add_R] = func_voring_comp(X_3,Y_3,Z_3,wake_vortex_ring_now{r,s},...
                   wake_vortex_ring_now{r,s}.GAMMA,0);
               [u4w_ind_add_R,v4w_ind_add_R,w4w_ind_add_R] = func_voring_comp(X_4,Y_4,Z_4,wake_vortex_ring_now{r,s},...
                   wake_vortex_ring_now{r,s}.GAMMA,0);
               [u1w_ind_add_L,v1w_ind_add_L,w1w_ind_add_L] = func_voring_comp(X_1,-Y_1,Z_1,wake_vortex_ring_now{r,s},...
                   wake_vortex_ring_now{r,s}.GAMMA,0);
               [u2w_ind_add_L,v2w_ind_add_L,w2w_ind_add_L] = func_voring_comp(X_2,-Y_2,Z_2,wake_vortex_ring_now{r,s},...
                   wake_vortex_ring_now{r,s}.GAMMA,0);
               [u3w_ind_add_L,v3w_ind_add_L,w3w_ind_add_L] = func_voring_comp(X_3,-Y_3,Z_3,wake_vortex_ring_now{r,s},...
                   wake_vortex_ring_now{r,s}.GAMMA,0);
               [u4w_ind_add_L,v4w_ind_add_L,w4w_ind_add_L] = func_voring_comp(X_4,-Y_4,Z_4,wake_vortex_ring_now{r,s},...
                   wake_vortex_ring_now{r,s}.GAMMA,0);
            u1w = u1w + u1w_ind_add_R+ u1w_ind_add_L;
            u2w = u2w + u2w_ind_add_R+ u2w_ind_add_L;
            u3w = u3w + u3w_ind_add_R+ u3w_ind_add_L;
            u4w = u4w + u4w_ind_add_R+ u4w_ind_add_L;
            v1w = v1w + v1w_ind_add_R- v1w_ind_add_L;
            v2w = v2w + v2w_ind_add_R- v2w_ind_add_L;
            v3w = v3w + v3w_ind_add_R- v3w_ind_add_L;
            v4w = v4w + v4w_ind_add_R- v4w_ind_add_L;
            w1w = w1w + w1w_ind_add_R+ w1w_ind_add_L;
            w2w = w2w + w2w_ind_add_R+ w2w_ind_add_L;
            w3w = w3w + w3w_ind_add_R+ w3w_ind_add_L;
            w4w = w4w + w4w_ind_add_R+ w4w_ind_add_L;
           end
           r = r - 1;
        end

%% Loop for NEAR WAKE LEADING EDGE wake vortex rings

        r = now_index-1;

        while r > now_index-3
           for s = 1:n_cols

               [LE_u1w_ind_add_R,LE_v1w_ind_add_R,LE_w1w_ind_add_R] = func_voring_comp(X_1,Y_1,Z_1,LE_wake_vortex_ring_now{r,s},...
                   LE_wake_vortex_ring_now{r,s}.GAMMA,0);
               [LE_u2w_ind_add_R,LE_v2w_ind_add_R,LE_w2w_ind_add_R] = func_voring_comp(X_2,Y_2,Z_2,LE_wake_vortex_ring_now{r,s},...
                   LE_wake_vortex_ring_now{r,s}.GAMMA,0);
               [LE_u3w_ind_add_R,LE_v3w_ind_add_R,LE_w3w_ind_add_R] = func_voring_comp(X_3,Y_3,Z_3,LE_wake_vortex_ring_now{r,s},...
                   LE_wake_vortex_ring_now{r,s}.GAMMA,0);
               [LE_u4w_ind_add_R,LE_v4w_ind_add_R,LE_w4w_ind_add_R] = func_voring_comp(X_4,Y_4,Z_4,LE_wake_vortex_ring_now{r,s},...
                   LE_wake_vortex_ring_now{r,s}.GAMMA,0);
               [LE_u1w_ind_add_L,LE_v1w_ind_add_L,LE_w1w_ind_add_L] = func_voring_comp(X_1,-Y_1,Z_1,LE_wake_vortex_ring_now{r,s},...
                   LE_wake_vortex_ring_now{r,s}.GAMMA,0);
               [LE_u2w_ind_add_L,LE_v2w_ind_add_L,LE_w2w_ind_add_L] = func_voring_comp(X_2,-Y_2,Z_2,LE_wake_vortex_ring_now{r,s},...
                   LE_wake_vortex_ring_now{r,s}.GAMMA,0);
               [LE_u3w_ind_add_L,LE_v3w_ind_add_L,LE_w3w_ind_add_L] = func_voring_comp(X_3,-Y_3,Z_3,LE_wake_vortex_ring_now{r,s},...
                   LE_wake_vortex_ring_now{r,s}.GAMMA,0);
               [LE_u4w_ind_add_L,LE_v4w_ind_add_L,LE_w4w_ind_add_L] = func_voring_comp(X_4,-Y_4,Z_4,LE_wake_vortex_ring_now{r,s},...
                   LE_wake_vortex_ring_now{r,s}.GAMMA,0);
            LE_u1w = LE_u1w + LE_u1w_ind_add_R+ LE_u1w_ind_add_L;
            LE_u2w = LE_u2w + LE_u2w_ind_add_R+ LE_u2w_ind_add_L;
            LE_u3w = LE_u3w + LE_u3w_ind_add_R+ LE_u3w_ind_add_L;
            LE_u4w = LE_u4w + LE_u4w_ind_add_R+ LE_u4w_ind_add_L;
            LE_v1w = LE_v1w + LE_v1w_ind_add_R- LE_v1w_ind_add_L;
            LE_v2w = LE_v2w + LE_v2w_ind_add_R- LE_v2w_ind_add_L;
            LE_v3w = LE_v3w + LE_v3w_ind_add_R- LE_v3w_ind_add_L;
            LE_v4w = LE_v4w + LE_v4w_ind_add_R- LE_v4w_ind_add_L;
            LE_w1w = LE_w1w + LE_w1w_ind_add_R+ LE_w1w_ind_add_L;
            LE_w2w = LE_w2w + LE_w2w_ind_add_R+ LE_w2w_ind_add_L;
            LE_w3w = LE_w3w + LE_w3w_ind_add_R+ LE_w3w_ind_add_L;
            LE_w4w = LE_w4w + LE_w4w_ind_add_R+ LE_w4w_ind_add_L;
           end
           r = r - 1;
        end


%% Loop for inductions of FAR WAKE wake vortons

        r = 1;

        while r < now_index-2
           if abs(r-i)<20
           for s = 1:n_cols

               [u1w_vind_add_R,v1w_vind_add_R,w1w_vind_add_R] = func_vortonind(X_1,Y_1,Z_1,FW_wake_vortons_now{r,s},sig);
               [u2w_vind_add_R,v2w_vind_add_R,w2w_vind_add_R] = func_vortonind(X_2,Y_2,Z_2,FW_wake_vortons_now{r,s},sig);
               [u3w_vind_add_R,v3w_vind_add_R,w3w_vind_add_R] = func_vortonind(X_3,Y_3,Z_3,FW_wake_vortons_now{r,s},sig);
               [u4w_vind_add_R,v4w_vind_add_R,w4w_vind_add_R] = func_vortonind(X_4,Y_4,Z_4,FW_wake_vortons_now{r,s},sig);
               [u1w_vind_add_L,v1w_vind_add_L,w1w_vind_add_L] = func_vortonind(X_1,-Y_1,Z_1,FW_wake_vortons_now{r,s},sig);
               [u2w_vind_add_L,v2w_vind_add_L,w2w_vind_add_L] = func_vortonind(X_2,-Y_2,Z_2,FW_wake_vortons_now{r,s},sig);
               [u3w_vind_add_L,v3w_vind_add_L,w3w_vind_add_L] = func_vortonind(X_3,-Y_3,Z_3,FW_wake_vortons_now{r,s},sig);
               [u4w_vind_add_L,v4w_vind_add_L,w4w_vind_add_L] = func_vortonind(X_4,-Y_4,Z_4,FW_wake_vortons_now{r,s},sig);

            u1w_v = u1w_v + u1w_vind_add_R+ u1w_vind_add_L;
            u2w_v = u2w_v + u2w_vind_add_R+ u2w_vind_add_L;
            u3w_v = u3w_v + u3w_vind_add_R+ u3w_vind_add_L;
            u4w_v = u4w_v + u4w_vind_add_R+ u4w_vind_add_L;
            v1w_v = v1w_v + v1w_vind_add_R- v1w_vind_add_L;
            v2w_v = v2w_v + v2w_vind_add_R- v2w_vind_add_L;
            v3w_v = v3w_v + v3w_vind_add_R- v3w_vind_add_L;
            v4w_v = v4w_v + v4w_vind_add_R- v4w_vind_add_L;
            w1w_v = w1w_v + w1w_vind_add_R+ w1w_vind_add_L;
            w2w_v = w2w_v + w2w_vind_add_R+ w2w_vind_add_L;
            w3w_v = w3w_v + w3w_vind_add_R+ w3w_vind_add_L;
            w4w_v = w4w_v + w4w_vind_add_R+ w4w_vind_add_L;
           end
           else
           end
           r = r + 1;
        end

%% Loop for inductions of FAR WAKE LEADING EDGE wake vortons

        r = 1;

        while r < now_index-2
           if abs(r-i)<20
           for s = 1:n_cols

               [LE_u1w_vind_add_R,LE_v1w_vind_add_R,LE_w1w_vind_add_R] = func_vortonind(X_1,Y_1,Z_1,LE_FW_wake_vortons_now{r,s},sig);
               [LE_u2w_vind_add_R,LE_v2w_vind_add_R,LE_w2w_vind_add_R] = func_vortonind(X_2,Y_2,Z_2,LE_FW_wake_vortons_now{r,s},sig);
               [LE_u3w_vind_add_R,LE_v3w_vind_add_R,LE_w3w_vind_add_R] = func_vortonind(X_3,Y_3,Z_3,LE_FW_wake_vortons_now{r,s},sig);
               [LE_u4w_vind_add_R,LE_v4w_vind_add_R,LE_w4w_vind_add_R] = func_vortonind(X_4,Y_4,Z_4,LE_FW_wake_vortons_now{r,s},sig);
               [LE_u1w_vind_add_L,LE_v1w_vind_add_L,LE_w1w_vind_add_L] = func_vortonind(X_1,-Y_1,Z_1,LE_FW_wake_vortons_now{r,s},sig);
               [LE_u2w_vind_add_L,LE_v2w_vind_add_L,LE_w2w_vind_add_L] = func_vortonind(X_2,-Y_2,Z_2,LE_FW_wake_vortons_now{r,s},sig);
               [LE_u3w_vind_add_L,LE_v3w_vind_add_L,LE_w3w_vind_add_L] = func_vortonind(X_3,-Y_3,Z_3,LE_FW_wake_vortons_now{r,s},sig);
               [LE_u4w_vind_add_L,LE_v4w_vind_add_L,LE_w4w_vind_add_L] = func_vortonind(X_4,-Y_4,Z_4,LE_FW_wake_vortons_now{r,s},sig);

            LE_u1w_v = LE_u1w_v + LE_u1w_vind_add_R+ LE_u1w_vind_add_L;
            LE_u2w_v = LE_u2w_v + LE_u2w_vind_add_R+ LE_u2w_vind_add_L;
            LE_u3w_v = LE_u3w_v + LE_u3w_vind_add_R+ LE_u3w_vind_add_L;
            LE_u4w_v = LE_u4w_v + LE_u4w_vind_add_R+ LE_u4w_vind_add_L;
            LE_v1w_v = LE_v1w_v + LE_v1w_vind_add_R- LE_v1w_vind_add_L;
            LE_v2w_v = LE_v2w_v + LE_v2w_vind_add_R- LE_v2w_vind_add_L;
            LE_v3w_v = LE_v3w_v + LE_v3w_vind_add_R- LE_v3w_vind_add_L;
            LE_v4w_v = LE_v4w_v + LE_v4w_vind_add_R- LE_v4w_vind_add_L;
            LE_w1w_v = LE_w1w_v + LE_w1w_vind_add_R+ LE_w1w_vind_add_L;
            LE_w2w_v = LE_w2w_v + LE_w2w_vind_add_R+ LE_w2w_vind_add_L;
            LE_w3w_v = LE_w3w_v + LE_w3w_vind_add_R+ LE_w3w_vind_add_L;
            LE_w4w_v = LE_w4w_v + LE_w4w_vind_add_R+ LE_w4w_vind_add_L;
           end
           else
           end
           r = r + 1;
        end


%% Loop for NEAR-WAKE - FAR-WAKE separation LINE


    for m = 1:n_cols
        X_1SEG = wake_vortex_ring_now{now_index-2,m}.X_3;
        X_2SEG = wake_vortex_ring_now{now_index-2,m}.X_4;
        Y_1SEG = wake_vortex_ring_now{now_index-2,m}.Y_3;
        Y_2SEG = wake_vortex_ring_now{now_index-2,m}.Y_4;
        Z_1SEG = wake_vortex_ring_now{now_index-2,m}.Z_3;
        Z_2SEG = wake_vortex_ring_now{now_index-2,m}.Z_4;
        GAMMASEG = wake_vortex_ring_now{now_index-3,m}.GAMMA;

         [u1_w_NWFW_add_R,v1_w_NWFW_add_R,w1_w_NWFW_add_R] = func_vortexl(X_1,Y_1,Z_1,X_1SEG,Y_1SEG,Z_1SEG,X_2SEG,Y_2SEG,Z_2SEG,GAMMASEG);
         [u1_w_NWFW_add_L,v1_w_NWFW_add_L,w1_w_NWFW_add_L] = func_vortexl(X_1,-Y_1,Z_1,X_1SEG,Y_1SEG,Z_1SEG,X_2SEG,Y_2SEG,Z_2SEG,GAMMASEG);
         u1_w_NWFW_tot = u1_w_NWFW_tot  + u1_w_NWFW_add_R +u1_w_NWFW_add_L;
         v1_w_NWFW_tot = v1_w_NWFW_tot + v1_w_NWFW_add_R -v1_w_NWFW_add_L;
         w1_w_NWFW_tot = w1_w_NWFW_tot + w1_w_NWFW_add_R +w1_w_NWFW_add_L;

         [u2_w_NWFW_add_R,v2_w_NWFW_add_R,w2_w_NWFW_add_R] = func_vortexl(X_2,Y_2,Z_2,X_1SEG,Y_1SEG,Z_1SEG,X_2SEG,Y_2SEG,Z_2SEG,GAMMASEG);
         [u2_w_NWFW_add_L,v2_w_NWFW_add_L,w2_w_NWFW_add_L] = func_vortexl(X_2,-Y_2,Z_2,X_1SEG,Y_1SEG,Z_1SEG,X_2SEG,Y_2SEG,Z_2SEG,GAMMASEG);
         u2_w_NWFW_tot = u2_w_NWFW_tot  + u2_w_NWFW_add_R +u2_w_NWFW_add_L;
         v2_w_NWFW_tot = v2_w_NWFW_tot + v2_w_NWFW_add_R -v2_w_NWFW_add_L;
         w2_w_NWFW_tot = w2_w_NWFW_tot + w2_w_NWFW_add_R +w2_w_NWFW_add_L;

         [u3_w_NWFW_add_R,v3_w_NWFW_add_R,w3_w_NWFW_add_R] = func_vortexl(X_3,Y_3,Z_3,X_1SEG,Y_1SEG,Z_1SEG,X_2SEG,Y_2SEG,Z_2SEG,GAMMASEG);
         [u3_w_NWFW_add_L,v3_w_NWFW_add_L,w3_w_NWFW_add_L] = func_vortexl(X_3,-Y_3,Z_3,X_1SEG,Y_1SEG,Z_1SEG,X_2SEG,Y_2SEG,Z_2SEG,GAMMASEG);
         u3_w_NWFW_tot = u3_w_NWFW_tot  + u3_w_NWFW_add_R +u3_w_NWFW_add_L;
         v3_w_NWFW_tot = v3_w_NWFW_tot + v3_w_NWFW_add_R -v3_w_NWFW_add_L;
         w3_w_NWFW_tot = w3_w_NWFW_tot + w3_w_NWFW_add_R +w3_w_NWFW_add_L;

         [u4_w_NWFW_add_R,v4_w_NWFW_add_R,w4_w_NWFW_add_R] = func_vortexl(X_4,Y_4,Z_4,X_1SEG,Y_1SEG,Z_1SEG,X_2SEG,Y_2SEG,Z_2SEG,GAMMASEG);
         [u4_w_NWFW_add_L,v4_w_NWFW_add_L,w4_w_NWFW_add_L] = func_vortexl(X_4,-Y_4,Z_4,X_1SEG,Y_1SEG,Z_1SEG,X_2SEG,Y_2SEG,Z_2SEG,GAMMASEG);
         u4_w_NWFW_tot = u4_w_NWFW_tot  + u4_w_NWFW_add_R +u4_w_NWFW_add_L;
         v4_w_NWFW_tot = v4_w_NWFW_tot + v4_w_NWFW_add_R -v4_w_NWFW_add_L;
         w4_w_NWFW_tot = w4_w_NWFW_tot + w4_w_NWFW_add_R +w4_w_NWFW_add_L;
    end

%% Loop for NEAR-WAKE - FAR-WAKE LEADING EDGE separation LINE


    for m = 1:n_cols
        X_1SEG = LE_wake_vortex_ring_now{now_index-2,m}.X_3;
        X_2SEG = LE_wake_vortex_ring_now{now_index-2,m}.X_4;
        Y_1SEG = LE_wake_vortex_ring_now{now_index-2,m}.Y_3;
        Y_2SEG = LE_wake_vortex_ring_now{now_index-2,m}.Y_4;
        Z_1SEG = LE_wake_vortex_ring_now{now_index-2,m}.Z_3;
        Z_2SEG = LE_wake_vortex_ring_now{now_index-2,m}.Z_4;
        GAMMASEG = LE_wake_vortex_ring_now{now_index-3,m}.GAMMA;

         [LE_u1_w_NWFW_add_R,LE_v1_w_NWFW_add_R,LE_w1_w_NWFW_add_R] = func_vortexl(X_1,Y_1,Z_1,X_1SEG,Y_1SEG,Z_1SEG,X_2SEG,Y_2SEG,Z_2SEG,GAMMASEG);
         [LE_u1_w_NWFW_add_L,LE_v1_w_NWFW_add_L,LE_w1_w_NWFW_add_L] = func_vortexl(X_1,-Y_1,Z_1,X_1SEG,Y_1SEG,Z_1SEG,X_2SEG,Y_2SEG,Z_2SEG,GAMMASEG);
		 [LE_u2_w_NWFW_add_R,LE_v2_w_NWFW_add_R,LE_w2_w_NWFW_add_R] = func_vortexl(X_2,Y_2,Z_2,X_1SEG,Y_1SEG,Z_1SEG,X_2SEG,Y_2SEG,Z_2SEG,GAMMASEG);
         [LE_u2_w_NWFW_add_L,LE_v2_w_NWFW_add_L,LE_w2_w_NWFW_add_L] = func_vortexl(X_2,-Y_2,Z_2,X_1SEG,Y_1SEG,Z_1SEG,X_2SEG,Y_2SEG,Z_2SEG,GAMMASEG);
		 [LE_u3_w_NWFW_add_R,LE_v3_w_NWFW_add_R,LE_w3_w_NWFW_add_R] = func_vortexl(X_3,Y_3,Z_3,X_1SEG,Y_1SEG,Z_1SEG,X_2SEG,Y_2SEG,Z_2SEG,GAMMASEG);
         [LE_u3_w_NWFW_add_L,LE_v3_w_NWFW_add_L,LE_w3_w_NWFW_add_L] = func_vortexl(X_3,-Y_3,Z_3,X_1SEG,Y_1SEG,Z_1SEG,X_2SEG,Y_2SEG,Z_2SEG,GAMMASEG);
		 [LE_u4_w_NWFW_add_R,LE_v4_w_NWFW_add_R,LE_w4_w_NWFW_add_R] = func_vortexl(X_4,Y_4,Z_4,X_1SEG,Y_1SEG,Z_1SEG,X_2SEG,Y_2SEG,Z_2SEG,GAMMASEG);
         [LE_u4_w_NWFW_add_L,LE_v4_w_NWFW_add_L,LE_w4_w_NWFW_add_L] = func_vortexl(X_4,-Y_4,Z_4,X_1SEG,Y_1SEG,Z_1SEG,X_2SEG,Y_2SEG,Z_2SEG,GAMMASEG);

         LE_u1_w_NWFW_tot =LE_u1_w_NWFW_tot + LE_u1_w_NWFW_add_R +LE_u1_w_NWFW_add_L;
         LE_v1_w_NWFW_tot =LE_v1_w_NWFW_tot + LE_v1_w_NWFW_add_R -LE_v1_w_NWFW_add_L;
         LE_w1_w_NWFW_tot =LE_w1_w_NWFW_tot + LE_w1_w_NWFW_add_R +LE_w1_w_NWFW_add_L;
         LE_u2_w_NWFW_tot =LE_u2_w_NWFW_tot + LE_u2_w_NWFW_add_R +LE_u2_w_NWFW_add_L;
         LE_v2_w_NWFW_tot =LE_v2_w_NWFW_tot + LE_v2_w_NWFW_add_R -LE_v2_w_NWFW_add_L;
         LE_w2_w_NWFW_tot =LE_w2_w_NWFW_tot + LE_w2_w_NWFW_add_R +LE_w2_w_NWFW_add_L;
         LE_u3_w_NWFW_tot =LE_u3_w_NWFW_tot + LE_u3_w_NWFW_add_R +LE_u3_w_NWFW_add_L;
         LE_v3_w_NWFW_tot =LE_v3_w_NWFW_tot + LE_v3_w_NWFW_add_R -LE_v3_w_NWFW_add_L;
         LE_w3_w_NWFW_tot =LE_w3_w_NWFW_tot + LE_w3_w_NWFW_add_R +LE_w3_w_NWFW_add_L;
         LE_u4_w_NWFW_tot =LE_u4_w_NWFW_tot + LE_u4_w_NWFW_add_R +LE_u4_w_NWFW_add_L;
         LE_v4_w_NWFW_tot =LE_v4_w_NWFW_tot + LE_v4_w_NWFW_add_R -LE_v4_w_NWFW_add_L;
         LE_w4_w_NWFW_tot =LE_w4_w_NWFW_tot + LE_w4_w_NWFW_add_R +LE_w4_w_NWFW_add_L;
    end


%%

        dX_1 = (u1+u1w+u1w_v+u1_w_NWFW_tot+LE_u1w+LE_u1w_v+LE_u1_w_NWFW_tot)*dt;
        dX_2 = (u2+u2w+u2w_v+u2_w_NWFW_tot+LE_u2w+LE_u2w_v+LE_u2_w_NWFW_tot)*dt;
        dX_3 = (u3+u3w+u3w_v+u3_w_NWFW_tot+LE_u3w+LE_u3w_v+LE_u3_w_NWFW_tot)*dt;
        dX_4 = (u4+u4w+u4w_v+u4_w_NWFW_tot+LE_u4w+LE_u4w_v+LE_u4_w_NWFW_tot)*dt;
        dY_1 = (v1+v1w+v1w_v+v1_w_NWFW_tot+LE_v1w+LE_v1w_v+LE_v1_w_NWFW_tot)*dt;
        dY_2 = (v2+v2w+v2w_v+v2_w_NWFW_tot+LE_v2w+LE_v2w_v+LE_v2_w_NWFW_tot)*dt;
        dY_3 = (v3+v3w+v3w_v+v3_w_NWFW_tot+LE_v3w+LE_v3w_v+LE_v3_w_NWFW_tot)*dt;
        dY_4 = (v4+v4w+v4w_v+v4_w_NWFW_tot+LE_v4w+LE_v4w_v+LE_v4_w_NWFW_tot)*dt;
        dZ_1 = (w1+w1w+w1w_v+w1_w_NWFW_tot+LE_w1w+LE_w1w_v+LE_w1_w_NWFW_tot)*dt;
        dZ_2 = (w2+w2w+w2w_v+w2_w_NWFW_tot+LE_w2w+LE_w2w_v+LE_w2_w_NWFW_tot)*dt;
        dZ_3 = (w3+w3w+w3w_v+w3_w_NWFW_tot+LE_w3w+LE_w3w_v+LE_w3_w_NWFW_tot)*dt;
        dZ_4 = (w4+w4w+w4w_v+w4_w_NWFW_tot+LE_w4w+LE_w4w_v+LE_w4_w_NWFW_tot)*dt;





        NW_wake_vortex_rings{i,j}.X_1 = wake_vortex_ring_now{i,j}.X_1 + dX_1 ;
        NW_wake_vortex_rings{i,j}.X_2 = wake_vortex_ring_now{i,j}.X_2 + dX_2 ;
        NW_wake_vortex_rings{i,j}.X_3 = wake_vortex_ring_now{i,j}.X_3 + dX_3 ;
        NW_wake_vortex_rings{i,j}.X_4 = wake_vortex_ring_now{i,j}.X_4 + dX_4 ;
        NW_wake_vortex_rings{i,j}.Y_1 = wake_vortex_ring_now{i,j}.Y_1 + dY_1 ;
        NW_wake_vortex_rings{i,j}.Y_2 = wake_vortex_ring_now{i,j}.Y_2 + dY_2 ;
        NW_wake_vortex_rings{i,j}.Y_3 = wake_vortex_ring_now{i,j}.Y_3 + dY_3 ;
        NW_wake_vortex_rings{i,j}.Y_4 = wake_vortex_ring_now{i,j}.Y_4 + dY_4 ;
        NW_wake_vortex_rings{i,j}.Z_1 = wake_vortex_ring_now{i,j}.Z_1 + dZ_1 ;
        NW_wake_vortex_rings{i,j}.Z_2 = wake_vortex_ring_now{i,j}.Z_2 + dZ_2 ;
        NW_wake_vortex_rings{i,j}.Z_3 = wake_vortex_ring_now{i,j}.Z_3 + dZ_3 ;
        NW_wake_vortex_rings{i,j}.Z_4 = wake_vortex_ring_now{i,j}.Z_4 + dZ_4 ;

     end


i = i - 1;
end







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cols = n_cols;
rows = n_rows;

FW_wake_vortons_buffer = cell(1,cols*(now_index-3));
FW_wake_vortons_now_buffer = cell(1,cols*(now_index-3));

k = 1;
i = 1;

while i < now_index-2
    for j = 1:cols

   FW_wake_vortons_buffer{k} = FW_wake_vortons{i,j};
   FW_wake_vortons_now_buffer{k}    = FW_wake_vortons_now{i,j};
   k = k + 1;
    end
   i = i + 1;
end



parfor cyc = 1:cols*(now_index-3)
        X_1 = FW_wake_vortons_buffer{cyc}.X;
        Y_1 = FW_wake_vortons_buffer{cyc}.Y;
        Z_1 = FW_wake_vortons_buffer{cyc}.Z;
        alfa_x = FW_wake_vortons_buffer{cyc}.A_X;
        alfa_y = FW_wake_vortons_buffer{cyc}.A_Y;
        alfa_z = FW_wake_vortons_buffer{cyc}.A_Z;

        alfa_vec = [alfa_x,alfa_y,alfa_z];

AX_w = 0;
AY_w = 0;
AZ_w = 0;
AX_wvr = 0;
AY_wvr = 0;
AZ_wvr = 0;
AX_v = 0;
AY_v = 0;
AZ_v = 0;
AX_w_NWFW_tot = 0;
AY_w_NWFW_tot = 0;
AZ_w_NWFW_tot = 0;
LE_AX_wvr = 0;
LE_AY_wvr = 0;
LE_AZ_wvr = 0;
LE_AX_v = 0;
LE_AY_v = 0;
LE_AZ_v = 0;
LE_AX_w_NWFW_tot = 0;
LE_AY_w_NWFW_tot = 0;
LE_AZ_w_NWFW_tot = 0;


u1 = 0;
v1 = 0;
w1 = 0;
u1w = 0;
v1w = 0;
w1w = 0;
u1w_v = 0;
v1w_v = 0;
w1w_v = 0;
u1_w_NWFW_tot = 0;
v1_w_NWFW_tot = 0;
w1_w_NWFW_tot = 0;
LE_u1w = 0;
LE_v1w = 0;
LE_w1w = 0;
LE_u1w_v = 0;
LE_v1w_v = 0;
LE_w1w_v = 0;
LE_u1_w_NWFW_tot = 0;
LE_v1_w_NWFW_tot = 0;
LE_w1_w_NWFW_tot = 0;

%% Loop for inducted velocity of bound vortex rings

        for l = 1:rows
           for m = 1:cols

               [u1_ind_add_R,v1_ind_add_R,w1_ind_add_R] = func_voring_comp(X_1,Y_1,Z_1,vortices_CELL{l,m},...
                   gamma_mat(l,m),0);
               [u1_ind_add_L,v1_ind_add_L,w1_ind_add_L] = func_voring_comp(X_1,-Y_1,Z_1,vortices_CELL{l,m},...
                   gamma_mat(l,m),0);


            u1 = u1 + u1_ind_add_R+ u1_ind_add_L;
            v1 = v1 + v1_ind_add_R- v1_ind_add_L;
            w1 = w1 + w1_ind_add_R+ w1_ind_add_L;

           end
        end

%% Loop for inducted velocity of leading edge wake vortex rings
        r = now_index-1;

        while r > now_index-3
           for s = 1:cols

               [LE_u1w_ind_add_R,LE_v1w_ind_add_R,LE_w1w_ind_add_R] = func_voring_comp(X_1,Y_1,Z_1,LE_wake_vortex_ring_now{r,s},...
                   wake_vortex_ring_now{r,s}.GAMMA,0);
               [LE_u1w_ind_add_L,LE_v1w_ind_add_L,LE_w1w_ind_add_L] = func_voring_comp(X_1,-Y_1,Z_1,LE_wake_vortex_ring_now{r,s},...
                   wake_vortex_ring_now{r,s}.GAMMA,0);

            LE_u1w = LE_u1w + LE_u1w_ind_add_R+ LE_u1w_ind_add_L;
            LE_v1w = LE_v1w + LE_v1w_ind_add_R- LE_v1w_ind_add_L;
            LE_w1w = LE_w1w + LE_w1w_ind_add_R+ LE_w1w_ind_add_L;

           end
           r = r - 1;
        end

%% Loop for inducted velocity of wake vortex rings
        r = now_index-1;

        while r > now_index-3
           for s = 1:cols

               [u1w_ind_add_R,v1w_ind_add_R,w1w_ind_add_R] = func_voring_comp(X_1,Y_1,Z_1,wake_vortex_ring_now{r,s},...
                   wake_vortex_ring_now{r,s}.GAMMA,0);
               [u1w_ind_add_L,v1w_ind_add_L,w1w_ind_add_L] = func_voring_comp(X_1,-Y_1,Z_1,wake_vortex_ring_now{r,s},...
                   wake_vortex_ring_now{r,s}.GAMMA,0);

            u1w = u1w + u1w_ind_add_R+ u1w_ind_add_L;
            v1w = v1w + v1w_ind_add_R- v1w_ind_add_L;
            w1w = w1w + w1w_ind_add_R+ w1w_ind_add_L;

           end
           r = r - 1;
        end


%% Loop for inducted velocity of wake vortons
        r = 1;

        while r < now_index-2
           for s = 1:cols


               [u1w_vind_add_R,v1w_vind_add_R,w1w_vind_add_R] = func_vortonind(X_1,Y_1,Z_1,FW_wake_vortons_now{r,s},sig);
               [u1w_vind_add_L,v1w_vind_add_L,w1w_vind_add_L] = func_vortonind(X_1,-Y_1,Z_1,FW_wake_vortons_now{r,s},sig);


            u1w_v = u1w_v + u1w_vind_add_R+ u1w_vind_add_L;
            v1w_v = v1w_v + v1w_vind_add_R- v1w_vind_add_L;
            w1w_v = w1w_v + w1w_vind_add_R+ w1w_vind_add_L;

           end
           r = r + 1;
        end

 %% Loop for inducted velocity of wake vortons
        r = 1;

        while r < now_index-2
           for s = 1:cols


               [LE_u1w_vind_add_R,LE_v1w_vind_add_R,LE_w1w_vind_add_R] = func_vortonind(X_1,Y_1,Z_1,LE_FW_wake_vortons_now{r,s},sig);
               [LE_u1w_vind_add_L,LE_v1w_vind_add_L,LE_w1w_vind_add_L] = func_vortonind(X_1,-Y_1,Z_1,LE_FW_wake_vortons_now{r,s},sig);


            LE_u1w_v = LE_u1w_v + LE_u1w_vind_add_R+ LE_u1w_vind_add_L;
            LE_v1w_v = LE_v1w_v + LE_v1w_vind_add_R- LE_v1w_vind_add_L;
            LE_w1w_v = LE_w1w_v + LE_w1w_vind_add_R+ LE_w1w_vind_add_L;

           end
           r = r + 1;
        end

%% Loop for NEAR-WAKE - FAR-WAKE separation line

    for m = 1:cols
        X_1SEG = wake_vortex_ring_now{now_index-2,m}.X_3;
        X_2SEG = wake_vortex_ring_now{now_index-2,m}.X_4;
        Y_1SEG = wake_vortex_ring_now{now_index-2,m}.Y_3;
        Y_2SEG = wake_vortex_ring_now{now_index-2,m}.Y_4;
        Z_1SEG = wake_vortex_ring_now{now_index-2,m}.Z_3;
        Z_2SEG = wake_vortex_ring_now{now_index-2,m}.Z_4;
        GAMMASEG = wake_vortex_ring_now{now_index-3,m}.GAMMA;

         [u1_w_NWFW_add_R,v1_w_NWFW_add_R,w1_w_NWFW_add_R] = func_vortexl(X_1,Y_1,Z_1,X_1SEG,Y_1SEG,Z_1SEG,X_2SEG,Y_2SEG,Z_2SEG,GAMMASEG);
         [u1_w_NWFW_add_L,v1_w_NWFW_add_L,w1_w_NWFW_add_L] = func_vortexl(X_1,-Y_1,Z_1,X_1SEG,Y_1SEG,Z_1SEG,X_2SEG,Y_2SEG,Z_2SEG,GAMMASEG);
         u1_w_NWFW_tot = u1_w_NWFW_tot + u1_w_NWFW_add_R +u1_w_NWFW_add_L;
         v1_w_NWFW_tot = v1_w_NWFW_tot + v1_w_NWFW_add_R -v1_w_NWFW_add_L;
         w1_w_NWFW_tot = w1_w_NWFW_tot + w1_w_NWFW_add_R +w1_w_NWFW_add_L;
    end

%% Loop for NEAR-WAKE - FAR-WAKE LEADING-EDGE separation line

    for m = 1:cols
        X_1SEG = LE_wake_vortex_ring_now{now_index-2,m}.X_3;
        X_2SEG = LE_wake_vortex_ring_now{now_index-2,m}.X_4;
        Y_1SEG = LE_wake_vortex_ring_now{now_index-2,m}.Y_3;
        Y_2SEG = LE_wake_vortex_ring_now{now_index-2,m}.Y_4;
        Z_1SEG = LE_wake_vortex_ring_now{now_index-2,m}.Z_3;
        Z_2SEG = LE_wake_vortex_ring_now{now_index-2,m}.Z_4;
        GAMMASEG = LE_wake_vortex_ring_now{now_index-3,m}.GAMMA;

        [LE_u1_w_NWFW_add_R,LE_v1_w_NWFW_add_R,LE_w1_w_NWFW_add_R] = func_vortexl(X_1,Y_1,Z_1,X_1SEG,Y_1SEG,Z_1SEG,X_2SEG,Y_2SEG,Z_2SEG,GAMMASEG);
        [LE_u1_w_NWFW_add_L,LE_v1_w_NWFW_add_L,LE_w1_w_NWFW_add_L] = func_vortexl(X_1,-Y_1,Z_1,X_1SEG,Y_1SEG,Z_1SEG,X_2SEG,Y_2SEG,Z_2SEG,GAMMASEG);
         LE_u1_w_NWFW_tot = LE_u1_w_NWFW_tot + LE_u1_w_NWFW_add_R +LE_u1_w_NWFW_add_L;
         LE_v1_w_NWFW_tot = LE_v1_w_NWFW_tot + LE_v1_w_NWFW_add_R -LE_v1_w_NWFW_add_L;
         LE_w1_w_NWFW_tot = LE_w1_w_NWFW_tot + LE_w1_w_NWFW_add_R +LE_w1_w_NWFW_add_L;
    end

%% Loop for inducted velocity of bound vortex rings


        for l = 1:rows
           for m = 1:cols

               vel_grad_mat_R = func_deriv_voring(X_1,Y_1,Z_1,vortices_CELL{l,m},...
                   gamma_mat(l,m));
               vel_grad_mat_L = func_deriv_voring(X_1,-Y_1,Z_1,vortices_CELL{l,m},...
                   gamma_mat(l,m));

            var_vect_R = alfa_vec * vel_grad_mat_R;
            var_vect_L = alfa_vec * vel_grad_mat_L;

            u1_ind_add_R = var_vect_R(1);
            v1_ind_add_R = var_vect_R(2);
            w1_ind_add_R = var_vect_R(3);

            u1_ind_add_L = var_vect_L(1);
            v1_ind_add_L = var_vect_L(2);
            w1_ind_add_L = var_vect_L(3);

            AX_w = AX_w + u1_ind_add_R+ u1_ind_add_L;
            AY_w = AY_w + v1_ind_add_R- v1_ind_add_L;
            AZ_w = AZ_w + w1_ind_add_R+ w1_ind_add_L;

           end
        end

%% Loop for inducted velocity of wake vortex rings

        r = now_index-1;

        while r > now_index-3
           for s = 1:cols

                vel_grad_mat_R = func_deriv_voring(X_1,Y_1,Z_1,wake_vortex_ring_now{r,s},...
                   wake_vortex_ring_now{r,s}.GAMMA);
                vel_grad_mat_L = func_deriv_voring(X_1,-Y_1,Z_1,wake_vortex_ring_now{r,s},...
                   wake_vortex_ring_now{r,s}.GAMMA);

            var_vect_R = alfa_vec * vel_grad_mat_R;
            var_vect_L = alfa_vec * vel_grad_mat_L;

            u1w_ind_add_R = var_vect_R(1);
            v1w_ind_add_R = var_vect_R(2);
            w1w_ind_add_R = var_vect_R(3);

            u1w_ind_add_L = var_vect_L(1);
            v1w_ind_add_L = var_vect_L(2);
            w1w_ind_add_L = var_vect_L(3);

            AX_wvr = AX_wvr + u1w_ind_add_R+ u1w_ind_add_L;
            AY_wvr = AY_wvr + v1w_ind_add_R- v1w_ind_add_L;
            AZ_wvr = AZ_wvr + w1w_ind_add_R+ w1w_ind_add_L;

           end
           r = r - 1;
        end
%% Loop for inducted velocity of LE wake vortex rings

        r = now_index-1;

        while r > now_index-3
           for s = 1:cols

                vel_grad_mat_R = func_deriv_voring(X_1,Y_1,Z_1,LE_wake_vortex_ring_now{r,s},...
                   LE_wake_vortex_ring_now{r,s}.GAMMA);
                vel_grad_mat_L = func_deriv_voring(X_1,-Y_1,Z_1,LE_wake_vortex_ring_now{r,s},...
                   LE_wake_vortex_ring_now{r,s}.GAMMA);

            var_vect_R = alfa_vec * vel_grad_mat_R;
            var_vect_L = alfa_vec * vel_grad_mat_L;

            u1w_ind_add_R = var_vect_R(1);
            v1w_ind_add_R = var_vect_R(2);
            w1w_ind_add_R = var_vect_R(3);

            u1w_ind_add_L = var_vect_L(1);
            v1w_ind_add_L = var_vect_L(2);
            w1w_ind_add_L = var_vect_L(3);

            LE_AX_wvr = LE_AX_wvr + u1w_ind_add_R+ u1w_ind_add_L;
            LE_AY_wvr = LE_AY_wvr + v1w_ind_add_R- v1w_ind_add_L;
            LE_AZ_wvr = LE_AZ_wvr + w1w_ind_add_R+ w1w_ind_add_L;

           end
           r = r - 1;
        end


%% Loop for inducted velocity of wake vortons

        r = 1;

        while r < now_index-2
           for s = 1:cols


               [u1w_vind_add_R,v1w_vind_add_R,w1w_vind_add_R] = func_vortonstre(FW_wake_vortons_buffer{cyc},FW_wake_vortons_now{r,s},sig);
               [u1w_vind_add_L,v1w_vind_add_L,w1w_vind_add_L] = func_vortonstre(FW_wake_vortons_buffer{cyc},FW_wake_vortons_now{r,s},sig);


            AX_v = AX_v + u1w_vind_add_R+ u1w_vind_add_L;
            AY_v = AY_v + v1w_vind_add_R- v1w_vind_add_L;
            AZ_v = AZ_v + w1w_vind_add_R+ w1w_vind_add_L;

           end
           r = r + 1;
        end

%% Loop for inducted velocity of wake vortons

        r = 1;

        while r < now_index-2
           for s = 1:cols


               [u1w_vind_add_R,v1w_vind_add_R,w1w_vind_add_R] = func_vortonstre(FW_wake_vortons_buffer{cyc},LE_FW_wake_vortons_now{r,s},sig);
               [u1w_vind_add_L,v1w_vind_add_L,w1w_vind_add_L] = func_vortonstre(FW_wake_vortons_buffer{cyc},LE_FW_wake_vortons_now{r,s},sig);


            LE_AX_v = LE_AX_v + u1w_vind_add_R+ u1w_vind_add_L;
            LE_AY_v = LE_AY_v + v1w_vind_add_R- v1w_vind_add_L;
            LE_AZ_v = LE_AZ_v + w1w_vind_add_R+ w1w_vind_add_L;

           end
           r = r + 1;
        end

%% Loop for NEAR-WAKE - FAR-WAKE separation line

    for m = 1:cols
        X_1SEG = wake_vortex_ring_now{now_index-2,m}.X_3;
        X_2SEG = wake_vortex_ring_now{now_index-2,m}.X_4;
        Y_1SEG = wake_vortex_ring_now{now_index-2,m}.Y_3;
        Y_2SEG = wake_vortex_ring_now{now_index-2,m}.Y_4;
        Z_1SEG = wake_vortex_ring_now{now_index-2,m}.Z_3;
        Z_2SEG = wake_vortex_ring_now{now_index-2,m}.Z_4;
        GAMMASEG = wake_vortex_ring_now{now_index-3,m}.GAMMA;

         [dudx_val_R,dudy_val_R,dudz_val_R,dvdx_val_R,dvdy_val_R,dvdz_val_R,dwdx_val_R,dwdy_val_R,dwdz_val_R] = func_deriv_vortexl(X_1,Y_1,Z_1,X_1SEG,Y_1SEG,Z_1SEG,X_2SEG,Y_2SEG,Z_2SEG,GAMMASEG);
         [dudx_val_L,dudy_val_L,dudz_val_L,dvdx_val_L,dvdy_val_L,dvdz_val_L,dwdx_val_L,dwdy_val_L,dwdz_val_L] = func_deriv_vortexl(X_1,-Y_1,Z_1,X_1SEG,Y_1SEG,Z_1SEG,X_2SEG,Y_2SEG,Z_2SEG,GAMMASEG);

            var_vect_R = alfa_vec * [dudx_val_R,dudy_val_R,dudz_val_R;...
                                    dvdx_val_R,dvdy_val_R,dvdz_val_R;...
                                    dwdx_val_R,dwdy_val_R,dwdz_val_R];

            var_vect_L = alfa_vec * [dudx_val_L,dudy_val_L,dudz_val_L;...
                                    dvdx_val_L,dvdy_val_L,dvdz_val_L;...
                                    dwdx_val_L,dwdy_val_L,dwdz_val_L];

            u1_w_NWFW_add_R = var_vect_R(1);
            v1_w_NWFW_add_R = var_vect_R(2);
            w1_w_NWFW_add_R = var_vect_R(3);

            u1_w_NWFW_add_L = var_vect_L(1);
            v1_w_NWFW_add_L = var_vect_L(2);
            w1_w_NWFW_add_L = var_vect_L(3);

         AX_w_NWFW_tot = AX_w_NWFW_tot + u1_w_NWFW_add_R +u1_w_NWFW_add_L;
         AY_w_NWFW_tot = AY_w_NWFW_tot + v1_w_NWFW_add_R -v1_w_NWFW_add_L;
         AZ_w_NWFW_tot = AZ_w_NWFW_tot + w1_w_NWFW_add_R +w1_w_NWFW_add_L;
    end

%% Loop for LE NEAR-WAKE - FAR-WAKE separation line

    for m = 1:cols
        X_1SEG = LE_wake_vortex_ring_now{now_index-2,m}.X_3;
        X_2SEG = LE_wake_vortex_ring_now{now_index-2,m}.X_4;
        Y_1SEG = LE_wake_vortex_ring_now{now_index-2,m}.Y_3;
        Y_2SEG = LE_wake_vortex_ring_now{now_index-2,m}.Y_4;
        Z_1SEG = LE_wake_vortex_ring_now{now_index-2,m}.Z_3;
        Z_2SEG = LE_wake_vortex_ring_now{now_index-2,m}.Z_4;
        GAMMASEG = LE_wake_vortex_ring_now{now_index-3,m}.GAMMA;

         [dudx_val_R,dudy_val_R,dudz_val_R,dvdx_val_R,dvdy_val_R,dvdz_val_R,dwdx_val_R,dwdy_val_R,dwdz_val_R] = func_deriv_vortexl(X_1,Y_1,Z_1,X_1SEG,Y_1SEG,Z_1SEG,X_2SEG,Y_2SEG,Z_2SEG,GAMMASEG);
         [dudx_val_L,dudy_val_L,dudz_val_L,dvdx_val_L,dvdy_val_L,dvdz_val_L,dwdx_val_L,dwdy_val_L,dwdz_val_L] = func_deriv_vortexl(X_1,-Y_1,Z_1,X_1SEG,Y_1SEG,Z_1SEG,X_2SEG,Y_2SEG,Z_2SEG,GAMMASEG);

            var_vect_R = alfa_vec * [dudx_val_R,dudy_val_R,dudz_val_R;...
                                    dvdx_val_R,dvdy_val_R,dvdz_val_R;...
                                    dwdx_val_R,dwdy_val_R,dwdz_val_R];

            var_vect_L = alfa_vec * [dudx_val_L,dudy_val_L,dudz_val_L;...
                                    dvdx_val_L,dvdy_val_L,dvdz_val_L;...
                                    dwdx_val_L,dwdy_val_L,dwdz_val_L];

            u1_w_NWFW_add_R = var_vect_R(1);
            v1_w_NWFW_add_R = var_vect_R(2);
            w1_w_NWFW_add_R = var_vect_R(3);

            u1_w_NWFW_add_L = var_vect_L(1);
            v1_w_NWFW_add_L = var_vect_L(2);
            w1_w_NWFW_add_L = var_vect_L(3);

         LE_AX_w_NWFW_tot = LE_AX_w_NWFW_tot + u1_w_NWFW_add_R +u1_w_NWFW_add_L;
         LE_AY_w_NWFW_tot = LE_AY_w_NWFW_tot + v1_w_NWFW_add_R -v1_w_NWFW_add_L;
         LE_AZ_w_NWFW_tot = LE_AZ_w_NWFW_tot + w1_w_NWFW_add_R +w1_w_NWFW_add_L;
    end
%% VORTONS UPDATING

        dA_X = (AX_w+AX_wvr+AX_v+AX_w_NWFW_tot+LE_AX_wvr+LE_AX_v+LE_AX_w_NWFW_tot)*dt;
        dA_Y = (AY_w+AY_wvr+AY_v+AY_w_NWFW_tot+LE_AY_wvr+LE_AY_v+LE_AY_w_NWFW_tot)*dt;
        dA_Z = (AZ_w+AZ_wvr+AZ_v+AZ_w_NWFW_tot+LE_AZ_wvr+LE_AZ_v+LE_AZ_w_NWFW_tot)*dt;

        dX_1 = (u1+u1w+u1w_v+u1_w_NWFW_tot+LE_u1w+LE_u1w_v+LE_u1_w_NWFW_tot)*dt;
        dY_1 = (v1+v1w+v1w_v+v1_w_NWFW_tot+LE_v1w+LE_v1w_v+LE_v1_w_NWFW_tot)*dt;
        dZ_1 = (w1+w1w+w1w_v+w1_w_NWFW_tot+LE_w1w+LE_w1w_v+LE_w1_w_NWFW_tot)*dt;

        FW_wake_vortons_buffer{cyc}.A_X = FW_wake_vortons_now_buffer{cyc}.A_X + dA_X ;
        FW_wake_vortons_buffer{cyc}.A_Y = FW_wake_vortons_now_buffer{cyc}.A_Y + dA_Y ;
        FW_wake_vortons_buffer{cyc}.A_Z = FW_wake_vortons_now_buffer{cyc}.A_Z + dA_Z ;
        FW_wake_vortons_buffer{cyc}.X = FW_wake_vortons_now_buffer{cyc}.X + dX_1 ;
        FW_wake_vortons_buffer{cyc}.Y = FW_wake_vortons_now_buffer{cyc}.Y + dY_1 ;
        FW_wake_vortons_buffer{cyc}.Z = FW_wake_vortons_now_buffer{cyc}.Z + dZ_1 ;



end


k = 1;
i = 1;

while i < now_index-2
    for j = 1:cols

    FW_wake_vortons{i,j} = FW_wake_vortons_buffer{k};
    k = k + 1;

    end
    i = i + 1;
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

i = now_index-1;

while i > now_index-3
     for j = 1:n_cols
        X_1 = LE_NW_wake_vortex_rings{i,j}.X_1;
        X_2 = LE_NW_wake_vortex_rings{i,j}.X_2;
        X_3 = LE_NW_wake_vortex_rings{i,j}.X_3;
        X_4 = LE_NW_wake_vortex_rings{i,j}.X_4;
        Y_1 = LE_NW_wake_vortex_rings{i,j}.Y_1;
        Y_2 = LE_NW_wake_vortex_rings{i,j}.Y_2;
        Y_3 = LE_NW_wake_vortex_rings{i,j}.Y_3;
        Y_4 = LE_NW_wake_vortex_rings{i,j}.Y_4;
        Z_1 = LE_NW_wake_vortex_rings{i,j}.Z_1;
        Z_2 = LE_NW_wake_vortex_rings{i,j}.Z_2;
        Z_3 = LE_NW_wake_vortex_rings{i,j}.Z_3;
        Z_4 = LE_NW_wake_vortex_rings{i,j}.Z_4;

        % Loop for inducted velocity of bound vortex rings
u1 = 0;
u2 = 0;
u3 = 0;
u4 = 0;
v1 = 0;
v2 = 0;
v3 = 0;
v4 = 0;
w1 = 0;
w2 = 0;
w3 = 0;
w4 = 0;
u1w = 0;
u2w = 0;
u3w = 0;
u4w = 0;
v1w = 0;
v2w = 0;
v3w = 0;
v4w = 0;
w1w = 0;
w2w = 0;
w3w = 0;
w4w = 0;
u1w_v = 0;
u2w_v = 0;
u3w_v = 0;
u4w_v = 0;
v1w_v = 0;
v2w_v = 0;
v3w_v = 0;
v4w_v = 0;
w1w_v = 0;
w2w_v = 0;
w3w_v = 0;
w4w_v = 0;
LE_u1w = 0;
LE_u2w = 0;
LE_u3w = 0;
LE_u4w = 0;
LE_v1w = 0;
LE_v2w = 0;
LE_v3w = 0;
LE_v4w = 0;
LE_w1w = 0;
LE_w2w = 0;
LE_w3w = 0;
LE_w4w = 0;
LE_u1w_v = 0;
LE_u2w_v = 0;
LE_u3w_v = 0;
LE_u4w_v = 0;
LE_v1w_v = 0;
LE_v2w_v = 0;
LE_v3w_v = 0;
LE_v4w_v = 0;
LE_w1w_v = 0;
LE_w2w_v = 0;
LE_w3w_v = 0;
LE_w4w_v = 0;
u1_w_NWFW_tot = 0;
v1_w_NWFW_tot = 0;
w1_w_NWFW_tot = 0;
u2_w_NWFW_tot = 0;
v2_w_NWFW_tot = 0;
w2_w_NWFW_tot = 0;
u3_w_NWFW_tot = 0;
v3_w_NWFW_tot = 0;
w3_w_NWFW_tot = 0;
u4_w_NWFW_tot = 0;
v4_w_NWFW_tot = 0;
w4_w_NWFW_tot = 0;
LE_u1_w_NWFW_tot = 0;
LE_v1_w_NWFW_tot = 0;
LE_w1_w_NWFW_tot = 0;
LE_u2_w_NWFW_tot = 0;
LE_v2_w_NWFW_tot = 0;
LE_w2_w_NWFW_tot = 0;
LE_u3_w_NWFW_tot = 0;
LE_v3_w_NWFW_tot = 0;
LE_w3_w_NWFW_tot = 0;
LE_u4_w_NWFW_tot = 0;
LE_v4_w_NWFW_tot = 0;
LE_w4_w_NWFW_tot = 0;

%% Loop for wing induction


        for l = 1:n_rows
           for m = 1:n_cols

               [u1_ind_add_R,v1_ind_add_R,w1_ind_add_R] = func_voring_comp(X_1,Y_1,Z_1,vortices_CELL{l,m},...
                   gamma_mat(l,m),0);
               [u2_ind_add_R,v2_ind_add_R,w2_ind_add_R] = func_voring_comp(X_2,Y_2,Z_2,vortices_CELL{l,m},...
                   gamma_mat(l,m),0);
               [u3_ind_add_R,v3_ind_add_R,w3_ind_add_R] = func_voring_comp(X_3,Y_3,Z_3,vortices_CELL{l,m},...
                   gamma_mat(l,m),0);
               [u4_ind_add_R,v4_ind_add_R,w4_ind_add_R] = func_voring_comp(X_4,Y_4,Z_4,vortices_CELL{l,m},...
                   gamma_mat(l,m),0);
               [u1_ind_add_L,v1_ind_add_L,w1_ind_add_L] = func_voring_comp(X_1,-Y_1,Z_1,vortices_CELL{l,m},...
                   gamma_mat(l,m),0);
               [u2_ind_add_L,v2_ind_add_L,w2_ind_add_L] = func_voring_comp(X_2,-Y_2,Z_2,vortices_CELL{l,m},...
                   gamma_mat(l,m),0);
               [u3_ind_add_L,v3_ind_add_L,w3_ind_add_L] = func_voring_comp(X_3,-Y_3,Z_3,vortices_CELL{l,m},...
                   gamma_mat(l,m),0);
               [u4_ind_add_L,v4_ind_add_L,w4_ind_add_L] = func_voring_comp(X_4,-Y_4,Z_4,vortices_CELL{l,m},...
                   gamma_mat(l,m),0);

            u1 = u1 + u1_ind_add_R+ u1_ind_add_L;
            u2 = u2 + u2_ind_add_R+ u2_ind_add_L;
            u3 = u3 + u3_ind_add_R+ u3_ind_add_L;
            u4 = u4 + u4_ind_add_R+ u4_ind_add_L;
            v1 = v1 + v1_ind_add_R- v1_ind_add_L;
            v2 = v2 + v2_ind_add_R- v2_ind_add_L;
            v3 = v3 + v3_ind_add_R- v3_ind_add_L;
            v4 = v4 + v4_ind_add_R- v4_ind_add_L;
            w1 = w1 + w1_ind_add_R+ w1_ind_add_L;
            w2 = w2 + w2_ind_add_R+ w2_ind_add_L;
            w3 = w3 + w3_ind_add_R+ w3_ind_add_L;
            w4 = w4 + w4_ind_add_R+ w4_ind_add_L;
           end
        end

%% Loop for NEAR WAKE wake vortex rings

        r = now_index-1;

        while r > now_index-3
           for s = 1:n_cols

               [u1w_ind_add_R,v1w_ind_add_R,w1w_ind_add_R] = func_voring_comp(X_1,Y_1,Z_1,wake_vortex_ring_now{r,s},...
                   wake_vortex_ring_now{r,s}.GAMMA,0);
               [u2w_ind_add_R,v2w_ind_add_R,w2w_ind_add_R] = func_voring_comp(X_2,Y_2,Z_2,wake_vortex_ring_now{r,s},...
                   wake_vortex_ring_now{r,s}.GAMMA,0);
               [u3w_ind_add_R,v3w_ind_add_R,w3w_ind_add_R] = func_voring_comp(X_3,Y_3,Z_3,wake_vortex_ring_now{r,s},...
                   wake_vortex_ring_now{r,s}.GAMMA,0);
               [u4w_ind_add_R,v4w_ind_add_R,w4w_ind_add_R] = func_voring_comp(X_4,Y_4,Z_4,wake_vortex_ring_now{r,s},...
                   wake_vortex_ring_now{r,s}.GAMMA,0);
               [u1w_ind_add_L,v1w_ind_add_L,w1w_ind_add_L] = func_voring_comp(X_1,-Y_1,Z_1,wake_vortex_ring_now{r,s},...
                   wake_vortex_ring_now{r,s}.GAMMA,0);
               [u2w_ind_add_L,v2w_ind_add_L,w2w_ind_add_L] = func_voring_comp(X_2,-Y_2,Z_2,wake_vortex_ring_now{r,s},...
                   wake_vortex_ring_now{r,s}.GAMMA,0);
               [u3w_ind_add_L,v3w_ind_add_L,w3w_ind_add_L] = func_voring_comp(X_3,-Y_3,Z_3,wake_vortex_ring_now{r,s},...
                   wake_vortex_ring_now{r,s}.GAMMA,0);
               [u4w_ind_add_L,v4w_ind_add_L,w4w_ind_add_L] = func_voring_comp(X_4,-Y_4,Z_4,wake_vortex_ring_now{r,s},...
                   wake_vortex_ring_now{r,s}.GAMMA,0);
            u1w = u1w + u1w_ind_add_R+ u1w_ind_add_L;
            u2w = u2w + u2w_ind_add_R+ u2w_ind_add_L;
            u3w = u3w + u3w_ind_add_R+ u3w_ind_add_L;
            u4w = u4w + u4w_ind_add_R+ u4w_ind_add_L;
            v1w = v1w + v1w_ind_add_R- v1w_ind_add_L;
            v2w = v2w + v2w_ind_add_R- v2w_ind_add_L;
            v3w = v3w + v3w_ind_add_R- v3w_ind_add_L;
            v4w = v4w + v4w_ind_add_R- v4w_ind_add_L;
            w1w = w1w + w1w_ind_add_R+ w1w_ind_add_L;
            w2w = w2w + w2w_ind_add_R+ w2w_ind_add_L;
            w3w = w3w + w3w_ind_add_R+ w3w_ind_add_L;
            w4w = w4w + w4w_ind_add_R+ w4w_ind_add_L;
           end
           r = r - 1;
        end

%% Loop for NEAR WAKE LEADING EDGE wake vortex rings

        r = now_index-1;

        while r > now_index-3
           for s = 1:n_cols

               [LE_u1w_ind_add_R,LE_v1w_ind_add_R,LE_w1w_ind_add_R] = func_voring_comp(X_1,Y_1,Z_1,LE_wake_vortex_ring_now{r,s},...
                   LE_wake_vortex_ring_now{r,s}.GAMMA,0);
               [LE_u2w_ind_add_R,LE_v2w_ind_add_R,LE_w2w_ind_add_R] = func_voring_comp(X_2,Y_2,Z_2,LE_wake_vortex_ring_now{r,s},...
                   LE_wake_vortex_ring_now{r,s}.GAMMA,0);
               [LE_u3w_ind_add_R,LE_v3w_ind_add_R,LE_w3w_ind_add_R] = func_voring_comp(X_3,Y_3,Z_3,LE_wake_vortex_ring_now{r,s},...
                   LE_wake_vortex_ring_now{r,s}.GAMMA,0);
               [LE_u4w_ind_add_R,LE_v4w_ind_add_R,LE_w4w_ind_add_R] = func_voring_comp(X_4,Y_4,Z_4,LE_wake_vortex_ring_now{r,s},...
                   LE_wake_vortex_ring_now{r,s}.GAMMA,0);
               [LE_u1w_ind_add_L,LE_v1w_ind_add_L,LE_w1w_ind_add_L] = func_voring_comp(X_1,-Y_1,Z_1,LE_wake_vortex_ring_now{r,s},...
                   LE_wake_vortex_ring_now{r,s}.GAMMA,0);
               [LE_u2w_ind_add_L,LE_v2w_ind_add_L,LE_w2w_ind_add_L] = func_voring_comp(X_2,-Y_2,Z_2,LE_wake_vortex_ring_now{r,s},...
                   LE_wake_vortex_ring_now{r,s}.GAMMA,0);
               [LE_u3w_ind_add_L,LE_v3w_ind_add_L,LE_w3w_ind_add_L] = func_voring_comp(X_3,-Y_3,Z_3,LE_wake_vortex_ring_now{r,s},...
                   LE_wake_vortex_ring_now{r,s}.GAMMA,0);
               [LE_u4w_ind_add_L,LE_v4w_ind_add_L,LE_w4w_ind_add_L] = func_voring_comp(X_4,-Y_4,Z_4,LE_wake_vortex_ring_now{r,s},...
                   LE_wake_vortex_ring_now{r,s}.GAMMA,0);
            LE_u1w = LE_u1w + LE_u1w_ind_add_R+ LE_u1w_ind_add_L;
            LE_u2w = LE_u2w + LE_u2w_ind_add_R+ LE_u2w_ind_add_L;
            LE_u3w = LE_u3w + LE_u3w_ind_add_R+ LE_u3w_ind_add_L;
            LE_u4w = LE_u4w + LE_u4w_ind_add_R+ LE_u4w_ind_add_L;
            LE_v1w = LE_v1w + LE_v1w_ind_add_R- LE_v1w_ind_add_L;
            LE_v2w = LE_v2w + LE_v2w_ind_add_R- LE_v2w_ind_add_L;
            LE_v3w = LE_v3w + LE_v3w_ind_add_R- LE_v3w_ind_add_L;
            LE_v4w = LE_v4w + LE_v4w_ind_add_R- LE_v4w_ind_add_L;
            LE_w1w = LE_w1w + LE_w1w_ind_add_R+ LE_w1w_ind_add_L;
            LE_w2w = LE_w2w + LE_w2w_ind_add_R+ LE_w2w_ind_add_L;
            LE_w3w = LE_w3w + LE_w3w_ind_add_R+ LE_w3w_ind_add_L;
            LE_w4w = LE_w4w + LE_w4w_ind_add_R+ LE_w4w_ind_add_L;
           end
           r = r - 1;
        end


%% Loop for inductions of FAR WAKE wake vortons

        r = 1;

        while r < now_index-2
           if abs(r-i)<20
           for s = 1:n_cols

               [u1w_vind_add_R,v1w_vind_add_R,w1w_vind_add_R] = func_vortonind(X_1,Y_1,Z_1,FW_wake_vortons_now{r,s},sig);
               [u2w_vind_add_R,v2w_vind_add_R,w2w_vind_add_R] = func_vortonind(X_2,Y_2,Z_2,FW_wake_vortons_now{r,s},sig);
               [u3w_vind_add_R,v3w_vind_add_R,w3w_vind_add_R] = func_vortonind(X_3,Y_3,Z_3,FW_wake_vortons_now{r,s},sig);
               [u4w_vind_add_R,v4w_vind_add_R,w4w_vind_add_R] = func_vortonind(X_4,Y_4,Z_4,FW_wake_vortons_now{r,s},sig);
               [u1w_vind_add_L,v1w_vind_add_L,w1w_vind_add_L] = func_vortonind(X_1,-Y_1,Z_1,FW_wake_vortons_now{r,s},sig);
               [u2w_vind_add_L,v2w_vind_add_L,w2w_vind_add_L] = func_vortonind(X_2,-Y_2,Z_2,FW_wake_vortons_now{r,s},sig);
               [u3w_vind_add_L,v3w_vind_add_L,w3w_vind_add_L] = func_vortonind(X_3,-Y_3,Z_3,FW_wake_vortons_now{r,s},sig);
               [u4w_vind_add_L,v4w_vind_add_L,w4w_vind_add_L] = func_vortonind(X_4,-Y_4,Z_4,FW_wake_vortons_now{r,s},sig);

            u1w_v = u1w_v + u1w_vind_add_R+ u1w_vind_add_L;
            u2w_v = u2w_v + u2w_vind_add_R+ u2w_vind_add_L;
            u3w_v = u3w_v + u3w_vind_add_R+ u3w_vind_add_L;
            u4w_v = u4w_v + u4w_vind_add_R+ u4w_vind_add_L;
            v1w_v = v1w_v + v1w_vind_add_R- v1w_vind_add_L;
            v2w_v = v2w_v + v2w_vind_add_R- v2w_vind_add_L;
            v3w_v = v3w_v + v3w_vind_add_R- v3w_vind_add_L;
            v4w_v = v4w_v + v4w_vind_add_R- v4w_vind_add_L;
            w1w_v = w1w_v + w1w_vind_add_R+ w1w_vind_add_L;
            w2w_v = w2w_v + w2w_vind_add_R+ w2w_vind_add_L;
            w3w_v = w3w_v + w3w_vind_add_R+ w3w_vind_add_L;
            w4w_v = w4w_v + w4w_vind_add_R+ w4w_vind_add_L;
           end
           else
           end
           r = r + 1;
        end

%% Loop for inductions of FAR WAKE LEADING EDGE wake vortons

        r = 1;

        while r < now_index-2
           if abs(r-i)<20
           for s = 1:n_cols

               [LE_u1w_vind_add_R,LE_v1w_vind_add_R,LE_w1w_vind_add_R] = func_vortonind(X_1,Y_1,Z_1,LE_FW_wake_vortons_now{r,s},sig);
               [LE_u2w_vind_add_R,LE_v2w_vind_add_R,LE_w2w_vind_add_R] = func_vortonind(X_2,Y_2,Z_2,LE_FW_wake_vortons_now{r,s},sig);
               [LE_u3w_vind_add_R,LE_v3w_vind_add_R,LE_w3w_vind_add_R] = func_vortonind(X_3,Y_3,Z_3,LE_FW_wake_vortons_now{r,s},sig);
               [LE_u4w_vind_add_R,LE_v4w_vind_add_R,LE_w4w_vind_add_R] = func_vortonind(X_4,Y_4,Z_4,LE_FW_wake_vortons_now{r,s},sig);
               [LE_u1w_vind_add_L,LE_v1w_vind_add_L,LE_w1w_vind_add_L] = func_vortonind(X_1,-Y_1,Z_1,LE_FW_wake_vortons_now{r,s},sig);
               [LE_u2w_vind_add_L,LE_v2w_vind_add_L,LE_w2w_vind_add_L] = func_vortonind(X_2,-Y_2,Z_2,LE_FW_wake_vortons_now{r,s},sig);
               [LE_u3w_vind_add_L,LE_v3w_vind_add_L,LE_w3w_vind_add_L] = func_vortonind(X_3,-Y_3,Z_3,LE_FW_wake_vortons_now{r,s},sig);
               [LE_u4w_vind_add_L,LE_v4w_vind_add_L,LE_w4w_vind_add_L] = func_vortonind(X_4,-Y_4,Z_4,LE_FW_wake_vortons_now{r,s},sig);

            LE_u1w_v = LE_u1w_v + LE_u1w_vind_add_R+ LE_u1w_vind_add_L;
            LE_u2w_v = LE_u2w_v + LE_u2w_vind_add_R+ LE_u2w_vind_add_L;
            LE_u3w_v = LE_u3w_v + LE_u3w_vind_add_R+ LE_u3w_vind_add_L;
            LE_u4w_v = LE_u4w_v + LE_u4w_vind_add_R+ LE_u4w_vind_add_L;
            LE_v1w_v = LE_v1w_v + LE_v1w_vind_add_R- LE_v1w_vind_add_L;
            LE_v2w_v = LE_v2w_v + LE_v2w_vind_add_R- LE_v2w_vind_add_L;
            LE_v3w_v = LE_v3w_v + LE_v3w_vind_add_R- LE_v3w_vind_add_L;
            LE_v4w_v = LE_v4w_v + LE_v4w_vind_add_R- LE_v4w_vind_add_L;
            LE_w1w_v = LE_w1w_v + LE_w1w_vind_add_R+ LE_w1w_vind_add_L;
            LE_w2w_v = LE_w2w_v + LE_w2w_vind_add_R+ LE_w2w_vind_add_L;
            LE_w3w_v = LE_w3w_v + LE_w3w_vind_add_R+ LE_w3w_vind_add_L;
            LE_w4w_v = LE_w4w_v + LE_w4w_vind_add_R+ LE_w4w_vind_add_L;
           end
           else
           end
           r = r + 1;
        end


%% Loop for NEAR-WAKE - FAR-WAKE separation LINE


    for m = 1:n_cols
        X_1SEG = wake_vortex_ring_now{now_index-2,m}.X_3;
        X_2SEG = wake_vortex_ring_now{now_index-2,m}.X_4;
        Y_1SEG = wake_vortex_ring_now{now_index-2,m}.Y_3;
        Y_2SEG = wake_vortex_ring_now{now_index-2,m}.Y_4;
        Z_1SEG = wake_vortex_ring_now{now_index-2,m}.Z_3;
        Z_2SEG = wake_vortex_ring_now{now_index-2,m}.Z_4;
        GAMMASEG = wake_vortex_ring_now{now_index-3,m}.GAMMA;

         [u1_w_NWFW_add_R,v1_w_NWFW_add_R,w1_w_NWFW_add_R] = func_vortexl(X_1,Y_1,Z_1,X_1SEG,Y_1SEG,Z_1SEG,X_2SEG,Y_2SEG,Z_2SEG,GAMMASEG);
         [u1_w_NWFW_add_L,v1_w_NWFW_add_L,w1_w_NWFW_add_L] = func_vortexl(X_1,-Y_1,Z_1,X_1SEG,Y_1SEG,Z_1SEG,X_2SEG,Y_2SEG,Z_2SEG,GAMMASEG);
         u1_w_NWFW_tot = u1_w_NWFW_tot  + u1_w_NWFW_add_R +u1_w_NWFW_add_L;
         v1_w_NWFW_tot = v1_w_NWFW_tot + v1_w_NWFW_add_R -v1_w_NWFW_add_L;
         w1_w_NWFW_tot = w1_w_NWFW_tot + w1_w_NWFW_add_R +w1_w_NWFW_add_L;

         [u2_w_NWFW_add_R,v2_w_NWFW_add_R,w2_w_NWFW_add_R] = func_vortexl(X_2,Y_2,Z_2,X_1SEG,Y_1SEG,Z_1SEG,X_2SEG,Y_2SEG,Z_2SEG,GAMMASEG);
         [u2_w_NWFW_add_L,v2_w_NWFW_add_L,w2_w_NWFW_add_L] = func_vortexl(X_2,-Y_2,Z_2,X_1SEG,Y_1SEG,Z_1SEG,X_2SEG,Y_2SEG,Z_2SEG,GAMMASEG);
         u2_w_NWFW_tot = u2_w_NWFW_tot  + u2_w_NWFW_add_R +u2_w_NWFW_add_L;
         v2_w_NWFW_tot = v2_w_NWFW_tot + v2_w_NWFW_add_R -v2_w_NWFW_add_L;
         w2_w_NWFW_tot = w2_w_NWFW_tot + w2_w_NWFW_add_R +w2_w_NWFW_add_L;

         [u3_w_NWFW_add_R,v3_w_NWFW_add_R,w3_w_NWFW_add_R] = func_vortexl(X_3,Y_3,Z_3,X_1SEG,Y_1SEG,Z_1SEG,X_2SEG,Y_2SEG,Z_2SEG,GAMMASEG);
         [u3_w_NWFW_add_L,v3_w_NWFW_add_L,w3_w_NWFW_add_L] = func_vortexl(X_3,-Y_3,Z_3,X_1SEG,Y_1SEG,Z_1SEG,X_2SEG,Y_2SEG,Z_2SEG,GAMMASEG);
         u3_w_NWFW_tot = u3_w_NWFW_tot  + u3_w_NWFW_add_R +u3_w_NWFW_add_L;
         v3_w_NWFW_tot = v3_w_NWFW_tot + v3_w_NWFW_add_R -v3_w_NWFW_add_L;
         w3_w_NWFW_tot = w3_w_NWFW_tot + w3_w_NWFW_add_R +w3_w_NWFW_add_L;

         [u4_w_NWFW_add_R,v4_w_NWFW_add_R,w4_w_NWFW_add_R] = func_vortexl(X_4,Y_4,Z_4,X_1SEG,Y_1SEG,Z_1SEG,X_2SEG,Y_2SEG,Z_2SEG,GAMMASEG);
         [u4_w_NWFW_add_L,v4_w_NWFW_add_L,w4_w_NWFW_add_L] = func_vortexl(X_4,-Y_4,Z_4,X_1SEG,Y_1SEG,Z_1SEG,X_2SEG,Y_2SEG,Z_2SEG,GAMMASEG);
         u4_w_NWFW_tot = u4_w_NWFW_tot  + u4_w_NWFW_add_R +u4_w_NWFW_add_L;
         v4_w_NWFW_tot = v4_w_NWFW_tot + v4_w_NWFW_add_R -v4_w_NWFW_add_L;
         w4_w_NWFW_tot = w4_w_NWFW_tot + w4_w_NWFW_add_R +w4_w_NWFW_add_L;
    end

%% Loop for NEAR-WAKE - FAR-WAKE LEADING EDGE separation LINE


    for m = 1:n_cols
        X_1SEG = LE_wake_vortex_ring_now{now_index-2,m}.X_3;
        X_2SEG = LE_wake_vortex_ring_now{now_index-2,m}.X_4;
        Y_1SEG = LE_wake_vortex_ring_now{now_index-2,m}.Y_3;
        Y_2SEG = LE_wake_vortex_ring_now{now_index-2,m}.Y_4;
        Z_1SEG = LE_wake_vortex_ring_now{now_index-2,m}.Z_3;
        Z_2SEG = LE_wake_vortex_ring_now{now_index-2,m}.Z_4;
        GAMMASEG = LE_wake_vortex_ring_now{now_index-3,m}.GAMMA;

         [LE_u1_w_NWFW_add_R,LE_v1_w_NWFW_add_R,LE_w1_w_NWFW_add_R] = func_vortexl(X_1,Y_1,Z_1,X_1SEG,Y_1SEG,Z_1SEG,X_2SEG,Y_2SEG,Z_2SEG,GAMMASEG);
         [LE_u1_w_NWFW_add_L,LE_v1_w_NWFW_add_L,LE_w1_w_NWFW_add_L] = func_vortexl(X_1,-Y_1,Z_1,X_1SEG,Y_1SEG,Z_1SEG,X_2SEG,Y_2SEG,Z_2SEG,GAMMASEG);
		 [LE_u2_w_NWFW_add_R,LE_v2_w_NWFW_add_R,LE_w2_w_NWFW_add_R] = func_vortexl(X_2,Y_2,Z_2,X_1SEG,Y_1SEG,Z_1SEG,X_2SEG,Y_2SEG,Z_2SEG,GAMMASEG);
         [LE_u2_w_NWFW_add_L,LE_v2_w_NWFW_add_L,LE_w2_w_NWFW_add_L] = func_vortexl(X_2,-Y_2,Z_2,X_1SEG,Y_1SEG,Z_1SEG,X_2SEG,Y_2SEG,Z_2SEG,GAMMASEG);
		 [LE_u3_w_NWFW_add_R,LE_v3_w_NWFW_add_R,LE_w3_w_NWFW_add_R] = func_vortexl(X_3,Y_3,Z_3,X_1SEG,Y_1SEG,Z_1SEG,X_2SEG,Y_2SEG,Z_2SEG,GAMMASEG);
         [LE_u3_w_NWFW_add_L,LE_v3_w_NWFW_add_L,LE_w3_w_NWFW_add_L] = func_vortexl(X_3,-Y_3,Z_3,X_1SEG,Y_1SEG,Z_1SEG,X_2SEG,Y_2SEG,Z_2SEG,GAMMASEG);
		 [LE_u4_w_NWFW_add_R,LE_v4_w_NWFW_add_R,LE_w4_w_NWFW_add_R] = func_vortexl(X_4,Y_4,Z_4,X_1SEG,Y_1SEG,Z_1SEG,X_2SEG,Y_2SEG,Z_2SEG,GAMMASEG);
         [LE_u4_w_NWFW_add_L,LE_v4_w_NWFW_add_L,LE_w4_w_NWFW_add_L] = func_vortexl(X_4,-Y_4,Z_4,X_1SEG,Y_1SEG,Z_1SEG,X_2SEG,Y_2SEG,Z_2SEG,GAMMASEG);

         LE_u1_w_NWFW_tot =LE_u1_w_NWFW_tot + LE_u1_w_NWFW_add_R +LE_u1_w_NWFW_add_L;
         LE_v1_w_NWFW_tot =LE_v1_w_NWFW_tot + LE_v1_w_NWFW_add_R -LE_v1_w_NWFW_add_L;
         LE_w1_w_NWFW_tot =LE_w1_w_NWFW_tot + LE_w1_w_NWFW_add_R +LE_w1_w_NWFW_add_L;
         LE_u2_w_NWFW_tot =LE_u2_w_NWFW_tot + LE_u2_w_NWFW_add_R +LE_u2_w_NWFW_add_L;
         LE_v2_w_NWFW_tot =LE_v2_w_NWFW_tot + LE_v2_w_NWFW_add_R -LE_v2_w_NWFW_add_L;
         LE_w2_w_NWFW_tot =LE_w2_w_NWFW_tot + LE_w2_w_NWFW_add_R +LE_w2_w_NWFW_add_L;
         LE_u3_w_NWFW_tot =LE_u3_w_NWFW_tot + LE_u3_w_NWFW_add_R +LE_u3_w_NWFW_add_L;
         LE_v3_w_NWFW_tot =LE_v3_w_NWFW_tot + LE_v3_w_NWFW_add_R -LE_v3_w_NWFW_add_L;
         LE_w3_w_NWFW_tot =LE_w3_w_NWFW_tot + LE_w3_w_NWFW_add_R +LE_w3_w_NWFW_add_L;
         LE_u4_w_NWFW_tot =LE_u4_w_NWFW_tot + LE_u4_w_NWFW_add_R +LE_u4_w_NWFW_add_L;
         LE_v4_w_NWFW_tot =LE_v4_w_NWFW_tot + LE_v4_w_NWFW_add_R -LE_v4_w_NWFW_add_L;
         LE_w4_w_NWFW_tot =LE_w4_w_NWFW_tot + LE_w4_w_NWFW_add_R +LE_w4_w_NWFW_add_L;
    end


      

        dX_1 = (u1+u1w+u1w_v+u1_w_NWFW_tot+LE_u1w+LE_u1w_v+LE_u1_w_NWFW_tot)*dt;
        dX_2 = (u2+u2w+u2w_v+u2_w_NWFW_tot+LE_u2w+LE_u2w_v+LE_u2_w_NWFW_tot)*dt;
        dX_3 = (u3+u3w+u3w_v+u3_w_NWFW_tot+LE_u3w+LE_u3w_v+LE_u3_w_NWFW_tot)*dt;
        dX_4 = (u4+u4w+u4w_v+u4_w_NWFW_tot+LE_u4w+LE_u4w_v+LE_u4_w_NWFW_tot)*dt;
        dY_1 = (v1+v1w+v1w_v+v1_w_NWFW_tot+LE_v1w+LE_v1w_v+LE_v1_w_NWFW_tot)*dt;
        dY_2 = (v2+v2w+v2w_v+v2_w_NWFW_tot+LE_v2w+LE_v2w_v+LE_v2_w_NWFW_tot)*dt;
        dY_3 = (v3+v3w+v3w_v+v3_w_NWFW_tot+LE_v3w+LE_v3w_v+LE_v3_w_NWFW_tot)*dt;
        dY_4 = (v4+v4w+v4w_v+v4_w_NWFW_tot+LE_v4w+LE_v4w_v+LE_v4_w_NWFW_tot)*dt;
        dZ_1 = (w1+w1w+w1w_v+w1_w_NWFW_tot+LE_w1w+LE_w1w_v+LE_w1_w_NWFW_tot)*dt;
        dZ_2 = (w2+w2w+w2w_v+w2_w_NWFW_tot+LE_w2w+LE_w2w_v+LE_w2_w_NWFW_tot)*dt;
        dZ_3 = (w3+w3w+w3w_v+w3_w_NWFW_tot+LE_w3w+LE_w3w_v+LE_w3_w_NWFW_tot)*dt;
        dZ_4 = (w4+w4w+w4w_v+w4_w_NWFW_tot+LE_w4w+LE_w4w_v+LE_w4_w_NWFW_tot)*dt;





        LE_NW_wake_vortex_rings{i,j}.X_1 = LE_wake_vortex_ring_now{i,j}.X_1 + dX_1 ;
        LE_NW_wake_vortex_rings{i,j}.X_2 = LE_wake_vortex_ring_now{i,j}.X_2 + dX_2 ;
        LE_NW_wake_vortex_rings{i,j}.X_3 = LE_wake_vortex_ring_now{i,j}.X_3 + dX_3 ;
        LE_NW_wake_vortex_rings{i,j}.X_4 = LE_wake_vortex_ring_now{i,j}.X_4 + dX_4 ;
        LE_NW_wake_vortex_rings{i,j}.Y_1 = LE_wake_vortex_ring_now{i,j}.Y_1 + dY_1 ;
        LE_NW_wake_vortex_rings{i,j}.Y_2 = LE_wake_vortex_ring_now{i,j}.Y_2 + dY_2 ;
        LE_NW_wake_vortex_rings{i,j}.Y_3 = LE_wake_vortex_ring_now{i,j}.Y_3 + dY_3 ;
        LE_NW_wake_vortex_rings{i,j}.Y_4 = LE_wake_vortex_ring_now{i,j}.Y_4 + dY_4 ;
        LE_NW_wake_vortex_rings{i,j}.Z_1 = LE_wake_vortex_ring_now{i,j}.Z_1 + dZ_1 ;
        LE_NW_wake_vortex_rings{i,j}.Z_2 = LE_wake_vortex_ring_now{i,j}.Z_2 + dZ_2 ;
        LE_NW_wake_vortex_rings{i,j}.Z_3 = LE_wake_vortex_ring_now{i,j}.Z_3 + dZ_3 ;
        LE_NW_wake_vortex_rings{i,j}.Z_4 = LE_wake_vortex_ring_now{i,j}.Z_4 + dZ_4 ;

     end


i = i - 1;
end







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




cols = n_cols;
rows = n_rows;

LE_FW_wake_vortons_buffer = cell(1,cols*(now_index-3));
LE_FW_wake_vortons_now_buffer = cell(1,cols*(now_index-3));

k = 1;
i = 1;

while i < now_index-2
    for j = 1:cols

   LE_FW_wake_vortons_buffer{k} = LE_FW_wake_vortons{i,j};
   LE_FW_wake_vortons_now_buffer{k}    = LE_FW_wake_vortons_now{i,j};
   k = k + 1;
    end
   i = i + 1;
end



parfor cyc = 1:cols*(now_index-3)

        X_1 = LE_FW_wake_vortons_buffer{cyc}.X;
        Y_1 = LE_FW_wake_vortons_buffer{cyc}.Y;
        Z_1 = LE_FW_wake_vortons_buffer{cyc}.Z;
        alfa_x = LE_FW_wake_vortons_buffer{cyc}.A_X;
        alfa_y = LE_FW_wake_vortons_buffer{cyc}.A_Y;
        alfa_z = LE_FW_wake_vortons_buffer{cyc}.A_Z;

        alfa_vec = [alfa_x,alfa_y,alfa_z];

AX_w = 0;
AY_w = 0;
AZ_w = 0;
AX_wvr = 0;
AY_wvr = 0;
AZ_wvr = 0;
AX_v = 0;
AY_v = 0;
AZ_v = 0;
AX_w_NWFW_tot = 0;
AY_w_NWFW_tot = 0;
AZ_w_NWFW_tot = 0;
LE_AX_wvr = 0;
LE_AY_wvr = 0;
LE_AZ_wvr = 0;
LE_AX_v = 0;
LE_AY_v = 0;
LE_AZ_v = 0;
LE_AX_w_NWFW_tot = 0;
LE_AY_w_NWFW_tot = 0;
LE_AZ_w_NWFW_tot = 0;
u1 = 0;
v1 = 0;
w1 = 0;
u1w = 0;
v1w = 0;
w1w = 0;
u1w_v = 0;
v1w_v = 0;
w1w_v = 0;
u1_w_NWFW_tot = 0;
v1_w_NWFW_tot = 0;
w1_w_NWFW_tot = 0;
LE_u1w = 0;
LE_v1w = 0;
LE_w1w = 0;
LE_u1w_v = 0;
LE_v1w_v = 0;
LE_w1w_v = 0;
LE_u1_w_NWFW_tot = 0;
LE_v1_w_NWFW_tot = 0;
LE_w1_w_NWFW_tot = 0;

%% Loop for inducted velocity of bound vortex rings

        for l = 1:rows
           for m = 1:cols

               [u1_ind_add_R,v1_ind_add_R,w1_ind_add_R] = func_voring_comp(X_1,Y_1,Z_1,vortices_CELL{l,m},...
                   gamma_mat(l,m),0);
               [u1_ind_add_L,v1_ind_add_L,w1_ind_add_L] = func_voring_comp(X_1,-Y_1,Z_1,vortices_CELL{l,m},...
                   gamma_mat(l,m),0);


            u1 = u1 + u1_ind_add_R+ u1_ind_add_L;
            v1 = v1 + v1_ind_add_R- v1_ind_add_L;
            w1 = w1 + w1_ind_add_R+ w1_ind_add_L;

           end
        end

%% Loop for inducted velocity of leading edge wake vortex rings
        r = now_index-1;

        while r > now_index-3
           for s = 1:cols

               [LE_u1w_ind_add_R,LE_v1w_ind_add_R,LE_w1w_ind_add_R] = func_voring_comp(X_1,Y_1,Z_1,LE_wake_vortex_ring_now{r,s},...
                   LE_wake_vortex_ring_now{r,s}.GAMMA,0);
               [LE_u1w_ind_add_L,LE_v1w_ind_add_L,LE_w1w_ind_add_L] = func_voring_comp(X_1,-Y_1,Z_1,LE_wake_vortex_ring_now{r,s},...
                   LE_wake_vortex_ring_now{r,s}.GAMMA,0);

            LE_u1w = LE_u1w + LE_u1w_ind_add_R+ LE_u1w_ind_add_L;
            LE_v1w = LE_v1w + LE_v1w_ind_add_R- LE_v1w_ind_add_L;
            LE_w1w = LE_w1w + LE_w1w_ind_add_R+ LE_w1w_ind_add_L;

           end
           r = r - 1;
        end

%% Loop for inducted velocity of wake vortex rings
        r = now_index-1;

        while r > now_index-3
           for s = 1:cols

               [u1w_ind_add_R,v1w_ind_add_R,w1w_ind_add_R] = func_voring_comp(X_1,Y_1,Z_1,wake_vortex_ring_now{r,s},...
                   wake_vortex_ring_now{r,s}.GAMMA,0);
               [u1w_ind_add_L,v1w_ind_add_L,w1w_ind_add_L] = func_voring_comp(X_1,-Y_1,Z_1,wake_vortex_ring_now{r,s},...
                   wake_vortex_ring_now{r,s}.GAMMA,0);

            u1w = u1w + u1w_ind_add_R+ u1w_ind_add_L;
            v1w = v1w + v1w_ind_add_R- v1w_ind_add_L;
            w1w = w1w + w1w_ind_add_R+ w1w_ind_add_L;

           end
           r = r - 1;
        end


%% Loop for inducted velocity of wake vortons
        r = 1;

        while r < now_index-2
           for s = 1:cols


               [u1w_vind_add_R,v1w_vind_add_R,w1w_vind_add_R] = func_vortonind(X_1,Y_1,Z_1,FW_wake_vortons_now{r,s},sig);
               [u1w_vind_add_L,v1w_vind_add_L,w1w_vind_add_L] = func_vortonind(X_1,-Y_1,Z_1,FW_wake_vortons_now{r,s},sig);


            u1w_v = u1w_v + u1w_vind_add_R+ u1w_vind_add_L;
            v1w_v = v1w_v + v1w_vind_add_R- v1w_vind_add_L;
            w1w_v = w1w_v + w1w_vind_add_R+ w1w_vind_add_L;

           end
           r = r + 1;
        end

 %% Loop for inducted velocity of LE wake vortons
        r = 1;

        while r < now_index-2
           for s = 1:cols


               [LE_u1w_vind_add_R,LE_v1w_vind_add_R,LE_w1w_vind_add_R] = func_vortonind(X_1,Y_1,Z_1,LE_FW_wake_vortons_now{r,s},sig);
               [LE_u1w_vind_add_L,LE_v1w_vind_add_L,LE_w1w_vind_add_L] = func_vortonind(X_1,-Y_1,Z_1,LE_FW_wake_vortons_now{r,s},sig);


            LE_u1w_v = LE_u1w_v + LE_u1w_vind_add_R+ LE_u1w_vind_add_L;
            LE_v1w_v = LE_v1w_v + LE_v1w_vind_add_R- LE_v1w_vind_add_L;
            LE_w1w_v = LE_w1w_v + LE_w1w_vind_add_R+ LE_w1w_vind_add_L;

           end
           r = r + 1;
        end

%% Loop for NEAR-WAKE - FAR-WAKE separation line

    for m = 1:cols
        X_1SEG = wake_vortex_ring_now{now_index-2,m}.X_3;
        X_2SEG = wake_vortex_ring_now{now_index-2,m}.X_4;
        Y_1SEG = wake_vortex_ring_now{now_index-2,m}.Y_3;
        Y_2SEG = wake_vortex_ring_now{now_index-2,m}.Y_4;
        Z_1SEG = wake_vortex_ring_now{now_index-2,m}.Z_3;
        Z_2SEG = wake_vortex_ring_now{now_index-2,m}.Z_4;
        GAMMASEG = wake_vortex_ring_now{now_index-3,m}.GAMMA;

         [u1_w_NWFW_add_R,v1_w_NWFW_add_R,w1_w_NWFW_add_R] = func_vortexl(X_1,Y_1,Z_1,X_1SEG,Y_1SEG,Z_1SEG,X_2SEG,Y_2SEG,Z_2SEG,GAMMASEG);
         [u1_w_NWFW_add_L,v1_w_NWFW_add_L,w1_w_NWFW_add_L] = func_vortexl(X_1,-Y_1,Z_1,X_1SEG,Y_1SEG,Z_1SEG,X_2SEG,Y_2SEG,Z_2SEG,GAMMASEG);
         u1_w_NWFW_tot = u1_w_NWFW_tot + u1_w_NWFW_add_R +u1_w_NWFW_add_L;
         v1_w_NWFW_tot = v1_w_NWFW_tot + v1_w_NWFW_add_R -v1_w_NWFW_add_L;
         w1_w_NWFW_tot = w1_w_NWFW_tot + w1_w_NWFW_add_R +w1_w_NWFW_add_L;
    end

%% Loop for NEAR-WAKE - FAR-WAKE LEADING-EDGE separation line

    for m = 1:cols
        X_1SEG = LE_wake_vortex_ring_now{now_index-2,m}.X_3;
        X_2SEG = LE_wake_vortex_ring_now{now_index-2,m}.X_4;
        Y_1SEG = LE_wake_vortex_ring_now{now_index-2,m}.Y_3;
        Y_2SEG = LE_wake_vortex_ring_now{now_index-2,m}.Y_4;
        Z_1SEG = LE_wake_vortex_ring_now{now_index-2,m}.Z_3;
        Z_2SEG = LE_wake_vortex_ring_now{now_index-2,m}.Z_4;
        GAMMASEG = LE_wake_vortex_ring_now{now_index-3,m}.GAMMA;

        [LE_u1_w_NWFW_add_R,LE_v1_w_NWFW_add_R,LE_w1_w_NWFW_add_R] = func_vortexl(X_1,Y_1,Z_1,X_1SEG,Y_1SEG,Z_1SEG,X_2SEG,Y_2SEG,Z_2SEG,GAMMASEG);
        [LE_u1_w_NWFW_add_L,LE_v1_w_NWFW_add_L,LE_w1_w_NWFW_add_L] = func_vortexl(X_1,-Y_1,Z_1,X_1SEG,Y_1SEG,Z_1SEG,X_2SEG,Y_2SEG,Z_2SEG,GAMMASEG);
         LE_u1_w_NWFW_tot = LE_u1_w_NWFW_tot + LE_u1_w_NWFW_add_R +LE_u1_w_NWFW_add_L;
         LE_v1_w_NWFW_tot = LE_v1_w_NWFW_tot + LE_v1_w_NWFW_add_R -LE_v1_w_NWFW_add_L;
         LE_w1_w_NWFW_tot = LE_w1_w_NWFW_tot + LE_w1_w_NWFW_add_R +LE_w1_w_NWFW_add_L;
    end

%% Loop for inducted velocity of bound vortex rings


        for l = 1:rows
           for m = 1:cols

               vel_grad_mat_R = func_deriv_voring(X_1,Y_1,Z_1,vortices_CELL{l,m},...
                   gamma_mat(l,m));
               vel_grad_mat_L = func_deriv_voring(X_1,-Y_1,Z_1,vortices_CELL{l,m},...
                   gamma_mat(l,m));

            var_vect_R = alfa_vec * vel_grad_mat_R;
            var_vect_L = alfa_vec * vel_grad_mat_L;

            u1_ind_add_R = var_vect_R(1);
            v1_ind_add_R = var_vect_R(2);
            w1_ind_add_R = var_vect_R(3);

            u1_ind_add_L = var_vect_L(1);
            v1_ind_add_L = var_vect_L(2);
            w1_ind_add_L = var_vect_L(3);

            AX_w = AX_w + u1_ind_add_R+ u1_ind_add_L;
            AY_w = AY_w + v1_ind_add_R- v1_ind_add_L;
            AZ_w = AZ_w + w1_ind_add_R+ w1_ind_add_L;

           end
        end

%% Loop for inducted velocity of wake vortex rings

        r = now_index-1;

        while r > now_index-3
           for s = 1:cols

                vel_grad_mat_R = func_deriv_voring(X_1,Y_1,Z_1,wake_vortex_ring_now{r,s},...
                   wake_vortex_ring_now{r,s}.GAMMA);
                vel_grad_mat_L = func_deriv_voring(X_1,-Y_1,Z_1,wake_vortex_ring_now{r,s},...
                   wake_vortex_ring_now{r,s}.GAMMA);

            var_vect_R = alfa_vec * vel_grad_mat_R;
            var_vect_L = alfa_vec * vel_grad_mat_L;

            u1w_ind_add_R = var_vect_R(1);
            v1w_ind_add_R = var_vect_R(2);
            w1w_ind_add_R = var_vect_R(3);

            u1w_ind_add_L = var_vect_L(1);
            v1w_ind_add_L = var_vect_L(2);
            w1w_ind_add_L = var_vect_L(3);

            AX_wvr = AX_wvr + u1w_ind_add_R+ u1w_ind_add_L;
            AY_wvr = AY_wvr + v1w_ind_add_R- v1w_ind_add_L;
            AZ_wvr = AZ_wvr + w1w_ind_add_R+ w1w_ind_add_L;

           end
           r = r - 1;
        end
%% Loop for inducted velocity of LE wake vortex rings

        r = now_index-1;

        while r > now_index-3
           for s = 1:cols

                vel_grad_mat_R = func_deriv_voring(X_1,Y_1,Z_1,LE_wake_vortex_ring_now{r,s},...
                   LE_wake_vortex_ring_now{r,s}.GAMMA);
                vel_grad_mat_L = func_deriv_voring(X_1,-Y_1,Z_1,LE_wake_vortex_ring_now{r,s},...
                   LE_wake_vortex_ring_now{r,s}.GAMMA);

            var_vect_R = alfa_vec * vel_grad_mat_R;
            var_vect_L = alfa_vec * vel_grad_mat_L;

            u1w_ind_add_R = var_vect_R(1);
            v1w_ind_add_R = var_vect_R(2);
            w1w_ind_add_R = var_vect_R(3);

            u1w_ind_add_L = var_vect_L(1);
            v1w_ind_add_L = var_vect_L(2);
            w1w_ind_add_L = var_vect_L(3);

            LE_AX_wvr = LE_AX_wvr + u1w_ind_add_R+ u1w_ind_add_L;
            LE_AY_wvr = LE_AY_wvr + v1w_ind_add_R- v1w_ind_add_L;
            LE_AZ_wvr = LE_AZ_wvr + w1w_ind_add_R+ w1w_ind_add_L;

           end
           r = r - 1;
        end


%% Loop for inducted velocity of wake vortons

        r = 1;

        while r < now_index-2
           for s = 1:cols


               [u1w_vind_add_R,v1w_vind_add_R,w1w_vind_add_R] = func_vortonstre(LE_FW_wake_vortons_buffer{cyc},FW_wake_vortons_now{r,s},sig);
               [u1w_vind_add_L,v1w_vind_add_L,w1w_vind_add_L] = func_vortonstre(LE_FW_wake_vortons_buffer{cyc},FW_wake_vortons_now{r,s},sig);


            AX_v = AX_v + u1w_vind_add_R+ u1w_vind_add_L;
            AY_v = AY_v + v1w_vind_add_R- v1w_vind_add_L;
            AZ_v = AZ_v + w1w_vind_add_R+ w1w_vind_add_L;

           end
           r = r + 1;
        end

%% Loop for inducted velocity of wake vortons

        r = 1;

        while r < now_index-2
           for s = 1:cols


               [u1w_vind_add_R,v1w_vind_add_R,w1w_vind_add_R] = func_vortonstre(LE_FW_wake_vortons_buffer{cyc},LE_FW_wake_vortons_now{r,s},sig);
               [u1w_vind_add_L,v1w_vind_add_L,w1w_vind_add_L] = func_vortonstre(LE_FW_wake_vortons_buffer{cyc},LE_FW_wake_vortons_now{r,s},sig);


            LE_AX_v = LE_AX_v + u1w_vind_add_R+ u1w_vind_add_L;
            LE_AY_v = LE_AY_v + v1w_vind_add_R- v1w_vind_add_L;
            LE_AZ_v = LE_AZ_v + w1w_vind_add_R+ w1w_vind_add_L;

           end
           r = r + 1;
        end

%% Loop for NEAR-WAKE - FAR-WAKE separation line

    for m = 1:cols
        X_1SEG = wake_vortex_ring_now{now_index-2,m}.X_3;
        X_2SEG = wake_vortex_ring_now{now_index-2,m}.X_4;
        Y_1SEG = wake_vortex_ring_now{now_index-2,m}.Y_3;
        Y_2SEG = wake_vortex_ring_now{now_index-2,m}.Y_4;
        Z_1SEG = wake_vortex_ring_now{now_index-2,m}.Z_3;
        Z_2SEG = wake_vortex_ring_now{now_index-2,m}.Z_4;
        GAMMASEG = wake_vortex_ring_now{now_index-3,m}.GAMMA;

         [dudx_val_R,dudy_val_R,dudz_val_R,dvdx_val_R,dvdy_val_R,dvdz_val_R,dwdx_val_R,dwdy_val_R,dwdz_val_R] = func_deriv_vortexl(X_1,Y_1,Z_1,X_1SEG,Y_1SEG,Z_1SEG,X_2SEG,Y_2SEG,Z_2SEG,GAMMASEG);
         [dudx_val_L,dudy_val_L,dudz_val_L,dvdx_val_L,dvdy_val_L,dvdz_val_L,dwdx_val_L,dwdy_val_L,dwdz_val_L] = func_deriv_vortexl(X_1,-Y_1,Z_1,X_1SEG,Y_1SEG,Z_1SEG,X_2SEG,Y_2SEG,Z_2SEG,GAMMASEG);

            var_vect_R = alfa_vec * [dudx_val_R,dudy_val_R,dudz_val_R;...
                                    dvdx_val_R,dvdy_val_R,dvdz_val_R;...
                                    dwdx_val_R,dwdy_val_R,dwdz_val_R];

            var_vect_L = alfa_vec * [dudx_val_L,dudy_val_L,dudz_val_L;...
                                    dvdx_val_L,dvdy_val_L,dvdz_val_L;...
                                    dwdx_val_L,dwdy_val_L,dwdz_val_L];

            u1_w_NWFW_add_R = var_vect_R(1);
            v1_w_NWFW_add_R = var_vect_R(2);
            w1_w_NWFW_add_R = var_vect_R(3);

            u1_w_NWFW_add_L = var_vect_L(1);
            v1_w_NWFW_add_L = var_vect_L(2);
            w1_w_NWFW_add_L = var_vect_L(3);

         AX_w_NWFW_tot = AX_w_NWFW_tot + u1_w_NWFW_add_R +u1_w_NWFW_add_L;
         AY_w_NWFW_tot = AY_w_NWFW_tot + v1_w_NWFW_add_R -v1_w_NWFW_add_L;
         AZ_w_NWFW_tot = AZ_w_NWFW_tot + w1_w_NWFW_add_R +w1_w_NWFW_add_L;
    end

%% Loop for LE NEAR-WAKE - FAR-WAKE separation line

    for m = 1:cols
        X_1SEG = LE_wake_vortex_ring_now{now_index-2,m}.X_3;
        X_2SEG = LE_wake_vortex_ring_now{now_index-2,m}.X_4;
        Y_1SEG = LE_wake_vortex_ring_now{now_index-2,m}.Y_3;
        Y_2SEG = LE_wake_vortex_ring_now{now_index-2,m}.Y_4;
        Z_1SEG = LE_wake_vortex_ring_now{now_index-2,m}.Z_3;
        Z_2SEG = LE_wake_vortex_ring_now{now_index-2,m}.Z_4;
        GAMMASEG = LE_wake_vortex_ring_now{now_index-3,m}.GAMMA;

         [dudx_val_R,dudy_val_R,dudz_val_R,dvdx_val_R,dvdy_val_R,dvdz_val_R,dwdx_val_R,dwdy_val_R,dwdz_val_R] = func_deriv_vortexl(X_1,Y_1,Z_1,X_1SEG,Y_1SEG,Z_1SEG,X_2SEG,Y_2SEG,Z_2SEG,GAMMASEG);
         [dudx_val_L,dudy_val_L,dudz_val_L,dvdx_val_L,dvdy_val_L,dvdz_val_L,dwdx_val_L,dwdy_val_L,dwdz_val_L] = func_deriv_vortexl(X_1,-Y_1,Z_1,X_1SEG,Y_1SEG,Z_1SEG,X_2SEG,Y_2SEG,Z_2SEG,GAMMASEG);

            var_vect_R = alfa_vec * [dudx_val_R,dudy_val_R,dudz_val_R;...
                                    dvdx_val_R,dvdy_val_R,dvdz_val_R;...
                                    dwdx_val_R,dwdy_val_R,dwdz_val_R];

            var_vect_L = alfa_vec * [dudx_val_L,dudy_val_L,dudz_val_L;...
                                    dvdx_val_L,dvdy_val_L,dvdz_val_L;...
                                    dwdx_val_L,dwdy_val_L,dwdz_val_L];

            u1_w_NWFW_add_R = var_vect_R(1);
            v1_w_NWFW_add_R = var_vect_R(2);
            w1_w_NWFW_add_R = var_vect_R(3);

            u1_w_NWFW_add_L = var_vect_L(1);
            v1_w_NWFW_add_L = var_vect_L(2);
            w1_w_NWFW_add_L = var_vect_L(3);

         LE_AX_w_NWFW_tot = LE_AX_w_NWFW_tot + u1_w_NWFW_add_R +u1_w_NWFW_add_L;
         LE_AY_w_NWFW_tot = LE_AY_w_NWFW_tot + v1_w_NWFW_add_R -v1_w_NWFW_add_L;
         LE_AZ_w_NWFW_tot = LE_AZ_w_NWFW_tot + w1_w_NWFW_add_R +w1_w_NWFW_add_L;
    end
%% VORTONS UPDATING

        dA_X = (AX_w+AX_wvr+AX_v+AX_w_NWFW_tot+LE_AX_wvr+LE_AX_v+LE_AX_w_NWFW_tot)*dt;
        dA_Y = (AY_w+AY_wvr+AY_v+AY_w_NWFW_tot+LE_AY_wvr+LE_AY_v+LE_AY_w_NWFW_tot)*dt;
        dA_Z = (AZ_w+AZ_wvr+AZ_v+AZ_w_NWFW_tot+LE_AZ_wvr+LE_AZ_v+LE_AZ_w_NWFW_tot)*dt;


        dX_1 = (u1+u1w+u1w_v+u1_w_NWFW_tot+LE_u1w+LE_u1w_v+LE_u1_w_NWFW_tot)*dt;
        dY_1 = (v1+v1w+v1w_v+v1_w_NWFW_tot+LE_v1w+LE_v1w_v+LE_v1_w_NWFW_tot)*dt;
        dZ_1 = (w1+w1w+w1w_v+w1_w_NWFW_tot+LE_w1w+LE_w1w_v+LE_w1_w_NWFW_tot)*dt;



        LE_FW_wake_vortons_buffer{cyc}.X = LE_FW_wake_vortons_now_buffer{cyc}.X + dX_1 ;
        LE_FW_wake_vortons_buffer{cyc}.Y = LE_FW_wake_vortons_now_buffer{cyc}.Y + dY_1 ;
        LE_FW_wake_vortons_buffer{cyc}.Z = LE_FW_wake_vortons_now_buffer{cyc}.Z + dZ_1
        LE_FW_wake_vortons_buffer{cyc}.A_X = LE_FW_wake_vortons_now_buffer{cyc}.A_X + dA_X ;
        LE_FW_wake_vortons_buffer{cyc}.A_Y = LE_FW_wake_vortons_now_buffer{cyc}.A_Y + dA_Y ;
        LE_FW_wake_vortons_buffer{cyc}.A_Z = LE_FW_wake_vortons_now_buffer{cyc}.A_Z + dA_Z ;

end


k = 1;
i = 1;

while i < now_index-2
    for j = 1:cols

    LE_FW_wake_vortons{i,j} = LE_FW_wake_vortons_buffer{k};
    k = k + 1;

    end
    i = i + 1;
end

end
