function  [X_NODE,Y_NODE,Z_NODE] = func_INPUT(C_F_mat,vortices_mat)
%func_output The function creates an output file to describe forces and
%moments on the undeformed wing
%   Detailed explanation goes here


global n_rows n_cols

 %% To Plot wing with forces



k1 = 1;

for i = 1:n_rows
    for j = 1:n_cols
        XC(i,j) = (vortices_mat{i,j}.X_C);
        YC(i,j) = (vortices_mat{i,j}.Y_C);
        ZC(i,j) = (vortices_mat{i,j}.Z_C);
        CFX(i,j)  = C_F_mat{i,j}(1);
        CFY(i,j)  = C_F_mat{i,j}(2);
        CFZ(i,j)  = C_F_mat{i,j}(3);
    X_1(k1) = vortices_mat{i,j}.X_1;
    X_2(k1) = vortices_mat{i,j}.X_2;
    X_3(k1) = vortices_mat{i,j}.X_3;
    X_4(k1) = vortices_mat{i,j}.X_4;
    Y_1(k1) = vortices_mat{i,j}.Y_1;
    Y_2(k1) = vortices_mat{i,j}.Y_2;
    Y_3(k1) = vortices_mat{i,j}.Y_3;
    Y_4(k1) = vortices_mat{i,j}.Y_4;
    Z_1(k1) = vortices_mat{i,j}.Z_1;
    Z_2(k1) = vortices_mat{i,j}.Z_2;
    Z_3(k1) = vortices_mat{i,j}.Z_3;
    Z_4(k1) = vortices_mat{i,j}.Z_4;
    
k1 = k1+1;
    end
end

k = 1;
for i = 1:n_rows+1
    for j = 1:n_cols+1
    
    NODE_LABEL(k) = i*2+j*200;
    
    
    
    if i == n_rows + 1
        if j == n_cols + 1
            
        X_NODE(k) = vortices_mat{i-1,j-1}.X_4;
        Y_NODE(k) = vortices_mat{i-1,j-1}.Y_4;
        Z_NODE(k) = vortices_mat{i-1,j-1}.Z_4;
        CFX_NODE(k)    =  0;
        CFY_NODE(k)    =  0;
        CFZ_NODE(k)    =  0;
        
        
        else
        if j == 1
        X_NODE(k) = vortices_mat{i-1,j}.X_3;
        Y_NODE(k) = vortices_mat{i-1,j}.Y_3;
        Z_NODE(k) = vortices_mat{i-1,j}.Z_3;
        CFX_NODE(k) = 0;
        CFY_NODE(k) = 0;
        CFZ_NODE(k) = 0;    
        else    
        X_NODE(k) = vortices_mat{i-1,j}.X_3;
        Y_NODE(k) = vortices_mat{i-1,j}.Y_3;
        Z_NODE(k) = vortices_mat{i-1,j}.Z_3;
        CFX_NODE(k) = 0;
        CFY_NODE(k) = 0;
        CFZ_NODE(k) = 0;
        end
        end
    else
    if j == n_cols+1
        X_NODE(k) = vortices_mat{i,j-1}.X_2;
        Y_NODE(k) = vortices_mat{i,j-1}.Y_2;
        Z_NODE(k) = vortices_mat{i,j-1}.Z_2;
        CFX_NODE(k)    =  0;
        CFY_NODE(k)    =  0;
        CFZ_NODE(k)    =  0;
        
    else
        if j == 1
        X_NODE(k) = vortices_mat{i,j}.X_1;
        Y_NODE(k) = vortices_mat{i,j}.Y_1;
        Z_NODE(k) = vortices_mat{i,j}.Z_1;
        CFX_NODE(k) = interp1(XC(:,j),CFX(:,j),X_NODE(k),'linear','extrap');
        CFY_NODE(k) = interp1(XC(:,j),CFY(:,j),X_NODE(k),'linear','extrap');
        CFZ_NODE(k) = interp1(XC(:,j),CFZ(:,j),X_NODE(k),'linear','extrap');    
        else
        X_NODE(k) = vortices_mat{i,j}.X_1;
        Y_NODE(k) = vortices_mat{i,j}.Y_1;
        Z_NODE(k) = vortices_mat{i,j}.Z_1;
        CFX_NODE(k) = (interp1(XC(:,j),CFX(:,j),X_NODE(k),'linear','extrap') + interp1(XC(:,j-1),CFX(:,j-1),X_NODE(k),'linear','extrap'))/2;
        CFY_NODE(k) = (interp1(XC(:,j),CFY(:,j),X_NODE(k),'linear','extrap') + interp1(XC(:,j-1),CFY(:,j-1),X_NODE(k),'linear','extrap'))/2;
        CFZ_NODE(k) = (interp1(XC(:,j),CFZ(:,j),X_NODE(k),'linear','extrap') + interp1(XC(:,j-1),CFZ(:,j-1),X_NODE(k),'linear','extrap'))/2;
        end
    end
    end
    
        MX_NODE(k) = 0;
        MY_NODE(k) = 0;
        MZ_NODE(k) = 0;
    k = k + 1;
    end
        
end


%  drawnow
%  figure()
%  patch([X_1;X_2;X_4;X_3],[Y_1;Y_2;Y_4;Y_3],[Z_1;Z_2;Z_4;Z_3],[1,1,1]);hold on;
%  patch([X_1;X_2;X_4;X_3],-[Y_1;Y_2;Y_4;Y_3],[Z_1;Z_2;Z_4;Z_3],[1,1,1]);hold on;
%  quiver3(X_NODE,Y_NODE,Z_NODE,CFX_NODE,CFY_NODE,CFZ_NODE,'k');hold on;
%  quiver3(X_NODE,-Y_NODE,Z_NODE,CFX_NODE,-CFY_NODE,CFZ_NODE,'k');hold on;
%  xlabel('X');ylabel('Y');zlabel('Z');
%  
 A = [NODE_LABEL;CFX_NODE;CFY_NODE;CFZ_NODE;MX_NODE;MY_NODE;MZ_NODE];
 
fileID = fopen('INPUT.dat','w');
fprintf(fileID,'%1$s %2$s %3$s %4$s %5$s %6$s %7$s \n','#Node_Label','Fx','Fy','Fz','Mx','My','Mz');
fprintf(fileID,'%4u %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f\n',A);
fclose(fileID);
end

