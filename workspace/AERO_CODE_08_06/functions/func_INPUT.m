function  func_INPUT(F_mat,vortices_mat)
%func_output The function creates an output file to describe forces and
%moments on the undeformed wing
%   Detailed explanation goes here


global n_rows n_cols rho


k1 = 1;

for j = 1:n_cols

	    NODE_LABEL(k1) = j*200;
		XP(j) = ((vortices_mat{1,j}.X_1) + (vortices_mat{n_rows,j}.X_3))/2;
        YP(j) = ((vortices_mat{1,j}.Y_1) + (vortices_mat{n_rows,j}.Y_2))/2;
        ZP(j) = ((vortices_mat{1,j}.Z_1) + (vortices_mat{n_rows,j}.Z_2))/2;
		FX(k1) = 0;
        FY(k1) = 0;
        FZ(k1) = 0;
		MX(k1) = 0;
		MY(k1) = 0;
		MZ(k1) = 0;
		
    for i = 1:n_rows
	

	
        RXCG = (vortices_mat{i,j}.X_C)-XP(j);
        RYCG = (vortices_mat{i,j}.Y_C)-YP(j);
        RZCG = (vortices_mat{i,j}.Z_C)-ZP(j);
        RXCN = RYCG;
		RYCN = -RXCG;
		RZCN = RZCG;
		
		
		FX(k1)  = FX(k1) + rho*F_mat{i,j}(2);
        FY(k1)  = FY(k1)-rho*F_mat{i,j}(1);
        FZ(k1)  = FZ(k1) + rho*F_mat{i,j}(3);
		
	    MX(k1)  = MX(k1) + RYCN*rho*F_mat{i,j}(3) + RZCN *rho*F_mat{i,j}(1);
		MY(k1)  = MY(k1) - RXCN*rho*F_mat{i,j}(3) + RZCN *rho*F_mat{i,j}(2);
 		MZ(k1)  = MZ(k1) - RXCN*rho*F_mat{i,j}(1) - RYCN *rho*F_mat{i,j}(2);
		
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
    

    end
	
	
	
k1 = k1+1;
end

    




%  drawnow
%  figure()
%  patch([X_1;X_2;X_4;X_3],[Y_1;Y_2;Y_4;Y_3],[Z_1;Z_2;Z_4;Z_3],[1,1,1]);hold on;
%  patch([X_1;X_2;X_4;X_3],-[Y_1;Y_2;Y_4;Y_3],[Z_1;Z_2;Z_4;Z_3],[1,1,1]);hold on;
%  quiver3(XC,YC,ZC,-FY,FX,FZ,'k');hold on;
%  quiver3(XC,-YC,ZC,-FY,-FX,FZ,'k');hold on;
%  xlabel('X');ylabel('Y');zlabel('Z');
 
 A = [NODE_LABEL;FX;FY;FZ;MX;MY;MZ];
 
fileID = fopen('INPUT.dat','w');
% fprintf(fileID,'%1$s %2$s %3$s %4$s %5$s %6$s %7$s \n','#Node_Label','Fx','Fy','Fz','Mx','My','Mz');
fprintf(fileID,'%4u %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f\n',A);
fclose(fileID);
end

