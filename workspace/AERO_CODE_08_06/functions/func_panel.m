function vortex_lattice = func_panel(M,N,index,patches,vortex_ring)


%% Pre-Allocation for speed

vortex_lattice = cell(M,N);
X1 = zeros(M,N);
X2 = zeros(M,N);
X3 = zeros(M,N);
X4 = zeros(M,N);
Y1 = zeros(M,N);
Y2 = zeros(M,N);
Y3 = zeros(M,N);
Y4 = zeros(M,N);
Z1 = zeros(M,N);
Z2 = zeros(M,N);
Z3 = zeros(M,N);
Z4 = zeros(M,N);
XC = zeros(M,N);
YC = zeros(M,N);
ZC = zeros(M,N);
Nx = zeros(M,N);
Ny = zeros(M,N);
Nz = zeros(M,N);

%% Chord-wise loop and Span-wise loop to form the lattice

for i = 1 : M
    
    % Span-wise loop
    
    for j = 1 : N
        
        % X-COORDINATE CORNER POINTS
        
        vortex_ring.X_1 = patches{1,index}.X_L - ...
            (patches{1,index}.X_L - patches{1,index}.X_R)*(j-1)/N + ...
            1/4*(patches{1,index}.c_L -(patches{1,index}.c_L - patches{1,index}.c_R)*(j-1)/N)/M+...
            (patches{1,index}.c_L -(patches{1,index}.c_L - patches{1,index}.c_R)*(j-1)/N)*(i-1)/M;
        X1(i,j) = vortex_ring.X_1; % CHECK
        
        vortex_ring.X_2 = patches{1,index}.X_R + ...
            (patches{1,index}.X_L - patches{1,index}.X_R)*(N-j)/N + ...
            1/4*(patches{1,index}.c_R +(patches{1,index}.c_L - patches{1,index}.c_R)*(N-j)/N)/M+...
            (patches{1,index}.c_R +(patches{1,index}.c_L - patches{1,index}.c_R)*(N-j)/N)*(i-1)/M;        
        X2(i,j) = vortex_ring.X_2; % CHECK 
        
        vortex_ring.X_3 = vortex_ring.X_1 + ...
            (patches{1,index}.c_L -(patches{1,index}.c_L - patches{1,index}.c_R)*(j-1)/N)/M;
        X3(i,j) = vortex_ring.X_3; % CHECK
        
        vortex_ring.X_4 = vortex_ring.X_2 + ...
            (patches{1,index}.c_R +(patches{1,index}.c_L - patches{1,index}.c_R)*(N-j)/N)/M;
        X4(i,j) = vortex_ring.X_4; % CHECK
       
        
        % Y-COORDINATE CORNER POINTS
        
        vortex_ring.Y_1 = patches{1,index}.Y_L - ...
            (patches{1,index}.Y_L - patches{1,index}.Y_R)*(j-1)/N;
        Y1(i,j) = vortex_ring.Y_1; % CHECK
        
        vortex_ring.Y_2 = patches{1,index}.Y_R + ...
            (patches{1,index}.Y_L - patches{1,index}.Y_R)*(N-j)/N; 
        Y2(i,j) = vortex_ring.Y_2; % CHECK 
        
        vortex_ring.Y_3 = vortex_ring.Y_1;
        Y3(i,j) = vortex_ring.Y_3; % CHECK
        
        vortex_ring.Y_4 = vortex_ring.Y_2;
        Y4(i,j) = vortex_ring.Y_4; % CHECK  
        
         % Z-COORDINATE CORNER POINTS
        
        vortex_ring.Z_1 = patches{1,index}.Z_L - ...
            (patches{1,index}.Z_L - patches{1,index}.Z_R)*(j-1)/N;
        Z1(i,j) = vortex_ring.Z_1; % CHECK
        
        vortex_ring.Z_2 = patches{1,index}.Z_R + ...
            (patches{1,index}.Z_L - patches{1,index}.Z_R)*(N-j)/N; 
        Z2(i,j) = vortex_ring.Z_2; % CHECK 
        
        vortex_ring.Z_3 = vortex_ring.Z_1;
        Z3(i,j) = vortex_ring.Z_3; % CHECK
        
        vortex_ring.Z_4 = vortex_ring.Z_2;
        Z4(i,j) = vortex_ring.Z_4; % CHECK  
                
        % X-COORDINATE COLLOCATION POINT
        
        vortex_ring.X_C = ((vortex_ring.X_3+vortex_ring.X_1)/2 + ...
            (vortex_ring.X_4+vortex_ring.X_2)/2)/2;
        XC(i,j) = vortex_ring.X_C; % CHECK
        
        % Y-COORDINATE COLLOCATION POINT
        
        vortex_ring.Y_C = (vortex_ring.Y_1+vortex_ring.Y_2)/2;
        YC(i,j) = vortex_ring.Y_C; % CHECK
        
        % Y-COORDINATE COLLOCATION POINT
        
        vortex_ring.Z_C = (vortex_ring.Z_1+vortex_ring.Z_2)/2;
        ZC(i,j) = vortex_ring.Z_C; % CHECK
        
        % NORMAL VECTOR ( CROSS PRODUCT OF DIAGONAL VECTORS )
        
        A1 = vortex_ring.X_3 - vortex_ring.X_2;
        A2 = vortex_ring.Y_3 - vortex_ring.Y_2;
        A3 = vortex_ring.Z_3 - vortex_ring.Z_2;
        B1 = vortex_ring.X_4 - vortex_ring.X_1;
        B2 = vortex_ring.Y_4 - vortex_ring.Y_1;
        B3 = vortex_ring.Z_4 - vortex_ring.Z_1;
        XX = A2*B3 - A3*B2;
        YY = B1*A3 - A1*B3;
        ZZ = A1*B2 - A2*B1;
        A = sqrt(XX^2+YY^2+ZZ^2);
        
        vortex_ring.N_x = XX/A;
        vortex_ring.N_y = YY/A;
        
        
        switch sign(ZZ/A)
            case 1             
        vortex_ring.N_z = ZZ/A;
        vortex_ring.Z_SIGN = sign(ZZ/A);
            case -1             
        vortex_ring.N_z = -ZZ/A;
        vortex_ring.Z_SIGN = sign(ZZ/A);
            case 0
        vortex_ring.N_z = ZZ/A;
        vortex_ring.Z_SIGN = sign(ZZ/A);
        end
        
        Nx(i,j) = vortex_ring.N_x;
        Ny(i,j) = vortex_ring.N_y;
        Nz(i,j) = vortex_ring.N_z;
        
        
        % PANEL AREA
        
        C1 = vortex_ring.X_3 - vortex_ring.X_1;
        C2 = vortex_ring.Y_3 - vortex_ring.Y_1;
        C3 = vortex_ring.Z_3 - vortex_ring.Z_1;
        D1 = vortex_ring.X_2 - vortex_ring.X_1;
        D2 = vortex_ring.Y_2 - vortex_ring.Y_1;
        D3 = vortex_ring.Z_2 - vortex_ring.Z_1;
        E11 = C2*D3 - C3*D2;
        E12 = D1*C3 - C1*D3;
        E13 = C1*D2 - C2*D1;
        E21 = B2*C3 - B3*C2;
        E22 = C1*B3 - B1*C3;
        E23 = B1*C2 - B2*C1;
        
        
        vortex_ring.S = 0.5*(sqrt(E11^2+E12^2+E13^2)+sqrt(E21^2+E22^2+E23^2));
        
        % TWIST ANGLE
        
        vortex_ring.TWIST = -(interp1([patches{1,index}.Y_L,patches{1,index}.Y_R],[patches{1,index}.twist_L,patches{1,index}.twist_R],...
         vortex_ring.Y_C,'linear','extrap'));
        
        
        % COLLECT PANELS INFORMATIONS INTO CELL
        
        vortex_lattice{i,j} = vortex_ring;
    end
    
end