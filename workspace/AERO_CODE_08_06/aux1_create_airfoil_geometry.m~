%% READ INPUT DATA FOR GIVEN AIRFOIL GEOMETRY 

close all;clc;clear all;

%% Adding paths

addpath([pwd,'/geometry']);
addpath([pwd,'/functions']);
addpath([pwd,'/data']);

%% Load geometry data ( Camber Line and Mesh )

% Number of points 

M = 100; 

filename = [pwd,'/geometry/airfoil_geometry_flatplate.dat'];
[yc,xx] = func_import_airfoil_geometry(filename,M);

%% To compute curvature of camber line

theta=acos(1-(2*xx));
NN=length(yc);


%% Calcoli per la derivazione numerica attraverso la funzione Pesider

Dc = zeros(NN,NN);

for i=2:NN-1
    xxst =xx(i-1:i+1);           % Stencil per la derivata prima.
    xxc = xx(i);                 % Punto di collocazione per la der. prima.
    wddc = func_der(xxst,xxc,1);    % Pesi della derivata prima.
    Dc(i,i-1:i+1) = wddc;       % Riempimento della matrice.
end


Dc(1,1:2)=func_der(xx(1:2),xx(1),1);
Dc(NN,NN-1:NN)=func_der(xx(NN-1:NN),xx(NN),1);  % Calcolo dei coefficienti ai limiti del dominio 

yyc=Dc*yc';  % Calcolo della derivata della linea media
plot(theta,yyc);axis equal;
xlabel('\theta');
ylabel('d(z_{c})/{d\theta}');


%% To save curvature in .mat file

save([pwd,'/data/data_airfoil.mat'],'yc','yyc','xx','-append');
