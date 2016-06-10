%% title:10 points in 3D using regularised Bio-Savart Law
clear all;
clc;
%% initialization
%the consequence of the vortex particles
n=10;
%time step length
dt=0.015;  t=100;
%3-D model in inertial frame, supposing the max time step is L*dt
core=1;        %constant
u=cell(n,t);   %3 vectors
w=cell(n,t);
r=cell(n,t);
%set the n particles, but the orders of u,w,r are not sure
for i=1:n
    r{i,1}=2*rand(1,3)-1;w{i,1}=2*rand(1,3)-1;u{i,1}=2*rand(1,3)-1;%randomly allocated particles
end

%% time-stepping calculation
%at the ith time step
for i=1:t
    %for the jth particle
    for j=1:n
        dw=0;du=0;
        %the kth(not itself) vortex particle's contribution
        for k=1:n
            if k~=j
                R=r{j,i}-r{k,i};
                %depending on the model,'exp' for example
                ro=norm(R/core);
                %different singularised models of f and g
%                 f=0;
%                 g=1;
                  f=sqrt(2/pi)*exp(-0.5*ro^2); 
                  g=erf(ro/sqrt(2))-ro*f;
                
%                  f=7.5/(ro^2+1)^3.5;
%                  g=(ro^3*(ro^2+2.5))/(ro^2+1)^2.5;
                
%                 f=3*exp(-ro^3); 
%                 g=1-f/3;
                du=du+g/(4*pi)*cross(R,w{k,i})/((norm(R))^3);
                dw=dw+dt*(1/(4*pi*core^3)*(-g*cross(w{j,i},w{k,i})/((norm(R/core))^3)+1/(norm(R))^2*(3*g/(norm(R/core)^3)-f)*(w{j,i}*R')*(cross(R,w{k,i}))));               
            end
        end
        %core variation is neglacted
        %refresh the u,r,w of jth particle,maybe Runge-Kutta method could be used
        u{j,i+1}=u{j,i}+du;
        r{j,i+1}=r{j,i}+u{j,i+1}*dt;
        w{j,i+1}=w{j,i}+dw;
    end
end
%% plot by each particle
% X=zeros(1,t);
% Y=zeros(1,t);
% Z=zeros(1,t);
% 
% for j=1:n
%     for i=1:t
%         X(i)=r{j,i}(1);
%         Y(i)=r{j,i}(2);
%         Z(i)=r{j,i}(3);
%     end
%     plot3(X,Y,Z);hold on;
% end
% title('trajectories of 10 points using Bio-Savart Law');
% grid on;
% xlabel('x');
% ylabel('y');
% zlabel('z');
%% plot by each timestep using writevideo
% X=zeros(1,n);
% Y=zeros(1,n);
% Z=zeros(1,n);
% %# preallocate
% nFrames = t;
% vidObj = VideoWriter('Vatistas.avi');
% vidObj.Quality = 100;
% vidObj.FrameRate = 10;
% open(vidObj);
% for i=1:t
%     for j=1:n
%         X(j)=r{j,i}(1);
%         Y(j)=r{j,i}(2);
%         Z(j)=r{j,i}(3);
%     end
%     scatter3(X,Y,Z,'r.');hold on;
%     axis tight;
%     writeVideo(vidObj, getframe(gca));
% end
% %movie(M);
% close(gcf);
% close(vidObj);
% winopen('Vatistas.avi');


%% plot by each timestep using getframe
X=zeros(1,n);
Y=zeros(1,n);
Z=zeros(1,n);
nFrames = t;
mov(1:nFrames) = struct('cdata',[], 'colormap',[]);
for i=1:t    
    for j=1:n
        X(j)=r{j,i}(1);
        Y(j)=r{j,i}(2);
        Z(j)=r{j,i}(3);
    end
    scatter3(X,Y,Z,'r.');hold on;
    if i==1
        title('trajectories of 10 points using Bio-Savart Law');
        grid on;
        xlabel('x');
        ylabel('y');
        zlabel('z');
    end
    mov(i)=getframe(gca);
end
% movie(mov);
% movie2avi(mov, 'BS Laws.avi', 'compression','None', 'fps',10);
% winopen('BS Laws.avi')


