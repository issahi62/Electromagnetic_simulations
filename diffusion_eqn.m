
%% INITIALIZE 
clc
close all 
clear 

%% DASHBOARD
L =1; % length of wire
T = 1; % max time; 

%% Parameters 
maxk = 2500; 
dt = T/maxk; 

n = 50; % number of space steps 
nint = 50; % wavefront intermediate points 

dx = L/n; 
cond = 1/4; % conductivity 

b = 2.*cond*dt/(dx*dx); 

%% Initial temp of the wire a sinus
x =zeros(1,n); 
u = zeros(1,n); 
for bj = 1:9
for i = 1:n+1
    x(i) = (i-1)*dx; 
    u(i,1) = sin(bj*pi*x(i)); 
    
end


% temperature at the boundary (T=0)
time = zeros(1, maxk); 
for k =1:maxk+1
    u(1,k)=0;
    u(n+1,k) =0; 
    time(k) = (k-1)*dt; 
end

% implementation of the explicit method 
for k = 1:maxk
    for i = 2:n
        u(i, k+1) = u(i,k)+ 0.5*b*(u(i-1,k)+u(i+1,k)-2.*u(i,k));
    end
end

figure(1); 
subplot(3,3,bj)
plot(x,u(:,1), '-', x,u(:,100), '-', x, u(:,300), '-', x, u(:,600), '-'); 
title('Temperature values'); 
xlabel('X'); 
ylabel('T'); 

figure(2)
subplot(3, 3, bj)
meshc(x, time, u.' ); 
%shading interp
title('Temperature values'); 
xlabel('X'); 
ylabel('T');
colormap('jet'); 
%view(0,-90);
pause(.5);
end
