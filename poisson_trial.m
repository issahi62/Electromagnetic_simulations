clear all hidden
close all 
clc

lx = 10;  % size of the potential 
ly =10; 
Nx =100; % number of grid size 
Ny =100; 

dx = lx/(Nx-1); 
dy = ly/(Ny-1); 

xa = 0:dx:lx; 
ya = 0:dy:ly; 

%% source term variable 
f = zeros(Nx, Ny); 
f (round(Nx/4), round(Ny/4)) = 3000; 
f (round(3*Nx/4), round(3*Ny/4)) = -3000; 
numiter = 1000; % number of iterations 
% initicial conditions 
p = zeros(Nx, Ny);
pn = zeros(Nx, Ny); 


%% boundary conditions

p(:, 1) =0; % right and left side of the box
p(:, Nx) =0;
p(1, :) =0; % up and down of the box
p(Ny, :) = 0;
    

row = 2:Nx-1; % 1 and Nx has a boundary condition of either 10 or 0; 
col = 2:Ny-1; 
for iter = 1:numiter
    pn=p; 
    p(row, col) = ((dy^2)*(pn(row+1, col)+...
                  pn(row-1, col))+...
                  (dx^2)*(pn(row, col+1)+...
                  pn(row, col-1))-...
                  (f(row, col)*dx^2*dy^2))/(2*(dx^2+dy^2)); 
   
end
%% plotting 
pcolor(xa, ya, p);
shading interp
colormap('jet'); 
title('laplace'); 
%imagesc(p)



