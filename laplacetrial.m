clear all hidden
close all 
clc

lx = 10;  % size of the potential 
ly =10; 
Nx =60; % number of grid size 
Ny =60; 

dx = lx/(Nx-1); 
dy = ly/(Ny-1); 

xa = 0:dx:lx; 
ya = 0:dy:ly; 

numiter = 1000; % number of iterations 
% initicial conditions 
p = zeros(Nx, Ny);
pn = zeros(Nx, Ny); 


row = 2:Nx-1; % 1 and Nx has a boundary condition of either 10 or 0; 
col = 2:Ny-1; 
for iter = 1:numiter
    pn=p; 
    p(row, col) = ((dy^2)*(pn(row+1, col)+...
                  pn(row-1, col))+...
                  (dx^2)*(pn(row, col+1)+...
                  pn(row, col-1)))...
                  /(2*(dx^2+dy^2)); 
              
              %% boundary conditions 
              
    p(:, 1) =10; % right and left side of the box
    p(:, Nx) =10; 
    p(1, :) =0; % up and down of the box
    p(Ny, :) = 0; 
end 

%% plotting 
surf(xa, ya, p);
shading interp
colormap('jet'); 
title('laplace'); 
%imagesc(p)



