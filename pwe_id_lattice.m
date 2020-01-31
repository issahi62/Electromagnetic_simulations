% INITIALISE 
clc
close all
clear 

%Igor photonic crystals.

%% DASHBOARD
t1 = .5; % thickness value
t2 = .5; 
e1 = 1; 
e2 = 1.52^2; % permittivity
a = t1+t2; % period 


numG = 50; % number of plane waves
N = numG;% due to zero symmetry. 

 % fourier coefficients of the lattice
 
 %% fourier coefficients
 counter_k = 0; % wave vector 
 counterG =1; 
 counterG1 = 1; 
 
 for G = -N*2*pi/a: 2*pi/a:N*2*pi/a
     for G1 = -N*2*pi/a: 2*pi/a:N*2*pi/a
         if (G-G1) == 0
             chi(counterG1, counterG) = 1/(a)*...
                 (1/e1*t1 + 1/e2*t2); 
         else 
             chi(counterG1, counterG) = 1i/(a)/(G-G1)*...
                 (1/e1*(exp(-1i*(G-G1)*t1)-1)+...
                  1/e2*(exp(-1i*(G-G1)*a)))-exp(-1i*(G-G1)*t1); 
         end
         counterG = counterG+1;
     end 
     counterG1 = counterG1+1;
     counterG =1; 
 end 

 
 counterG =1; 
 counterG1 = 1; 
 %% formation of Q matrix
 %k is bolch wave vectors
 for k = -pi/a : .2*pi/a : pi/a
     for G = -N*2*pi/a: 2*pi/a:N*2*pi/a
         for G1 = -N*2*pi/a: 2*pi/a:N*2*pi/a
             M(counterG1, counterG) = chi(counterG1, counterG)*(k+G1)*(k+G);
             counterG = counterG+1;
         end 
         counterG1 = counterG1+1; 
         counterG =1;
     end 
     counterG1 =1;
     
     V = eig(M); 
     counter_k = counter_k+1; 
     dispersion(counter_k, :) = sqrt(sort(abs(V))); 
     ks(counter_k)=k; 
 end 
 dispersion = dispersion*a/(2*pi);
 hold on;
 for u =1:8 
     plot(ks, dispersion(:, u), 'LineWidth', 2); 
 end 
 xlabel('wave vector'); 
 ylabel('a/lambda'); 
 xlim([-pi/a, pi/a]); 
 
