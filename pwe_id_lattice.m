 %%%%TE and TM modes for 1D lattice with field distribution %%%%%%
 %%
% INITIALISE 
clc
close all
clear 

%Igor photonic crystals.
TM = 0;
TE = 1; 
%% DASHBOARD
t1 = .5; % thickness value
t2 = .5; 
e1 = 1.52^2; 
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
     for G1 = -N*2*pi: 2*pi/a:N*2*pi/a
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
 for k = -pi/a : 1*pi/a : pi/a
     for G = -N*2*pi/a: 2*pi/a:N*2*pi/a
         for G1 = -N*2*pi/a: 2*pi/a:N*2*pi/a
             M(counterG1, counterG) = chi(counterG1, counterG)*(k+G1)*(k+G);
             Mtm(counterG1, counterG) = chi(counterG1, counterG)*(k+G).^2;
              counterG = counterG+1;
         end 
         counterG1 = counterG1+1; 
         counterG =1;
     end 
     counterG1 =1;
     
    [D, V] = eig(M); 
     [D2, V2] = eig(Mtm); 
     counter_k = counter_k+1; 
     dispersion(counter_k, :) = sqrt(sort(abs(diag(V)))); 
     
     dispersion2(counter_k, :) = sqrt(sort(abs(diag(V2))));
     eigvectors_pan(:, :, counter_k) = D; 
     eigvectors_pansTE(:, :, counter_k) = D2; 
     ks(counter_k)=k; 
 end 
 
 %%
 dispersion = dispersion*a/(2*pi);
 hold on;
 if TM == 1
 for u = 1:8
 plot(ks, dispersion(:, u), 'LineWidth', 2); 
 end
 end 
 if TE == 1
 for u =1:8 
     %plot(ks, dispersion(:, u), 'LineWidth', 2);
     plot(ks, dispersion2(:, u), 'LineWidth', 2); 
 end 
 end 
 xlabel('wave vector'); 
 ylabel('a/lambda'); 
 xlim([-pi/a, pi/a]); 
 
 if TM == 1 && TE==1
     disp('error'); 
     
 end 
 x = linspace(-1, 1); 
 y = linspace(-1,1); 
 if TM == 1
 figure(3);
 for j = 1:9
 PSI = eigvectors_pan(:, j, 1);
 PSI(end) =[]; 
 PSI = reshape(PSI, [10, 10]); 
 PSI = invFFT2D(PSI,100,100);
 E(:, :, j) = PSI / max(PSI(:)); 
 hold on
 subplot(3, 3, j)
 pcolor(x,y,abs(E(:, :, j)));
 shading flat
 %colormap(jet)
 colorbar
 xlabel('x'); 
 ylabel('y'); 
 end
 end
 
 if TE == 1
 figure(2);
 for j = 1:9
 PSI = eigvectors_pansTE(:, j,  1);
 PSI(end) =[]; 
 PSI = reshape(PSI, [10, 10]); 
 PSI = invFFT2D(PSI,100,100);
 E(:, :, j) = PSI / max(PSI(:)); 
 hold on
 subplot(3, 3, j)
 pcolor(x,y, real(E(:, :, j)));
 shading flat
 colorbar
 %colormap(jet)
  xlabel('x'); 
 ylabel('y'); 
 end
 end
 %colormap('jet'); 
 
 
 
 function [Vxy] = invFFT2D(Vk2D,Ny,Nx)

Nkx=length(Vk2D(1,:));
Nky=length(Vk2D(:,1));

Nx1=Nx/2-floor(Nkx/2);
Nx2=Nx/2+ceil(Nkx/2);
Ny1=Ny/2-floor(Nky/2);
Ny2=Ny/2+ceil(Nky/2);

Vk2D00=zeros(Ny,Nx);
Vk2D00( Ny1+1:Ny2 , Nx1+1:Nx2)=Vk2D;
Vxy=ifft2(ifftshift(Vk2D00));

end
