x=-1:0.01:1;
[c,h] =meshgrid(x);
%analytical wave equation solved. 
A = sin(pi.*c).*cos(pi.*h);
subplot(2,1 ,1);
imagesc(x,x,A)
colorbar
title('wave solution'); 
%% spherical waves 
A1 = abs(1./c).*sin(pi.*c).*cos(pi.*h);

subplot(2, 1, 2); 
imagesc(x,x,A1)
title('spherical waves'); 
colormap(jet); 
colorbar

