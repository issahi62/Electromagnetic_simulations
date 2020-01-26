[c,h] =meshgrid(-1:0.01:1);
for t = 1:100
A = exp((h.^2+c.^2)./2); 
imagesc(A);
colormap('jet');
drawnow; 
pause(0.1);
end 
imagesc(abs(fftshift(fftn(A))))