x = linspace(0,2*pi,100);
y = linspace(-pi,pi,100);
r = sqrt(x.^2 + y.^2);
t = linspace(-1,1,100);
k = 2*pi;
w_0 = 1;
Y = exp(1i*k.*r - w_0.*t);
subplot(2,1,1)
imagesc(real(Y))
colormap(jet)
subplot(2,1,2)
imagesc(abs(Y))
colormap(jet)



[x,y] = meshgrid(-pi:.01:pi);
 
Y = test3(x,y);

subplot(2,2,1)
imagesc(real(Y))
axis equal tight
colormap(jet)
subplot(2,2,2)
imagesc(abs(Y))
axis equal tight
colormap(jet)

