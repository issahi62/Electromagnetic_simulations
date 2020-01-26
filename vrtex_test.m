fieldValues(h)=laguerreGaussian(R,T,Z,pValue,lValue,lambda,beamWaist);
 % calculate the generalized LG beams
E(h) = fieldValues(h); % electric field envelope
k0 = 2*pi/lambda; % propagation vector prefactor
X = kx*lambda*L/2/pi; % propagated X grid, 2D
Y = ky*lambda*L/2/pi; % propagated Y grid, 2D
Xar = kxar*lambda*L/2/pi; % cross section of X Grid, 1D
dX = dkx*lambda*L/2/pi; % propagated dX spacing
dY = dky*lambda*L/2/pi; % propagated dY spacing
Ek1 = dx*dy*fftshift(fft2(E)); % FT of pattern
EX1 = (1/sqrt(-1)/L/lambda)*exp(sqrt(-1)*k0*L)...
 *exp(sqrt(-1)*k0*(X.^2+Y.^2)/2/L).*Ek1;
 % Full Fraunhofer propagation
% Plot
% scrsz = get(0,'ScreenSize'); % gets screen size
% figure('OuterPosition',[1 .05*scrsz(4) scrsz(3)/1 .95*scrsz(4)])
figure;
subplot(1,2,1);
imagesc(abs(E));
axis square;
title('amp');
subplot(1,2,2);
imagesc(angle(E));
title('phase');
axis square;
imagesc(Xar,Xar,abs(EX1).^2);
axis equal
set(gca,'FontSize',15);
set(gca,'color',[0 0 0]);
title(['m = ', num2str(m)])
xlim([-1.5*L 1.5*L]);
ylim([-1.5*L 1.5*L]);



