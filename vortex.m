

AU=2.1e-6; % pixel size
W=1e-3; % size of sample
n=round(W/AU); % samples
l=2; % charge
gamma=0; % arbritrary phase in mode converter
delta=1; % Mode conversion
a=1; % Mode converter deformation
CH=1; % +/- change in phi
N=1.468; % index of plate
mask=zeros(n); % make proper sized zeros called mask
maskt=zeros(n); % maskt is the size of the mode converter
mask2=zeros(n); % Output mask
dlam=1; % change in wavelngth
lam=dlam*532e-9; % Reconstruction wavelength
f=3e8/lam; % Frequency
k=2*pi/lam; % wave number
w=2*pi*f; % angular frequency
t=0; % Time variable
z=0; % direction of propagation
I=1:n; % Sizing
x=I-n/2; % x vector
y=n/2-I; % y vector
[X,Y]=meshgrid(x,y); % vector of x,y
phi=atan2(X,Y); % vortex phase
r=sqrt(x.^2+y.^2); % not needed

A=exp(-1i.*l.*phi); % vortex input
A1=exp(-1i.*(w.*t-k.*z)).*A; %input wave
t=exp(-1i.*a.*((delta.*(CH.*phi))-gamma)); % Mode converter
% t=phi.*(N-1).*l.*lam./(2*pi); % SPP
t2=(delta.*(CH.*phi)).*a-gamma-k.*X; % ignore
A2=t.*A1; % output
mask=(1/(2*pi)).*mod(angle(A),2*pi); % phase profile of vortex
maskt=(1/(2*pi)).*exp(-1i.*mod((t2),2*pi)); % profile of mode converter
mask2=(1/(2*pi)).*mod(angle(A2),2*pi); % angular profile of output beam
% Plot input wave
figure
imagesc(((mask)))
colormap(jet)
colorbar