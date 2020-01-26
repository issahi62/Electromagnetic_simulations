
% Setup Grid
xmax = 250.00; % x grid size in microns
ymax = xmax; % y grid size in microns
N = 2*1024+1; % number of x & y grid points
v = 0:1:(N-1); % integer vector from 0 to N-1
dx = xmax/(N-1); % spacing in x
dy = ymax/(N-1); % spacing in y
[x,y] = meshgrid(v,v); % meshgrid containing NxN points (integer values)
x = -xmax/2 + x*dx; % shifts mesh to cover +/- xmax/2, 2D-grid
y = -ymax/2 + y*dy; % shifts mesh to cover +/- ymax/2, 2D-grid
xar = -xmax/2 + v*dx; % cross section to cover +/- xmax/2, 1D
yar = -ymax/2 + v*dy; % cross section to cover +/- ymax/2, 1D
dkx = 2*pi/xmax; % spacing in kx (inverse microns)
dky = 2*pi/ymax; % spacing in ky (inverse microns)
kxmax = 2*pi/dx; % grid size in microns (from kx1 to kxN)
kymax = 2*pi/dy; % grid size in microns (from ky1 to kyN)
[kx,ky] = meshgrid(v,v); % meshgrid containing NxN points (integer values)
kx = -kxmax/2 + kx*dkx; % shifts mesh to cover +/- kxmax/2, 2D-grid
ky = -kymax/2 + ky*dky; % shifts mesh to cover +/- kymax/2, 2D-grid
kxar = -kxmax/2 + v*dkx; % cross section to cover +/- kxmax/2, 1D
kyar = -kymax/2 + v*dky; % cross section to cover +/- kymax/2, 1D
% System Properties
wx = 15; % field width in x in microns
wy = wx; % field width in y in microns
E0 = 1; % peak field strength (arbitrary units)
lambda = 1; % wavelength in microns
L = 1.00e4; % distance to observation screen in microns

l = 6; % azimuthal mode index
% convert to polar coordinates
[theta,rho]=cart2pol(x,y);
% Generate aperture
E = zeros(N); % 0 matrix to setup aperture & propagate
aperture = 5;
switch aperture
case 1
% Triangular aperture
deg = 30*pi/180; % angle of equilateral triangle in radians
h = find (y>(x./(tan(deg))-(wx)*2/3)&y>-x./tan(deg)-(wx)*2/3&y<((wx)*1/3));
case 2
% Square Aperture
Squareaperture = (abs(x)>=0 &abs(x)<wx).*(abs(y)>0 &abs(y)<=wy);
h = find (Squareaperture);
case 3
% Circular Aperture
R = 20;
Circularaperture = sqrt(abs(x).^2 + abs(y).^2) < R ;
h = find(Circularaperture);
case 4
% Annular Triangle with Radial spokes Aperture
wxInner=wx;
wx1=wx;
dw1=2;
deg = 30*pi/180; % angle of equilateral triangle in radians
h1= find (y>(x/tan(deg)-(wx1+dw1)*2/3)&y>-x/tan(deg)&(wx1+dw1)*2/3&y<((wx1+dw1)*1/3)&~(y>(x/tan(deg)-(wx1-dw1)*2/3)&y>-x/tan(deg)-(wx1-dw1)*2/3&y<((wx1-dw1)*1/3)));
% Create cross
lineLength=10; %microns
lineThickness=0.1*lineLength;
apAngle=pi/6;
h2=find((rho.*cos(theta-apAngle)>=0 & rho.*cos(theta-apAngle)<=lineLength &abs(rho.*sin(theta-apAngle))<=lineThickness/2) |(rho.*cos(theta+2*pi/3-apAngle)>=0 & rho.*cos(theta+2*pi/3-apAngle)<=lineLength & abs(rho.*sin(theta+2*pi/3-apAngle))<=lineThickness/2)|(rho.*cos(theta-2*pi/3-apAngle)>=0 & rho.*cos(theta-2*pi/3-apAngle)<=lineLength & abs(rho.*sin(theta-2*pi/3-apAngle))<=lineThickness/2));
h=union(h1,h2);
case 5
% Multi Annular Triangle Aperture
% First Annular Triangle Aperture
wx1=wx;
dw1=2;
deg = 30*pi/180; % angle of equilateral triangle in radians
h1 = find (y>(x/tan(deg)-(wx1+dw1)*2/3)&y>-x/tan(deg)-...
(wx1+dw1)*2/3&y<((wx1+dw1)*1/3)&...
~(y>(x/tan(deg)-(wx1-dw1)*2/3)&y>-x/tan(deg)-(wx1-dw1)*2/3&y<((wx1-...
dw1)*1/3)));
% Second Annular Triangle Aperture
wx2=2*wx1;
dw1=2;
h2 = find (y>(x/tan(deg)-(wx2+dw1)*2/3)&y>-x/tan(deg)-...
(wx2+dw1)*2/3&y<((wx2+dw1)*1/3)&...
~(y>(x/tan(deg)-(wx2-dw1)*2/3)&y>-x/tan(deg)-(wx2-dw1)*2/3&y<((wx2-...
dw1)*1/3)));
h = union(h1,h2);
end
% find martix elements inside used aperture
E(h) = 1; % set aperture values to 1
Eaperture = E; % store aperture
% Generate Fraunhofer propagation pattern
w0 = 7;
p = 5; % radial mode index
T = theta(h);
R = rho(h);
Z=[];
beamWaist=w0; % beam waist
pValue = p; % radial index value
lValue = l; % azimuthal index value
