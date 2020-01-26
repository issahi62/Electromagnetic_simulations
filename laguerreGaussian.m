% fieldValues=laguerreGaussian(R,T,Z,pValue,lValue,lambda,beamWaist)
%
% Input:
% R: radial coordinate grid [m]
% T: azimutal coordinate grid [rad]
% Z: axial coordinate grid [m]
% pValue: the p value of the beam (radial)
% lValue: the l value of the beam (azimutal phase index, must be integer)
% lambda: the wavelength to simulate [m]
% beamWaist: the beam waist [m]
%
% Example:

%
function fieldValues=laguerreGaussian(R,T,Z,pValue,lValue,lambda,beamWaist)
 if (nargin<1)
 xRange = [-768/2:1:768/2];
 [X,Y] = meshgrid(xRange,xRange);
 [T,R] = cart2pol(X,Y);
 Z = [];
 pValue = 4;
 lValue = 3;
 end
 if (nargin<6)
 lambda = 500e-9;
end
 if (nargin<7)
     beamWaist = 100;
 end
 if (isempty(Z))
     
    Z = zeros(size(R));
 end

 rayleighRange = pi*beamWaist^2/lambda;
 W = beamWaist*sqrt(1+(Z./rayleighRange).^2);
 radiusofCurvature = (Z+(Z==0)).*(1+(rayleighRange./Z).^2);
 gouyPhase = atan2(Z,rayleighRange);
 fieldValues = (1./W).*(R.*sqrt(2)./W).^abs(lValue).*exp(-R.^2./W.^2)...
.*L(abs(lValue),pValue,2*R.^2./W.^2).*exp(1i.*(2*pi/lambda).*R.^2./(2*radiusofCurvature))...
.*exp(1i*lValue.*T).*exp(-1i*(2*pValue+abs(lValue)+1).*gouyPhase);
 %Normalize
 fieldValues=fieldValues./sqrt(sum(abs(fieldValues(:)).^2));
 if (nargout==0)
 figure;
 colormap gray;
 imagesc(abs(fieldValues).^2);
 %clear fieldValues;
 end
end