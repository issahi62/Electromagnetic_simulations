%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% CONSTANTS DEFINITION %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% SLM parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%
slmx=600; slmy=600;
rad=300;
ordvor=1;
rotoff=45;
rotang =45; 
offx=-55; offy=0;
gsc=230/2/pi;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% Beam controlling parameters %%%%%%%%%%%%%%%%%
BeamPhase=[40 40];
% Rotation of beams in a plane perpendicular to the propagation axis%
BeamAngle=[0 180];
BeamFocus=[0 0];
%%% Moves beams in a plane perpendicular to the propagation axis%%%%%
BeamDeflectRad=[.5 2];
%%% Ads topological charge to the beams %%%%
BeamVortex=[1 -1];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% VARIABLE DEFINITIONS %%%%%%%%%%%%%%%%%%%
BeamAngle=BeamAngle/180*pi;
rotoff=rotoff/180*pi;
x=repmat((1:slmx),[slmy,1])-slmx/2;
y=repmat((1:slmy)',[1,slmx])-slmy/2;
x=x/max(max(x))/2;
y=y/max(max(y))/2;
sq=x.^2+y.^2;
sq=sq*2*pi;
[theta, r]=cart2pol(x/rad,y/rad);
idx=(r<=1);
BaseMask=zeros(slmy,slmx);
BaseMask(idx)=1;
gkx=2*pi*x;
gky=2*pi*y;
z=x+1i*y;
z=z*exp(1i*(rotang+rotoff));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% WAVE FRONT CORRECTION (IF NEEDED) %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% An image file containing the correction should be provided %%%%%
% optcorr=zeros(slmy,slmx);
% rim=optcorr<=1;
% a=exist(folder,'file');
% if a==2
%     rim2=imread("CORRECTION FILE",'bmp')';
% if size(rim2,2)>=slmy && size(rim2,1)>=slmx
% ii=round(0.5*(size(rim2,1)-slmx))+1;
% ie=ii+slmx-1;
% ji=round(0.5*(size(rim2,2)-slmy))+1;
% je=ji+slmy-1;
% optcorr=(double(rim2(ii:ie,ji:je))*2*pi/255)';
% end
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% BEAMS GENERATION %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% We can generate one or two beams %%%%%%%%%%
mask=zeros(slmy,slmx)*(1+1i*1);
for i=1:1:size(BeamAngle,2)
z1=z*exp(1i*BeamAngle(1,i));
x1=real(z1);
y1=imag(z1);
AiryBeamam=(x1.^3+y1.^3)/3;
AiryBeam=AiryBeam/max(max(AiryBeam))*0.5*BeamPhase(1,i)*pi;
vortex=angle((x1+1i*y1).^(BeamVortex(1,i)));
focus=BeamFocus(1,i)*sq;
defl=BeamDeflectRad(1,i)*exp(1i*(BeamAngle(1,i)+rotang));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% GENERATES EITHER AIRY-VORTEX OR LAGUERRE-GAUSS BEAMS %%%%%%
%%%%%% (ONE AT THE TIME) %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Generates Airy-vortex beams %%%%%%
mask=mask+exp(1i*(AiryBeam+vortex+focus+real(defl)*...
gkx-imag(defl)*gky+4*pi*(rand(1,1)-0.5)));
%%%%%% Generates Laguerre-Gaussian beams %%%%%%%%
mask=mask+exp(1i*(vortex+focus+...
real(defl)*gkx-imag(defl)*gky+4*pi*(rand(1,1)-0.5)));
end
mask=mask.*exp(1i*(offx*gkx+offy*gky+optcorr));
tmask=mask/max(max(abs(mask)));
%%%%%%%% Mask displayed on the slm %%%%%%%%%%%%%
maskout=(angle(mask+(1-abs(mask).*exp(1i*(-34*gkx+45*gky))))+pi)*gsc;