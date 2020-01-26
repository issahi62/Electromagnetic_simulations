%% DEVELOPED BY Issah 
% replace the lattice_name (sq5224nm) with  your lattice name
%% INITIALIZE 
clc
close all 
clear hidden 

%% FILES FOR SIMULATION 
s = get(0, 'ScreenSize');
figure1 = figure('Color', 'white', 'Position', [0 0 s(3) s(4)]);
%load sq5224nm.mat
Snorm = sq5224nm; 


set(0, 'defaultaxesfontname', 'Helvetica');
set(0, 'defaultaxesfontsize', 15);
set(0, 'defaulttextfontname', 'Helvetica');
set(0, 'defaulttextfontsize', 15);
set(0, 'defaultlinelinewidth', 2);
% Codes for image processing

%% DASHBOARD 
N = 1000; % size to be taken
r = 8; % radius
SY =0;% shift in either the x or y direction
SX =0;
hev = 4.135e-15;% plancks constant
c = 2.99e8;% speed of light
n =1.52; % refractive index;
wav = Snorm(1:N, 1).*1e-9; % wavelength
wav = hev*c./wav;


%% PROCESSSING USING IMAGEADJUST
S = Snorm(1:N, 2:end);
I = 1:size(S,1);
I = linspace(-2, 2, size(I,2)); 
KK = norm(S); 
S = S./KK;
S = imadjust(S);
S(S<-1)=-1;
S(isnan(S))=0;
S(S>1)=1;

%% signal noise removal
fSn = fftshift(fftn(S)); 
SfS=size(S);
FI = SfS; 
[X, Y] = meshgrid(1 : SfS(2), 1 : SfS(1));
Z =sqrt((X-SfS(2)/2).^2+(Y-SfS(1)/2).^2); % circle
CS = fSn;
CS(Z>r)=0;
[mxV,indMax]=max(CS(:));
[ny,nx]=ind2sub(SfS, indMax);
Rt=sqrt((X-nx).^2+(Y-ny).^2);
%%
%set values outside side maximum zero
Ft=CS;
Ft(Rt>r)=0;
Fc=circshift(Ft,[round(FI(1)/2 - ny)+SY, round(FI(2)/2 - nx)+SX]);
FFtf = ifftshift(Fc); 
% Actual image (Complete)
InF = (ifft2(FFtf));% inverfourier
S2 = (InF)';

%% VISULIZE IMAGE
subplot(1,3,1)
imagesc((I), wav, (medfilt2(imadjust(S')))); 
shading flat
colorbar
xlabel('$k_{||} (\mu^{-1} m)$', 'Interpreter', 'latex'); 
ylabel('E(eV)', 'Interpreter', 'latex'); 
axis xy;
title('Experimental'); 
caxis([0,1]);
subplot(1,3,2); 

imagesc(I, wav, imsharpen(histeq(abs(S2).^2)));
shading flat
colorbar
xlabel('$k_{||} (\mu^{-1} {m})$', 'Interpreter', 'latex'); 
ylabel('E(eV)', 'Interpreter', 'latex'); 
colormap(brewermap([],'YlGnBu'))
title('Simulated'); 
%caxis([0 1]);
axis xy;
subplot(1,3,3)
imagesc(angle(S2).^2);
shading flat
colorbar
colormap('parula')
title('phase'); 
axis xy;


%% creating Lines on image 
annotation(figure1,'line',[0.41796875 0.578125],...
    [0.758695652173913 0.165217391304348],'Color',[1 0 0],'LineWidth',4,...
    'LineStyle',':');
annotation(figure1,'line',[0.5796875 0.4203125],...
    [0.710869565217391 0.143478260869565],'Color',[1 0 0],'LineWidth',4,...
    'LineStyle',':');
% Create line
annotation(figure1,'line',[0.1328125 0.30078125],...
    [0.756521739130435 0.152173913043478],'Color',[1 0 0],'LineWidth',4,...
    'LineStyle',':');

% Create line
annotation(figure1,'line',[0.29765625 0.140625],...
    [0.692478260869565 0.147826086956522],'Color',[1 0 0],'LineWidth',4,...
    'LineStyle',':');


