format compact
close all
clear all
clc
%% Reading in simulated object
I0=double(imread('F-16_250px.jpg')); % 400x400 pixels
% I0=double(rgb2gray(I0));
%% Parameter Setup
% length in cm
M=512; % # pixels, 1D
dx=10e-4; % pixel pitch 
w=633e-7; % wavelength 
z=5; % propagation distance
theta=4.1; % reference beam angle; degree
% res=w*z/M/dx % sampling distance (smallest resolvable element)
%% Simulating Target Phase
[x,y]=meshgrid(1:length(I0),1:length(I0));
PHASE=15*sin(((x-200).^2+(y-225).^2)/10000)+0.002*((x-37).^2+(y-100).^2);
I0=I0.*exp(1i*PHASE);
%% Noise
snr=.01; % White gaussian noise
wellDepth=40000;
G_AD=wellDepth/M % e-/px (Gain)
F_Ill=1 % 400e-/px
Ref_Ill=80 % 32000e-/px
%% Showing CCD FOV
I_FOV=zeros(M);
I_FOV((M/2-length(I0)/2+1):(M/2+length(I0)/2),(M/2-length(I0)/2+1):(M/2+length(I0)/2))=I0;
figure();
imshow(mat2gray(abs(I_FOV))); title('CCD FOV (Aimed at Target)')
axis off
%% Object field propagation to CCD
pad1=5; % Padding
r=1:pad1*M;
c=1:pad1*M;
[C, R]=meshgrid(c, r);
I=zeros(pad1*M);
I((M*pad1/2-length(I_FOV)/2+1):(M*pad1/2+length(I_FOV)/2),(M*pad1/2-length(I_FOV)/2+1):(M*pad1/2+length(I_FOV)/2))=I_FOV;
A0=fftshift(ifft2(fftshift(I)));
deltaf=1/pad1/M/dx;
p=exp(-2i*pi*z.*((1/w)^2-((R-M*pad1/2-1).*deltaf).^2-((C-M*pad1/2-1).*deltaf).^2).^0.5); % eq (4-20), Goodman
Az=A0.*p;
EO=fftshift(fft2(fftshift(Az))); % OBJ field at CCD (real domain)
EO=EO((M*pad1/2-M/2+1):(M*pad1/2+M/2),(M*pad1/2-M/2+1):(M*pad1/2+M/2)); % reduce diffraction-plane size
%% Interference at the hologram plane
% zero-padding in the spectrum domain
pad2=4;
r2=1:pad2*M;
c2=1:pad2*M;
Az=fftshift(ifft2(fftshift(EO)));
Az2=zeros(pad2*M);
Az2((M*pad2/2-M/2+1):(M*pad2/2+M/2),(M*pad2/2-M/2+1):(M*pad2/2+M/2))=Az;
EOf=fftshift(fft2(fftshift(Az2))); % OBJ field at CCD (real domain)
AV=(min(min(abs(EOf)))+max(max(abs(EOf))))/2; % ref wave amplitude
EOf=(F_Ill/AV)*EOf; % Scale amplitude to 400 e-/px
[C2, R2]=meshgrid(c2, r2);
Ref=Ref_Ill*exp(1i*2*pi*sind(theta)*dx/4.*(R2-M*pad2/2-1)/w+1i*2*pi*sind(theta)*dx/4.*(C2-M*pad2/2-1)/w); % eq 3.5a - DH textbook, with scaled amplitude set to 32000 e-/px
%IH=(EOf+Ref).*conj(EOf+Ref)/G_AD+awgn(real(Ref),snr,'measured')/G_AD; % |F+R|^2 + noise
IH=(EOf+Ref).*conj(EOf+Ref)/G_AD+awgn(abs(Ref),snr)/G_AD; % |F+R|^2 + noise
scale=.5; % pad3/pad2
IH=IH((M*pad2/2-M*scale*pad2/2+1):(M*pad2/2+M*scale*pad2/2),(M*pad2/2-M*scale*pad2/2+1):(M*pad2/2+M*scale*pad2/2));
figure; imshow(mat2gray(IH)); title('Hologram')
SP=fftshift(fft2(fftshift(IH)));
figure; imshow(500.*mat2gray(abs(SP)));title('Hologram spectrum')
%% Windowing
hold on
Window1=150;%200;%180;%150;
Window2=321;%247;%271;%321;
WindowPixels=(Window2-Window1)^2;
W=Window1:Window2; % Window
rectangle('Position',[W(1) W(1) W(length(W))-W(1) W(length(W))-W(1)],'EdgeColor','r');
SP_W=SP(W,W);
% padding windowed image term
pad_w=2;
SP_WP=zeros(pad_w*M); % padded windowed spectrum
SP_WP((pad_w*M/2-length(SP_W)/2+1):(pad_w*M/2+length(SP_W)/2),(pad_w*M/2-length(SP_W)/2+1):(pad_w*M/2+length(SP_W)/2))=SP_W;
figure; imshow(10*mat2gray(abs(SP_WP))); title('Filtered Image Term (1st Order)')
IH_W=fftshift(ifft2(fftshift(SP_WP)));
% IH_PHASE=IH_W((M*pad_w/2-M/2+1):(M*pad_w/2+M/2),(M*pad_w/2-M/2+1):(M*pad_w/2+M/2)); % reducing windowed hologram back to M x M
figure; imshow(mat2gray(abs(IH_W))); title('Filtered Hologram')
% figure; imshow(mat2gray(abs(IH_PHASE))); title('Filtered Hologram, M x M')

%% Reconstruction (Fresnel diffraction)
r3=1:length(IH_W);
c3=1:length(IH_W);
[C3, R3]=meshgrid(c3, r3);
EN=((C3-M*pad_w/2-1).^2+(R3-M*pad_w/2-1).^2).^0.5;
RR=EN.*dx/4;
QP=exp(1i*pi/w/z.*(RR.^2)); % Quadratic phase exponential
% QP_factor=1;
QP_factor=exp(1i*2*pi*z/w)*exp(1i*pi*((C3-M*pad_w/2-1).^2+(R3-M*pad_w/2-1).^2)/w/z)/1i/w/z;
FTS=QP_factor.*fftshift(fft2(fftshift(IH_W.*QP))); % Goodman eq 4-17
I2=FTS.*conj(FTS);
figure; imshow(5.*mat2gray(I2));
title(['Reconstructed image (Illumination ratio, Ref:F=' num2str(Ref_Ill/F_Ill) ':1)'])
%title(['Reconstructed image (Window Pixels=' num2str(WindowPixels) ')'])
%title('Reconstructed image')
axis off

%% Extracting Phase
WRAPPED_PHASE=angle(IH_W);
UNWRAPPED_PHASE=unwrap(WRAPPED_PHASE,[],1);
UNWRAP_PHASE=unwrap(UNWRAPPED_PHASE,[],2);
OP=angle(EO); UOP=unwrap(OP,[],1); Orig_PHASE=unwrap(UOP,[],2);
figure;subplot(222)
mesh(1:M,1:M,UOP);title('Original phase (radian)'); axis square
subplot(223);
x_f=1/pad_w:1/pad_w:M; y_f=x_f;
mesh(x_f,y_f,WRAPPED_PHASE);title('Recovered Phase (Wrapped)')
axis square
subplot(224); 
mesh(x_f,y_f,UNWRAPPED_PHASE);title('Recovered Phase (Unwrapped)')
axis square

%% Other Figures
figure(); subplot(221)
imshow(mat2gray(abs(I_FOV))); title('CCD FOV (Aimed at Target)')
subplot(223); imshow(mat2gray(abs(IH_W))); title('Filtered Hologram')
subplot(222); imshow(100.*mat2gray(abs(SP)));title('Hologram spectrum')
hold on; rectangle('Position',[W(1) W(1) W(length(W))-W(1) W(length(W))-W(1)],'EdgeColor','r');
hold off
subplot(224); imshow(5.*mat2gray(I2)); title('Reconstructed image')