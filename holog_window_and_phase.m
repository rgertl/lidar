format compact
close all
clear all
%% Reading in simulated object
I0=imread('airplane.jpg'); % 256x256 pixels, 8bit image
I0=double(rgb2gray(I0));
%% Simulating Target Phase
[x,y]=meshgrid(1:256,1:256);
PHASE=15*sin(((x-200).^2+(y-225).^2)/10000)+0.002*((x-37).^2+(y-100).^2);
I0=I0.*exp(1i*PHASE);

%% Parameter Setup
% length in cm
M=256; % # pixels, 1D
dx=10e-4; % pixel pitch (10 um)
w=633e-8; % wavelength (633 nm) %..actually 63.3 nm
% w=633e-9;
z=20; % propagation distance
% z=4*dx*(D+M*dx)/w % eq 4.50
theta=0.4; % reference beam angle; degree
% angle=asind(3*w/(8*dx)) % eq 4.51
res=w*z/M/dx % sampling distance (smallest resolvable element)
D=72*dx;
snr=.1; % White gaussian noise
G_AD=9.76; % e-/px (Gain)
F_Ill=1; % e-/px
Ref_Ill=1; % e-/px
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
I((M*pad1/2-length(I0)/2+1):(M*pad1/2+length(I0)/2),(M*pad1/2-length(I0)/2+1):(M*pad1/2+length(I0)/2))=I0;
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
AV=(min(min(abs(EOf)))+max(max(abs(EOf))))/2; % ref wave amplitude
[C2, R2]=meshgrid(c2, r2);
Ref=Ref_Ill*exp(1i*2*pi*sind(theta)*dx/4.*(R2-M*pad2/2-1)/w+1i*2*pi*sind(theta)*dx/4.*(C2-M*pad2/2-1)/w); % eq 3.5a - DH textbook, with scaled amplitude set to 32000 e-/px
IH=(EOf+Ref).*conj(EOf+Ref)/G_AD+awgn(real(Ref),snr,'measured')/G_AD; % |F+R|^2 + noise
scale=.5; % pad3/pad2
IH=IH((M*pad2/2-M*scale*pad2/2+1):(M*pad2/2+M*scale*pad2/2),(M*pad2/2-M*scale*pad2/2+1):(M*pad2/2+M*scale*pad2/2));
figure; imshow(mat2gray(IH)); title('Hologram')
SP=fftshift(fft2(fftshift(IH)));
figure; imshow(100.*mat2gray(abs(SP)));title('Hologram spectrum')
%% Windowing
hold on
W=51:180; % Window
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
% pad3=2;
% r3=1:pad3*M;
% c3=1:pad3*M;
% [C3, R3]=meshgrid(c3, r3);
% THOR=((R3-M*pad3/2-1).^2+(C3-M*pad3/2-1).^2).^0.5;
r3=1:length(IH_W);
c3=1:length(IH_W);
[C3, R3]=meshgrid(c3, r3);
THOR=(R3.^2+C3.^2).^0.5;
RR=THOR.*dx/4;
QP=exp(1i*pi/w/z.*(RR.^2)); % Quadratic phase exponential
QP_factor=1;
QP_factor=exp(1i*2*pi*z/w)*exp(1i*pi*(R3.^2+C3.^2)/w/z)/1i/w/z;
FTS=QP_factor.*fftshift(fft2(fftshift(IH_W.*QP))); % Goodman eq 4-17
I2=FTS.*conj(FTS);
figure; imshow(5.*mat2gray(I2));
title('Reconstructed image')
axis off

%% Extracting Phase
% Sfreq = (-1/2:1/(pad_w*M):1/2-1/(pad_w*M));
% [Sx,Sy] = meshgrid(Sfreq,Sfreq); 
WRAPPED_PHASE=angle(IH_W);
UNWRAPPED_PHASE=unwrap(WRAPPED_PHASE,[],1);
UNWRAP_PHASE=unwrap(UNWRAPPED_PHASE,[],2);
figure;subplot(221)
mesh(x,y,PHASE);title('Original phase (radian)')
subplot(223); 
mesh(C3,R3,WRAPPED_PHASE);title('Recovered Phase (Wrapped)')
subplot(224); 
mesh(C3,R3,UNWRAPPED_PHASE);title('Recovered Phase (Unwrapped)')

%% Extracting Phase
% % Sfreq = (-1/2:1/(pad_w*M):1/2-1/(pad_w*M));
% % [Sx,Sy] = meshgrid(Sfreq,Sfreq); 
% WRAPPED_PHASE=angle(IH_PHASE);
% UNWRAPPED_PHASE=unwrap(WRAPPED_PHASE,[],1);
% UNWRAP_PHASE=unwrap(UNWRAPPED_PHASE,[],2);
% figure;subplot(221)
% mesh(x,y,PHASE);title('Original phase (radian)')
% subplot(223); 
% mesh(x,y,WRAPPED_PHASE);title('Recovered Phase (Wrapped)')
% subplot(224); 
% mesh(x,y,UNWRAPPED_PHASE);title('Recovered Phase (Unwrapped)')
%% Other Figures
figure(); subplot(221)
imshow(mat2gray(abs(I_FOV))); title('CCD FOV (Aimed at Target)')
subplot(223); imshow(mat2gray(abs(IH_W))); title('Filtered Hologram')
subplot(222); imshow(500.*mat2gray(abs(SP)));title('Hologram spectrum')
hold on; rectangle('Position',[W(1) W(1) W(length(W))-W(1) W(length(W))-W(1)],'EdgeColor','r');
hold off
subplot(224); imshow(5.*mat2gray(I2)); title('Reconstructed image')