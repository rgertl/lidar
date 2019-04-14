%% Modification Comments
% 4/14/19 - SB - I still didn't have much luck trying to apply the Z and
% theta formulas given in the textbook. Still using the values they
% provide. I've been trying to get the windowing to work (see
% oa_holog_windowing.m), but have run into an issue with that as well. I
% emailed the professor to ask about it. Haven't yet added in the
% digitization stuff we discussed last time; the ref wave amplitude is
% set on line 60. 
% 4/7/19 - SB - changed all of the indexing to dynamically update based on
% chosen values for M (# pixels), and the 3 padding factors, pad1, pad2,
% and pad3.
% Textbook uses units of cm. I'm still trying to figure out how to change
% these values and incorporate the formulas from the text book (in my CCD
% notes on the google drive). I've added a bunch of different JPEGs for us 
% to experiment with. 

format compact
close all
%% Reading in simulated target
I0=imread('airplane.jpg'); % 256x256 pixels, 8bit image
I0=double(rgb2gray(I0));

%% Parameter Setup
% length in cm
M=256; % # pixels, 1D
dx=10e-4; % pixel pitch (10 um)
w=633e-8; % wavelength (633 nm) %..actually 63.3 nm
% w=633e-9;
z=20; % propagation distance
% z=4*dx*(D+M*dx)/w % eq 4.50
angle=0.3; % reference beam angle; degree
% angle=asind(3*w/(8*dx)) % eq 4.51
res=w*z/M/dx % sampling distance (smallest resolvable element)
%% Showing CCD FOV
I_FOV=zeros(M);
I_FOV((M/2-length(I0)/2+1):(M/2+length(I0)/2),(M/2-length(I0)/2+1):(M/2+length(I0)/2))=I0;
figure(); imshow(mat2gray(I_FOV)); title('CCD FOV (Aimed at Target)')
%% Object beam propagation
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
[C2, R2]=meshgrid(c2, r2);
Ref=AV*exp(1i*2*pi*sind(theta)*dx/4.*(R2-M*pad2/2-1)/w+1i*2*pi*sind(theta)*dx/4.*(C2-M*pad2/2-1)/w);
IH=(EOf+Ref).*conj(EOf+Ref); % |F+R|^2
scale=.5; % pad3/pad2
IH=IH((M*pad2/2-M*scale*pad2/2+1):(M*pad2/2+M*scale*pad2/2),(M*pad2/2-M*scale*pad2/2+1):(M*pad2/2+M*scale*pad2/2));
figure; imshow(mat2gray(IH));title('Hologram')
SP=fftshift(fft2(fftshift(IH)));
figure; imshow(100.*mat2gray(abs(SP)));title('Hologram spectrum')
%% Windowing
% hold on
% W=81:223; % Window
% rectangle('Position',[W(1) W(1) W(length(W))-W(1) W(length(W))-W(1)],'EdgeColor','r');
% SP=SP(W,W);
% IH_W=fftshift(ifft2(fftshift(SP)));
% IH_W=IH_W.*conj(IH_W)
% figure; imshow(mat2gray(abs(IH_W)))

%% Reconstruction (Fresnel diffraction)
pad3=2;
r3=1:pad3*M;
c3=1:pad3*M;
[C3, R3]=meshgrid(c3, r3);
THOR=((R3-M*pad3/2-1).^2+(C3-M*pad3/2-1).^2).^0.5;
% r3=1:length(IH_W);
% c3=1:length(IH_W);
% [C3, R3]=meshgrid(c3, r3);
% THOR=((R3-length(IH_W)/2-1).^2+(C3-length(IH_W)/2-1).^2).^0.5;
RR=THOR.*dx/4;
QP=exp(1i*pi/w/z.*(RR.^2)); % Quadratic phase exponential
FTS=fftshift(fft2(fftshift(IH.*QP))); % Goodman eq 4-17
% FTS=fftshift(fft2(fftshift(IH_W.*QP))); % Goodman eq 4-17
I2=FTS.*conj(FTS);
figure; imshow(5.*mat2gray(I2)); title('Reconstructed image')

%% Other Figures
figure(); subplot(221)
imshow(mat2gray(I_FOV)); title('CCD FOV (Aimed at Target)')
subplot(222); imshow(mat2gray(IH));
title('Hologram')
subplot(223); imshow(100.*mat2gray(abs(SP)));title('Hologram spectrum')
subplot(224); imshow(5.*mat2gray(I2)); title('Reconstructed image')