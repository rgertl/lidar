close all;
clear all
%% Reading input bitmap file
I0=imread('airplane.jpg'); % 256x256 pixels, 8bit image
I0=double(rgb2gray(I0));
figure()
imshow(mat2gray(I0));axis off
% parameter setup
M=256;
deltax=0.001; % pixel pitch 0.001 cm (10 um)
w=633*10^-8; % wavelength 633 nm
z=20; % z=M*deltax^2/w; % propagation distance
% Step 1: simulation of propagation
r=1:5*M;
c=1:5*M;
[C, R]=meshgrid(c, r);
I=zeros(5*M);
I(513:768,513:768)=I0;
A0=fftshift(ifft2(fftshift(I)));
deltaf=1/5/M/deltax;
p=exp(-2i*pi*z.*((1/w)^2-((R-641).*deltaf).^2-((C-641).*deltaf).^2).^0.5);
Az=A0.*p;
EO=fftshift(fft2(fftshift(Az)));
EO=EO(513:768,513:768); % reduce diffraction-plane size
% Step 2: interference at the hologram plane
% zero-padding in the spectrum domain
Az=fftshift(ifft2(fftshift(EO)));
Az2=zeros(4*M);
Az2(385:640,385:640)=Az;
EOf=fftshift(fft2(fftshift(Az2)));
AV=(min(min(abs(EOf)))+max(max(abs(EOf))))/2;
angle=0.3; % reference beam angle; degree
r2=1:4*M;
c2=1:4*M;
[C2, R2]=meshgrid(c2, r2);
Ref=AV*exp(1i*2*pi*sind(angle)*deltax/4.*(R2-2*M-1)/w+1i*2*pi*sind(angle)*deltax/4.*(C2-2*M-1)/w);
IH=(EOf+Ref).*conj(EOf+Ref);
IH=IH(257:768,257:768); % reduce the hologram size
figure; imshow(mat2gray(IH));
title('Hologram')
axis off
SP=fftshift(ifft2(fftshift(IH)));
figure; imshow(50.*mat2gray(abs(SP)));
title('Hologram spectrum')
axis off
% Step 3: reconstruction (Fresnel diffraction)
r3=1:2*M;
c3=1:2*M;
[C3, R3]=meshgrid(c3, r3);
THOR=((R3-M-1).^2+(C3-M-1).^2).^0.5;
A=THOR.*deltax/4;
QP=exp(1i*pi/w/z.*(A.^2));
FTS=fftshift(fft2(fftshift(IH.*QP)));
I2=FTS.*conj(FTS);
figure; imshow(5.*mat2gray(I2));
title('Reconstructed image')
axis off