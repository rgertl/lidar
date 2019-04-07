%% Modification Comments
% 4/7/19 - SB - changed all of the indexing to dynamically update based on
% chosen values for M (# pixels), and the 3 padding factors, pad1, pad2,
% and pad3.
% Textbook uses units of cm. I'm still trying to figure out how to change
% these values and incorporate the formulas from the text book (in my CCD
% notes on the google drive). I think it will just be a function of playing
% with the padding (lines 31,46,68)....or maybe we just need the image to 
% be bigger. When I try to scale up the number of pixels, the results
% become extremely small (only a few pixels). 

%% Reading in simulated object
I0=imread('airplane.jpg'); % 256x256 pixels, 8bit image
I0=double(rgb2gray(I0));
% figure()
% imshow(mat2gray(I0));
axis off
%% parameter setup 
% D=.01; % [meters]
% M=512; % # pixels, 1D
M=256; % p. 113 - off-axis CCD must be 4x bandwidth
dx=10e-4; % pixel pitch (10 um)
w=633e-8; % wavelength (633 nm) %..actually 63.3 nm
% w=633e-9;
z=20; % z=M*deltax^2/w; % propagation distance
% z=4*dx*(D+M*dx)/w % eq 4.50
angle=0.3; % reference beam angle; degree
% angle=asind(3*w/(8*dx)) % eq 4.51
res=w*z/M/dx % sampling distance (smallest resolvable element)
%% Step 1: simulation of propagation
pad1=5; % Padding
r=1:pad1*M;
c=1:pad1*M;
[C, R]=meshgrid(c, r);
I=zeros(5*M);
% I(513:768,513:768)=I0;
I((M*pad1/2-length(I0)/2+1):(M*pad1/2+length(I0)/2),(M*pad1/2-length(I0)/2+1):(M*pad1/2+length(I0)/2))=I0;
A0=fftshift(ifft2(fftshift(I)));
deltaf=1/pad1/M/dx;
p=exp(-2i*pi*z.*((1/w)^2-((R-M*pad1/2-1).*deltaf).^2-((C-M*pad1/2-1).*deltaf).^2).^0.5);
Az=A0.*p;
EO=fftshift(fft2(fftshift(Az)));
EO=EO((M*pad1/2-M/2+1):(M*pad1/2+M/2),(M*pad1/2-M/2+1):(M*pad1/2+M/2)); % reduce diffraction-plane size
%% Step 2: interference at the hologram plane
% zero-padding in the spectrum domain
pad2=4;
r2=1:pad2*M;
c2=1:pad2*M;
Az=fftshift(ifft2(fftshift(EO)));
Az2=zeros(4*M);
Az2((M*pad2/2-M/2+1):(M*pad2/2+M/2),(M*pad2/2-M/2+1):(M*pad2/2+M/2))=Az;
EOf=fftshift(fft2(fftshift(Az2)));
AV=(min(min(abs(EOf)))+max(max(abs(EOf))))/2;
[C2, R2]=meshgrid(c2, r2);
Ref=AV*exp(1i*2*pi*sind(angle)*dx/4.*(R2-M*pad2/2-1)/w+1i*2*pi*sind(angle)*dx/4.*(C2-M*pad2/2-1)/w);
IH=(EOf+Ref).*conj(EOf+Ref);
scale=.5; % crop hologram, 0 = no hologram, 1 = full size (condition: scale*pad2 = integer)
% IH=IH(257:768,257:768); % reduce the hologram size
IH=IH((M*pad2/2-M*scale*pad2/2+1):(M*pad2/2+M*scale*pad2/2),(M*pad2/2-M*scale*pad2/2+1):(M*pad2/2+M*scale*pad2/2));
figure; imshow(mat2gray(IH));
title('Hologram')
axis off
SP=fftshift(ifft2(fftshift(IH)));
figure; imshow(50.*mat2gray(abs(SP)));
title('Hologram spectrum')
axis off
%% Step 3: reconstruction (Fresnel diffraction)
pad3=2;
r3=1:pad3*M;
c3=1:pad3*M;
[C3, R3]=meshgrid(c3, r3);
THOR=((R3-M*pad3/2-1).^2+(C3-M*pad3/2-1).^2).^0.5;
A=THOR.*dx/4;
QP=exp(1i*pi/w/z.*(A.^2));
FTS=fftshift(fft2(fftshift(IH.*QP)));
I2=FTS.*conj(FTS);
figure; imshow(5.*mat2gray(I2));
title('Reconstructed image')
axis off