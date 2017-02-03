% clear all;
% set(0,'DefaultFigureWindowStyle','Docked');
% warning off;
global M N O Fx Fy fGreen E U_m I_m denominator prop_wavenumber

baseDir='/Users/eckert/Documents/MATLAB/CIFinalProject/';
filename='born_cellHighRes_1um_new.mat';

NA = 0.25;
mag = 12.5;
lambda = 0.532;
%ps=1.25/mag;
ps = 6.5/mag;

psz = ps;
wavenumber = 2*pi/lambda;

F2 = @(x) fftshift(fft2(ifftshift(x)));
Ft2 = @(x) fftshift(ifft2(ifftshift(x)));


%% Object Creation

F = @(x) fftn(x);
IF = @(x) ifftn(x);

matrix_res=0.1;
[Bead1,cellsize_z,matrix_res]=cell3D(matrix_res); %Create object

Bead=fftshift(Bead1);

RI=1.33;
% RI = 1.59;
% del_n = 0.01;
% radius = 2;
% Bead = zeros(size(X));
% Bead((X/radius/0.8).^2+(Y/radius/1.2).^2 + (Z/radius/0.6).^2 <= 1) = RI+del_n;
% Bead(X.^2+(Y-0.5).^2 + (Z-0.5).^2 <= (0.1*radius)^2) = RI+2*del_n;
% Bead(X.^2+(Y+0.5).^2 + (Z+0.5).^2 <= (0.2*radius)^2) = RI+0.5*del_n;
% Bead(Bead==0) = RI;

M = size(Bead,2);
N = size(Bead,1);
O = size(Bead,3);



% x = -(N-mod(N,2))/2:1:(N-mod(N,2))/2-(mod(N,2)==0);
% x = ps*x;
% y = -(M-mod(M,2))/2:1:(M-mod(M,2))/2-(mod(M,2)==0);
% y = ps*y;
% z = -(O-mod(O,2))/2:1:(O-mod(O,2))/2-(mod(O,2)==0);
% z = psz*z;


x = -(N-mod(N,2))/2:1:(N-mod(N,2))/2-(mod(N,2)==0);
x = matrix_res*x;
y = -(M-mod(M,2))/2:1:(M-mod(M,2))/2-(mod(M,2)==0);
y = matrix_res*y;
z = -(O-mod(O,2))/2:1:(O-mod(O,2))/2-(mod(O,2)==0);
z = matrix_res*z;

x = ifftshift(x);
y = ifftshift(y);
z = ifftshift(z);


[X,Y,Z] = meshgrid(x,y,z);

V = wavenumber^2*(RI^2-Bead.^2);
Vs = fftshift(V);

figure;imagesc(fftshift(z),fftshift(y),squeeze(fftshift(Bead(:,1,:))));axis image;colormap jet;
tic
%% Green's Function

M2 = 2*M;N2 = 2*N;O2 = 2*O;
dx = x(2)-x(1);dy = y(2)-y(1);dz = z(2)-z(1);
% x2 = -(N2-mod(N2,2))/2:1:(N2-mod(N2,2))/2-(mod(N2,2)==0);
% x2 = ps*x2;
% y2 = -(M2-mod(M2,2))/2:1:(M2-mod(M2,2))/2-(mod(M2,2)==0);
% y2 = ps*y2;
% z2 = -(O2-mod(O2,2))/2:1:(O2-mod(O2,2))/2-(mod(O2,2)==0);
% z2 = psz*z2;
x2 = -(N2-mod(N2,2))/2:1:(N2-mod(N2,2))/2-(mod(N2,2)==0);
x2 = matrix_res*x2;
y2 = -(M2-mod(M2,2))/2:1:(M2-mod(M2,2))/2-(mod(M2,2)==0);
y2 = matrix_res*y2;
z2 = -(O2-mod(O2,2))/2:1:(O2-mod(O2,2))/2-(mod(O2,2)==0);
z2 = matrix_res*z2;

x2 = ifftshift(x2); y2 = ifftshift(y2); z2 = ifftshift(z2);

[X2,Y2,Z2] = meshgrid(x2,y2,z2);
[X3,Y3,Z3] = meshgrid(x2,y2,fftshift(z));
clear x2 y2

Green = exp(1i*2*pi*RI/lambda*sqrt(X2.^2+Y2.^2+Z2.^2))./(4*pi*sqrt(X2.^2+Y2.^2+Z2.^2+eps));
R = (dx*dy*dz*3/4/pi)^(1/3);
Green(1,1,1) = (exp(1i*2*pi*RI/lambda*R)*(1-1i*2*pi*RI/lambda*R)-1)/(2*pi*RI/lambda)^2/(dx*dy*dz);
Green = -1*Green;
fGreen = F(Green)*dx*dy*dz;

%% Ewald Sphere Coordinates

%dfx = 1/N2/ps; dfy = 1/M2/ps;
dfx = 1/(N2.*matrix_res); dfy = 1/(M2.*matrix_res);
fx = dfx*(-(N2-mod(N2,2))/2:1:(N2-mod(N2,2))/2-(mod(N2,2)==0));
fy = dfy*(-(M2-mod(M2,2))/2:1:(M2-mod(M2,2))/2-(mod(M2,2)==0));
fx = ifftshift(fx);
fy = ifftshift(fy);
[Fx,Fy] = meshgrid(fx,fy);
Fx = repmat(Fx,[1,1,O]);  Fy = repmat(Fy,[1,1,O]);


P = (Fx(:,:,1).^2 + Fy(:,:,1).^2 <= (NA/lambda)^2); %define pupil
G = -1i.*P./sqrt((RI/lambda)^2-(Fx(:,:,1).^2+Fy(:,:,1).^2))/4/pi; %2D green's*P
G(isnan(G))=0;
prop_wavenumber = 1i*2*pi*(P.*sqrt((1/lambda)^2-(Fx(:,:,1).^2+Fy(:,:,1).^2))); %z wavevector
P = repmat(P,[1,1,O]);
G = repmat(G,[1,1,O]);
%E = psz*G.*exp(-1i*2*pi*Z3.*P.*sqrt((RI/lambda)^2-(Fx.^2+Fy.^2)));
E = matrix_res*G.*exp(-1i*2*pi*Z3.*P.*sqrt((RI/lambda)^2-(Fx.^2+Fy.^2)));

%clear P G

fx2=(1/(N*matrix_res))*(-(N-mod(N,2))/2:1:(N-mod(N,2))/2-(mod(N,2)==0));
[Fx2,Fy2]=meshgrid(fx2);
P2 = (Fx2(:,:,1).^2 + Fy2(:,:,1).^2 <= (NA/lambda)^2); %define pupil


%% different forward models
angle_illum=-0.45:0.15:0.45;
%Illumination angles
[angleU,angleV]=meshgrid(angle_illum,angle_illum); %Angles (radians)
angleU=angleU(:); %Flatten
angleV=angleV(:);
freqCent=sin([angleV angleU])./lambda; %Get k-space coordinates
freqCent(:,3)=sqrt(1-sin(angleV).^2-sin(angleU).^2)./lambda; %Get z k-space coordinate
numIllum=size(freqCent,1); %Number of illumination angles
I=zeros(N,M,numIllum);
UDS=I;
I2=I;
z_plane=0;
for ii=1:numIllum
fx_in = sin(angleU(ii))/lambda;
fy_in = sin(angleV(ii))/lambda;

U_in = exp(1i*2*pi*(fx_in*X+fy_in*Y+sqrt((RI/lambda)^2-(fx_in^2+fy_in^2))*Z));
U_in = fftshift(U_in);
U_in_air = zeros(M,N,length(z_plane));
if (fx_in^2 + fy_in^2)*lambda^2 <= NA^2
    for zIdx = 1:length(z_plane)
        U_in_air(:,:,zIdx) = exp(1i*2*pi*(fx_in*X(:,:,1)+fy_in*Y(:,:,1)+sqrt((1/lambda)^2-(fx_in^2+fy_in^2))*z_plane(zIdx)));
        U_in_air(:,:,zIdx) = fftshift(U_in_air(:,:,zIdx));
    end
end

mode = 'BA';
U_temp = zeros(size(Fx));
U_farfield = zeros(size(Vs));
%U_k = U_k(:,:,:,end);

switch mode
    case 'BA'
        U_temp(1:M,1:N,:) = U_in.*Vs;
    case 'RA'
        U_temp(1:M,1:N,:) = U_in.*Vs;      
    case 'MS'
        U_temp(1:M,1:N,:) = U_k.*Vs;
end

U_temp = sum(fft2(U_temp).*E,3);

for zIdx = 1:O
    temp = ifft2(U_temp.*exp(z(zIdx)*prop_wavenumber));
    U_farfield(:,:,zIdx) = temp(1:M,1:N);
end

U_farfield = fftshift(U_farfield,3);
U_0_air = exp(1i*2*pi*(fx_in*X+fy_in*Y+sqrt((1/lambda)^2-(fx_in^2+fy_in^2))*Z));
U_0_air = fftshift(U_0_air);

if strcmp(mode,'RA')
    U_farfield = U_0_air.*exp(U_farfield./U_0_air);
else
    U_farfield = U_0_air + U_farfield;
end

%U_ff_temp=Ft(F(U_farfield(:,:,round(O/2+1))).*P2);
UDS(:,:,ii)=U_farfield(:,:,round(O/2+1));
I(:,:,ii)=abs(U_farfield(:,:,round(O/2+1))).^2;
%I(:,:,ii)=abs(U_ff_temp).^2;
end
t_sim=toc;
n_r=RI;
save([baseDir filename],'I','t_sim','n_r','angle_illum')

