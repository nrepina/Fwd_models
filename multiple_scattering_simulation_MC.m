clear all;
%set(0,'DefaultFigureWindowStyle','Docked');
warning off;
global M N O Fx Fy fGreen E U_m I_m denominator prop_wavenumber

baseDir='/Users/eckert/Documents/MATLAB/CIFinalProject/';
filename='recborn_cellHighRes_1um_new.mat';

NA = 0.25;
mag = 12.5;
lambda = 0.532;
ps = 6.5/mag;    psz = ps;
wavenumber = 2*pi/lambda;

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


P = (Fx(:,:,1).^2 + Fy(:,:,1).^2 <= (NA/lambda)^2);
G = -1i*P./sqrt((RI/lambda)^2-(Fx(:,:,1).^2+Fy(:,:,1).^2))/4/pi;
G(isnan(G))=0;
prop_wavenumber = 1i*2*pi*(P.*sqrt((1/lambda)^2-(Fx(:,:,1).^2+Fy(:,:,1).^2)));
P = repmat(P,[1,1,O]);
G = repmat(G,[1,1,O]);
%E = psz*G.*exp(-1i*2*pi*Z3.*P.*sqrt((RI/lambda)^2-(Fx.^2+Fy.^2)));
E = matrix_res*G.*exp(-1i*2*pi*Z3.*P.*sqrt((RI/lambda)^2-(Fx.^2+Fy.^2)));
clear P G

fx2=(1/(N*matrix_res))*(-(N-mod(N,2))/2:1:(N-mod(N,2))/2-(mod(N,2)==0));
[Fx2,Fy2]=meshgrid(fx2);
P2 = (Fx2(:,:,1).^2 + Fy2(:,:,1).^2 <= (NA/lambda)^2); %define pupil


%% Nestrov Iterative Forward model (farfield)


angle_illum=-0.45:0.15:0.45;
%Illumination angles
[angleU,angleV]=meshgrid(angle_illum,angle_illum); %Angles (radians)
angleU=angleU(:); %Flatten
angleV=angleV(:);
freqCent=sin([angleV angleU])./lambda; %Get k-space coordinates
freqCent(:,3)=sqrt(1-sin(angleV).^2-sin(angleU).^2)./lambda; %Get z k-space coordinate
numIllum=size(freqCent,1); %Number of illumination angles
I=zeros(N,M,numIllum);


angle_step = 1;
NA_illu = 0;
z_plane = [0];%fftshift(z);
max_iteration = 100;
tol = 5e-7;

figure;
Usave=I;
for mIdx = 1:numIllum
    
%     if angle_step > 1
%         x_step = mod(mIdx-1,angle_step);
%         NA_x = -NA_illu + x_step*2*NA_illu/(angle_step-1);
%         y_step = floor((mIdx-1)/angle_step);
%         NA_y = -NA_illu + y_step*2*NA_illu/(angle_step-1);
%     else
%         NA_x = NA_illu;
%         NA_y = NA_illu;
%     end

    theta_x = angleU(mIdx); theta_y = angleV(mIdx);
    fx_in = sin(theta_x)/lambda;
    fy_in = sin(theta_y)/lambda;

    U_in = exp(1i*2*pi*(fx_in*X+fy_in*Y+sqrt((RI/lambda)^2-(fx_in^2+fy_in^2))*Z));
    U_in = fftshift(U_in);
    U_in_air = zeros(M,N,length(z_plane));
    if (fx_in^2 + fy_in^2)*lambda^2 <= NA^2
        for zIdx = 1:length(z_plane)
            U_in_air(:,:,zIdx) = exp(1i*2*pi*(fx_in*X(:,:,1)+fy_in*Y(:,:,1)+sqrt((1/lambda)^2-(fx_in^2+fy_in^2))*z_plane(zIdx)));
            U_in_air(:,:,zIdx) = fftshift(U_in_air(:,:,zIdx));
        end
    end
    U_ini = U_in;

    [U_m, U_k, alpha, beta] = Nestrov_forward_farfield(Vs,U_ini,U_in,U_in_air,z_plane,max_iteration,tol);
    I(:,:,mIdx) = abs(U_m).^2;
    Usave(:,:,mIdx)=U_m;
    fprintf('illumination %d/%d done!\n',mIdx,angle_step^2);
    
%     subplot 131
%     imagesc(fftshift(x),fftshift(y),abs(U_m(:,:,1)).^2);
%     axis image;colormap gray;caxis([0.5,1.5]);
%     title('Intensity','fontsize',24);
%     subplot 132
%     imagesc(fftshift(x),fftshift(y),squeeze(angle(U_m(:,:,1)./U_in_air(:,:,1))));
%     title('Phase','fontsize',24);
%     axis image;colormap gray;caxis([0,1]);
%     subplot 133
%     plot(fftshift(y),squeeze(angle(U_m(:,N/2+1,1)./U_in_air(:,N/2+1,1))),['b','o'],fftshift(y),sum(2*pi/lambda*(fftshift(squeeze(Bead(:,1,:)))-RI)*psz,2),['r','--']);
%     axis square;
%     title('1D Phase','fontsize',24);
%     drawnow;
end
t_sim=toc;
n_r=RI;
save([baseDir filename],'I','n_r','angle_illum','Usave')

totSz=100; %um
oSz=size(Usave);
padAmt=round((totSz-oSz(1:2).*matrix_res)/(2 * matrix_res));

%Uair(:,:,zIdx) = exp(1i*2*pi*(fx_in*X(:,:,1)+fy_in*Y(:,:,1)+sqrt((1/lambda)^2-(fx_in^2+fy_in^2))*z_plane(zIdx)));
obj=padarray(Usave,padAmt,0);
objSz=size(obj);
cent=floor(objSz(1:2)/2)+1;
NsampR=round(totSz/(dpix_c/mag)).*ones(1,2); %Number of pixels on camera to capture the total size
objF=F(obj);
objF2=objF(cent(1)-floor(NsampR(1)/2):cent(1)+floor(NsampR(1)/2)-1, ...
    cent(2)-floor(NsampR(2)/2):cent(2)+floor(NsampR(2)/2)-1,:);
obj2=Ft(objF2);

I2=abs(obj2).^2;
