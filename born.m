%% Born simulation
baseDir='/Users/eckert/Documents/MATLAB/CIFinalProject/';
filename='born_cellHighRes_1um.mat';

totSz=9; %um
n_r=1.33; %surrounding refractive index
matrix_res=0.1; %um
angle_illum=-0.45:0.15:0.45; %radians

lambda = 0.532/n_r; %um % wavelength of illumination, assume monochromatic
NA =0.25; % numerical aperture of the objective
mag = 12.5; %Magnification of the system
dpix_c = 6.5; %um %pixel size on the sensor plane %PCO = 6.5um
%% Load in object defined by refractive index

[object,cellsize_z,matrix_res]=cell3D(matrix_res); %Create object
%[object,cellsize_z,matrix_res]=simulObj(matrix_res); %Create object
oSz=size(object);
padAmt=round((totSz.*ones(1,2)-(oSz(1:2).*matrix_res))/(2 * matrix_res));
tic
%% Translate n(r) to scattering potential, V(r)
k0_V=2*pi/lambda;
objectV=((k0_V).^2).*( object.^2 - n_r.^2 );
%objectV=((k0_V).^2).*( object.^2 - 1);
background=0;
objectV=padarray(objectV, padAmt, background);
objSz=size(objectV);
%% Set up illumination and imaging systems
%Illumination angles
[angleU,angleV]=meshgrid(angle_illum,angle_illum); %Angles (radians)
angleU=angleU(:); %Flatten
angleV=angleV(:);
freqCent=sin([angleV angleU])./lambda; %Get k-space coordinates
freqCent(:,3)=sqrt(1-sin(angleV).^2-sin(angleU).^2)./lambda; %Get z k-space coordinate
numIllum=size(freqCent,1); %Number of illumination angles

illumNA=sqrt(sin(angleU).^2 + sin(angleV).^2);
BF=illumNA<=NA; %Get index of brightfield images

NsampR=round(totSz/(dpix_c/mag)).*ones(1,2); %Number of pixels on camera to capture the total size

FOV=objSz.*matrix_res; %Physical extent of object
du=1./FOV;%Corresponds to sampling in Fourier space

%Phsyical sampling locations of simulated object (um)
[x_obj,y_obj,z_obj]=meshgrid(matrix_res.*(-floor(objSz(2)/2):floor(objSz(2)/2)-1),matrix_res.*(-floor(objSz(1)/2):floor(objSz(1)/2)-1),matrix_res.*(-floor(objSz(3)/2):floor(objSz(3)/2)-1));
%Sampling locations in k-space of simulated object (1/um)
[u_obj,v_obj,w_obj]=meshgrid(du(2).*(-floor(objSz(2)/2):floor(objSz(2)/2)-1),du(1).*(-floor(objSz(1)/2):floor(objSz(1)/2)-1),du(3).*(-floor(objSz(3)/2):floor(objSz(3)/2)-1));
%radial locations (real space) (um)
r_obj=sqrt(x_obj.^2 + y_obj.^2 + z_obj.^2);
%coordinate of center
cent=floor(objSz/2)+1;

%Matrix for Fresnel propagator
k0=pi*lambda*(u_obj(:,:,1).^2 + v_obj(:,:,1).^2);
%Pupil in k-space 
pupil=sqrt(u_obj(:,:,1).^2 + v_obj(:,:,1).^2)<=NA/lambda;

I=zeros([objSz(1:2) numIllum]);
%I=zeros([NsampR numIllum]); %Initialize stack

%% Set up 3D Green's function for point source
G=(1/(4*pi)) .* exp(1i.*k0_V.*r_obj)./r_obj; %3D Green's function

%Deal with r=0
R=((3/(4*pi)) .*matrix_res.^3).^(1/3); % Radius of sphere around center
G0=(1/(k0_V^2)).*((1-1i*k0_V*R).*exp(1i.*k0_V.*R) - 1) *(matrix_res.^-3);
%Integral over sphere around center must == G(r=0) * dxdydz

G(cent(1),cent(2),cent(3))=G0; %Set origin equal to calculated value

%% Simulate
z0=-cellsize_z/2; %Focus halfway through volume
%***Check where this is with definition of slices!
H0 =exp(-1i*k0*z0);

sub=abs(mod(NsampR,2)-1);

downsamp= @(x) x(cent(1)-floor(NsampR(1)/2):cent(1)+floor(NsampR(1)/2)-sub(1), ...
    cent(2)-floor(NsampR(2)/2):cent(2)+floor(NsampR(2)/2)-sub(2));

F = @(x) fftshift(fft2(ifftshift(x)));
Ft = @(x) fftshift(ifft2(ifftshift(x)));
F3 = @(x) fftshift(fft(ifftshift(x),[],3));
Ft3 = @(x) fftshift(ifft(ifftshift(x),[],3));
objFilt=I;
Gkernel=F3(G).*matrix_res.^3;
for ii=1:numIllum
%     %Incident field
Ein=exp(1i.*2*pi.*(x_obj.*freqCent(ii,2) + y_obj.*freqCent(ii,1) + z_obj.*freqCent(ii,3)));
%     Ef=0;
Ei=Ein;
%     for jj=1:20  
EV=Ei.*objectV; %Multiply incident field by scattering potential
%     
    Es=Ft3(F3(EV).*Gkernel); %Convolve with psf to get scattered field
    Ef=Ein+Es; 
     Ei=Ef;
%     end
%     %Image
%     objOut=Ft(F(Ef(:,:,end)).*H0); %Backpropagate final slice output to focus plane
%     objFilt=Ft(F(objOut).*pupil);%Ft(downsamp(F(objOut).*pupil)).*NsampR(1).*NsampR(2)./(objSz(1).*objSz(2)); %Filter by the pupil
% 
%     objFilt2=Ft(F(Ef(:,:,cent(3))).*pupil);%Ft(downsamp(F(objOut).*pupil)).*NsampR(1).*NsampR(2)./(objSz(1).*objSz(2)); %Filter by the pupil
%     I(:,:,ii)=abs(objFilt).^2; %Get intensity
%     I2(:,:,ii)=abs(objFilt2).^2; %Get intensity

%     
%     
    %Image
    Obj=F(EV);
 
    wolfKernel=1i.*exp(1i.*2.*pi.*z_obj.*sqrt((lambda^-2)-abs(u_obj(:,:,1)).^2 - abs(v_obj(:,:,1)).^2))./(4.*pi.*sqrt((lambda^-2)-abs(u_obj(:,:,1)).^2 - abs(v_obj(:,:,1)).^2));
    ObjFilt3=Obj.*wolfKernel.*repmat(pupil,[1 1 objSz(3)]);
    ObjFilt=sum(ObjFilt3,3).*matrix_res;
%    objFilt(:,:,ii)=Ft(downsamp(ObjFilt)).*NsampR(1).*NsampR(2)./(objSz(1).*objSz(2));
    if BF(ii)
        %Brightfield
        objFilt=Ft(ObjFilt)+ Ein(:,:,cent(3)); 
        %objFilt2= objFilt(:,:,ii) + Ft(downsamp(F(Ein(:,:,cent(3))))).*matrix_res.^2; %Filter by the pupil
    else
        %Darkfield
        objFilt=Ft(ObjFilt); 
        %objFilt2=objFilt(:,:,ii);
    end
    I(:,:,ii)=abs(objFilt).^2; %Get intensity
    
%     
end
t_sim=toc;
%save([baseDir filename],'I','t_sim','totSz','n_r','dSlice_targ','numSlice','cellsize_z','angle_illum')
% end

