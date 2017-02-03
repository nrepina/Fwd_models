%% Set up FDTD data
baseDir='/Users/eckert/Desktop/fdtd-cell1/';
filename='fdtd_cell1_fixedAngles_HighRes_wCam.mat';

fileList=dir([baseDir '*.h5']);
numFiles=length(fileList);

totSz=100; %um
cellsize=30; %um, simulated matrix size
matrix_res=0.1; %um, simulated resolution
padAmt=round((totSz-cellsize)/(2 * matrix_res));
n_r=1.33;
matrix_res=0.1;

info5=h5info([baseDir fileList(1).name]);
objSz=circshift(info5.Datasets(1).Dataspace(1).Size,2);
objSz=objSz-20;

 objSz(1:2)=objSz(1:2)+padAmt.*2; %After padding
%% Set up illumination
angle_illum=(-0.45:0.15:0.45); %radians
lambda = 0.532/n_r; %um % wavelength of illumination, assume monochromatic
NA =0.25; % numerical aperture of the objective
mag = 12.5; %Magnification of the system
dpix_c = 6.5; %um %pixel size on the sensor plane %PCO = 6.5um

[angleU,angleV]=meshgrid(angle_illum,angle_illum);
%Switch angles because in FDTD simulation, angle-x scans faster than
%angle_y
angleU=angleU(:);
angleV=angleV(:);
freqCent=sin([angleV angleU])./lambda;
numIllum=size(freqCent,1);

NsampR=round(totSz/(dpix_c/mag)).*ones(1,2); %Number of pixels on camera to capture the total size

FOV=objSz.*matrix_res; %Physical extent of object
du=1./FOV;%Corresponds to sampling in Fourier space

k0_V=2*pi/lambda;

%Phsyical sampling locations of simulated object
[x_obj,y_obj]=meshgrid(matrix_res.*(-floor(objSz(2)/2):floor(objSz(2)/2)-1),matrix_res.*(-floor(objSz(1)/2):floor(objSz(1)/2)-1));
%Sampling locations in k-space of simulated object
[u_obj,v_obj]=meshgrid(du(2).*(-floor(objSz(2)/2):floor(objSz(2)/2)-1),du(1).*(-floor(objSz(1)/2):floor(objSz(1)/2)-1));
%Matrix for Fresnel propagator
k0=pi*lambda*(u_obj.^2 + v_obj.^2);
%Pupil in k-space 
pupil=sqrt(u_obj.^2 + v_obj.^2)<=NA/lambda;

I=zeros([NsampR numIllum]); %Initialize stack

tic
%% Load and process
cent=floor(objSz/2)+1;
downsamp= @(x) x(cent(1)-floor(NsampR(1)/2):cent(1)+floor(NsampR(1)/2)-1, ...
    cent(2)-floor(NsampR(2)/2):cent(2)+floor(NsampR(2)/2)-1);

z_sel=10.5; %10.5 um from beginning (1.5 from end)
z0=-4.5; %um, refocus to middle of cell
H0 =exp(1i.*k0_V.*z0).*exp(-1i*k0*z0);

%I=zeros([objSz(1:2) numIllum]);
for ii = 33;%1:numIllum
    %Load
    Er=h5read([baseDir fileList(ii).name],'/ex.r');
    Ei=h5read([baseDir fileList(ii).name],'/ex.i');
    
    Ex=Er + 1i.*Ei;
    Ex=permute(Ex,[2 3 1]); %Put z in third dim
    Ex=Ex(10:309,10:309,10:129); %Get rid of PML
    %sliderDisplayIm(abs(Ex));
    temp(:,:,ii)=abs(permute(Ex(:,150,:),[3 1 2]));
    Ex_final=Ex(:,:,round(z_sel/matrix_res));
    
    Ex_large=exp(1i.*2*pi.*(x_obj.*freqCent(ii,2) + y_obj.*freqCent(ii,1) ));
    % Create image for each angle
    background=abs(Ex_final(80,80));
    %Ex_large(padAmt:end-padAmt-1,padAmt:end-padAmt-1)=Ex_final;
    Ex_large=padarray(Ex_final, [padAmt padAmt], background);
    %temp(:,:,ii)=F(Ex_large);
    obj=Ft(downsamp(F(Ex_large).*H0.*pupil)); %propagate back, filter by pupil, downsample
    %obj=Ft(F(Ex_final).*H0.*pupil); %propagate back, filter by pupil
    I(:,:,ii)=abs(obj).^2;
    
end
t_sim=toc;

save([baseDir filename],'I','t_sim','totSz','n_r','dSlice_targ','numSlice','cellsize_z','angle_illum')

