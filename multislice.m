%% Multislice simulation
baseDir='/Users/eckert/Documents/MATLAB/CIFinalProject/';
filename='multislice_cellHighRes_1um_wCam';

totSz=100; %um
dSlice_targ=1; %um; target thickness of slices
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

tic
%% Slice it up for multislice
%del(phi) = (2*pi/lambda)*OPL
%OPL=integral(n*ds)
oSz=size(object);

padAmt=round((totSz-oSz(1:2).*matrix_res)/(2 * matrix_res));
numSlice=round(cellsize_z/dSlice_targ);

%Indices to slice at
sliceE1=round(linspace(1,oSz(3)+1,numSlice+1));
sliceE2=sliceE1(2:end)-1;
sliceE1(end)=[];
nds=object.*matrix_res;

dslice=(sliceE2-sliceE1+1).*matrix_res; %distance between each slice;

OPL=zeros([oSz(1:2) numSlice]);
for jj=1:numSlice
    OPL(:,:,jj)=sum(nds(:,:,sliceE1(jj):sliceE2(jj)),3); %integrate over n*ds for each slice to get OPL
end

ph = (2*pi/lambda).*OPL; %phase of each slice
ampl=ones(size(ph)); %Completely transmissive

obj=ampl.*exp(1i.*ph); %create amplitude, phase slices

background=1.*exp(1i.*(2*pi/lambda).*n_r.*dslice(1)); %Calculate value to pad with

%Create larger space
obj=padarray(obj,padAmt,background);

objSz=size(obj);
objSz=objSz(1:2);
%% Set up illumination and imaging systems

%Illumination angles
[angleU,angleV]=meshgrid(angle_illum,angle_illum); %Angles (radians)
angleU=angleU(:); %Flatten
angleV=angleV(:);
freqCent=sin([angleV angleU])./lambda; %Get k-space coordinates
numIllum=size(freqCent,1); %Number of illumination angles

NsampR=round(totSz/(dpix_c/mag)).*ones(1,2); %Number of pixels on camera to capture the total size

FOV=objSz.*matrix_res; %Physical extent of object
du=1./FOV;%Corresponds to sampling in Fourier space

%Phsyical sampling locations of simulated object
[x_obj,y_obj]=meshgrid(matrix_res.*(-floor(objSz(2)/2):floor(objSz(2)/2)-1),matrix_res.*(-floor(objSz(1)/2):floor(objSz(1)/2)-1));
%Sampling locations in k-space of simulated object
[u_obj,v_obj]=meshgrid(du(2).*(-floor(objSz(2)/2):floor(objSz(2)/2)-1),du(1).*(-floor(objSz(1)/2):floor(objSz(1)/2)-1));
%Matrix for Fresnel propagator
k0=pi*lambda*(u_obj.^2 + v_obj.^2);
%Pupil in k-space 
pupil=sqrt(u_obj.^2 + v_obj.^2)<=NA/lambda;

I=zeros([NsampR numIllum]); %Initialize stack
%I=zeros([objSz(1:2) numIllum]);
%% Simulate
z0=-cellsize_z/2; %Focus halfway through volume
%***Check where this is with definition of slices!
H0 =exp(-1i*k0*z0);

cent=floor(objSz/2)+1;
downsamp= @(x) x(cent(1)-floor(NsampR(1)/2):cent(1)+floor(NsampR(1)/2)-1, ...
    cent(2)-floor(NsampR(2)/2):cent(2)+floor(NsampR(2)/2)-1);

F = @(x) fftshift(fft2(ifftshift(x)));
Ft = @(x) fftshift(ifft2(ifftshift(x)));
    
for ii=1:numIllum
    i0 = exp(1i.*2.*pi.*(x_obj.*freqCent(ii,2)+y_obj.*freqCent(ii,1)));

    phi=zeros([objSz numSlice]); %incident
    psi=phi; %outgoing

    phi(:,:,1)=i0;
    psi(:,:,1)=phi(:,:,1).*obj(:,:,1);

    for jj=2:numSlice
        H =exp(-1i*k0*dslice(jj)); 
        phi(:,:,jj)=Ft(F(psi(:,:,jj-1)).*H); %propagate
        psi(:,:,jj)=phi(:,:,jj).*obj(:,:,jj); %multiply by new slice
    end

    %Image
    objOut=Ft(F(psi(:,:,numSlice)).*H0.*pupil); %Backpropagate final slice output to focus plane
    %objFilt=objOut;
    objFilt=Ft(downsamp(F(objOut).*pupil)).*NsampR(1).*NsampR(2)./(objSz(1).*objSz(2)); %Filter by the pupil

    I(:,:,ii)=abs(objFilt).^2; %Get intensity
end
t_sim=toc;

% NsampR=[200 200];
% downsamp= @(x) x(cent(1)-floor(NsampR(1)/2):cent(1)+floor(NsampR(1)/2)-1, ...
% cent(2)-floor(NsampR(2)/2):cent(2)+floor(NsampR(2)/2)-1,:);
% I_crop=downsamp(I);

save([baseDir filename '.mat'],'I','I_crop','t_sim','totSz','n_r','dSlice_targ','numSlice','cellsize_z','angle_illum')

if 0
%% Reconstruct
illumNA=n_r.*sqrt(sin(angleU).^2 + sin(angleV).^2);
[illumNA_sorted,idx]=sort(illumNA);

I_r=I(:,:,idx);
freqCent_r=freqCent(idx,:);

[nObj, w_NA, x,y, x0,y0, u0, v0,k0,u,v,k,xObj,dx_obj, xImg, dx_img]=defineSystem(lambda, NA, mag, dpix_c, NsampR, illumNA_sorted);

con=NsampR(1).*dpix_c./mag; %Conversion factor (pixels/(1/um))

%Define processing options - included in paramFile
%Further parameters that depend on image parameters are also defined below
%Iteration parameters
opts.tol = 1; %tol: maximum change of error allowed in two consecutive iterations
opts.maxIter = 35; %maxIter: maximum iterations 
opts.minIter = 2; %minIter: minimum iterations
opts.monotone = 1; %monotone (1, default): if monotone, error has to monotonically dropping when iters>minIter
opts.iters = 1; %Initialize iteration counter
%Display parameters
opts.display = 1; %display: display results (0: no (default) 1: yes)
opts.mode = 'real'; %mode: display in 'real' space or 'fourier' space
%Propagation type
opts.Prop_mode=0; %Propagation mode; if 0, Fresnel propagation; otherwise, angular spectrum
%Regularization & other processing controls
opts.OP_alpha = 100; %OP_alpha: regularization parameter for Object
opts.OP_beta = 1; %OP_beta: regularization parameter for Pupil
opts.poscalibrate =0;%'sa';%0; %poscalibrate: flag for LED position correction using
    % '0': no correction
    % 'sa': simulated annealing method      
    % 'ga': genetic algorithm
    % caution: takes considerably longer time to compute a single iteration with correction      
opts.calbratetol = 1e-1; %calbratetol: parameter in controlling error tolerance in sa (simulated annealing)

sumI=sum(I_r,3)/size(I_r,3);
opts.scale = ones(size(idx)); %scale: LED brightness map (but all galvo presumed to be same value)
opts.O0 = imresize(sqrt(sumI),nObj); %O0: initial guess for obj (in real space) %Initialize to normalized coherent sum
%opts.O0 = upsamp(opts.O0);
opts.P0 = w_NA; %P0: initial guess for P
opts.Ps = w_NA;
opts.H0 =defineProp(lambda, z0, k0, opts.Prop_mode,con);
%H0: known portion of the aberration function, 
% e.g. sample with a known defocus-induced quadratic aberration
% function can be defined here

opts.F = F; %F: operator of Fourier transform 
opts.Ft = Ft; %Ft: operator of inverse Fourier transform

%spatial coordiates for object space
opts.x_obj=x;
opts.y_obj=y;
opts.xObj=xObj;
opts.dx_obj=dx_obj;
opts.con=con; %Save the conversion factor

rStart=tic;
[objR,Pupil,errorR,scaleR,freqCent_cal] = ptychReconstruct(I,nObj,freqCent,opts);
reconTime=toc(rStart);

save([baseDir filename '_ptych.mat'],'objR','Pupil')
end

figure(); imagesc(abs(objR(500:650,500:650))); axis image; axis off; set(gcf,'color','w'); colorbar; colormap gray


























