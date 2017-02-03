oldSz=20;
objSz=size(UDS);
NsampR=round(oldSz/(dpix_c/mag)).*ones(1,2);
cent=floor(objSz/2)+1;
downsamp= @(x) x(cent(1)-floor(NsampR(1)/2):cent(1)+floor(NsampR(1)/2)-1, ...
    cent(2)-floor(NsampR(2)/2):cent(2)+floor(NsampR(2)/2)-1,:);

F = @(x) fftshift(fft2(ifftshift(x)));
Ft = @(x) fftshift(ifft2(ifftshift(x)));

UDS_d=Ft(downsamp(F(UDS)));


newSz=100;
NsampR=round(newSz/(dpix_c/mag)).*ones(1,2); %Number of pixels on camera to capture the total size
Ilarge=zeros([NsampR numIllum]);
oldN=size(UDS_d);
newAmt=round((NsampR-oldN(1:2))/2);

Ilarge(newAmt:end-newAmt-1,newAmt:end-newAmt-1,:)=abs(UDS_d).^2;
