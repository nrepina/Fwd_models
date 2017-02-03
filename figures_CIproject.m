if 0
filename='lowResCells.gif';

figH=figure(); axis image; caxis([1.31 1.5]);
axis off; set(figH,'color','w'); hold on;

for ii=1:size(object,3)
    figure(figH);
    imagesc(object(:,:,ii))
    drawnow
    im=frame2im(getframe(figH));
    im=imresize(im(64:748,239:922,:),0.3);
    [iInd,CM]=rgb2ind(im,256);
    if ii==1
        imwrite(iInd,CM,filename,'gif','LoopCount',inf,'DelayTime',1/30);
    else
        imwrite(iInd,CM,filename,'gif','WriteMode','append','DelayTime',1/30);
    end
    
end
end

%%
if 0
filename='fdtd_2_results.gif';

figH=figure(); axis image;
axis off; set(figH,'color','w'); hold on;

for ii=1:size(I,3)
    figure(figH);
    imagesc(I(:,:,ii))
    drawnow
    im=frame2im(getframe(figH));
    im=im(303:503,479:679,:);
    [iInd,CM]=rgb2ind(im,256);
    if ii==1
        imwrite(iInd,CM,filename,'gif','LoopCount',inf,'DelayTime',1/3);
    else
        imwrite(iInd,CM,filename,'gif','WriteMode','append','DelayTime',1/3);
    end
    
end
end
%%
if 0
filename='fdtd_stack_2.gif';

figH=figure(); axis image; caxis([0 1])
axis off; set(figH,'color','w'); hold on; 

for ii=1:size(temp,3)
    figure(figH);
    imagesc(flipud(temp(:,:,ii)))
    drawnow
    im=frame2im(getframe(figH));
    %im=imresize(im(64:748,239:922,:),0.3);
    im=im(233:579,147:1014,:);
    [iInd,CM]=rgb2ind(im,256);
    if ii==1
        imwrite(iInd,CM,filename,'gif','LoopCount',inf,'DelayTime',1/10);
    else
        imwrite(iInd,CM,filename,'gif','WriteMode','append','DelayTime',1/10);
    end
    
end
end
%%
if 0
filename='born_single.gif';

figH=figure(); axis image; caxis([1 1.1])
axis off; set(figH,'color','w'); hold on; 

for ii=1:size(temp2,3)
    figure(figH);
    imagesc(temp2(:,:,ii))
    drawnow
    im=frame2im(getframe(figH));
    im=imresize(im(64:748,239:922,:),0.3);
    %im=im(303:503,479:679,:);
    [iInd,CM]=rgb2ind(im,256);
    if ii==1
        imwrite(iInd,CM,filename,'gif','LoopCount',inf,'DelayTime',1/30);
    else
        imwrite(iInd,CM,filename,'gif','WriteMode','append','DelayTime',1/30);
    end
    
end
end

%%
if 1
    temp=I(:,:,25);
[n,nbin]=hist(temp(:),100);
[~,mI]=max(n);
figH=figure(); imagesc(temp./nbin(mI)); axis image; colormap gray; colorbar; set(figH,'color','w')
caxis([0 2])
    
    
end

