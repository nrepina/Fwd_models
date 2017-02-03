%%
close all; clear all
%% Cell3D high res plot
res=0.05;
cell_highres= cell3D(res);
alpha_map=cell3D_alpha(res);
%alpha_map=gradient(cell_highres);
%max1=max(alpha_map(:))
%alpha_map_scaled=0.01*alpha_map./max1;


vol3d('cdata', cell_highres, 'texture', '3D', 'Alpha', alpha_map);
camlight(300,-80); lighting gouraud; light('Color', 'w'); axis equal;
view(3); 
colormap('summer');
colorbar;
set(gcf, 'Position', [100 100 800 800]);
set(gca, 'XLim', [0 180], 'YLim', [0 180], 'ZLim', [0 180]) ;
tickx = get(gca,'xtick')
ticky = get(gca,'ytick')
tickz = get(gca,'ztick')
set(gca,'ytick', tickx, 'ztick', tickx, 'xticklabel',(tickx.*res), 'yticklabel',(tickx.*res), 'zticklabel',(tickx.*res));
set(gca,'FontSize',20);
export_fig cell_highres.png -m1 -transparent

%% Cell3D thick plot
res=0.1;
cell_thick= cell3D_2(res);
alpha_map=cell3D_2_alpha(res);

vol3d('cdata', cell_thick, 'texture', '3D', 'Alpha', alpha_map);
camlight(300,-70); lighting gouraud; light('Color', 'w'); axis equal;
view(3); 
colormap('summer');
colorbar;
set(gcf, 'Position', [100 100 800 800]);
set(gca, 'XLim', [0 90], 'YLim', [0 90], 'ZLim', [0 390]) ;
tickx = get(gca,'xtick');
ticky = get(gca,'ytick');
tickz = get(gca,'ztick');
%set(gca,'xtick', ticky, 'xticklabel',(ticky.*res), 'yticklabel',(ticky.*res), 'zticklabel',(tickz.*res));
set(gca,'ytick', tickx, 'xticklabel',(tickx.*res), 'yticklabel',(tickx.*res), 'zticklabel',(tickz.*res));
set(gca,'FontSize',20);
export_fig cell_thick.png -m2 -transparent


%% Cell2D high res plot
res=0.1;
cell_highres= cell3D(res);
cell_highres_46=cell_highres(:,:,46);

%color = [0, 0.565, 0.482
    %0.569, 0.753, 0.478
    %0.996, 0.984, 0.545]

image(cell_highres_46, 'CDataMapping','scaled')
colormap('summer')
colorbar
axis equal
set(gcf, 'Position', [100 100 800 800]);
tickx = get(gca,'xtick');
ticky = get(gca,'ytick');
set(gca,'FontSize',20)
set(gca,'xtick', ticky, 'xticklabel',(ticky.*res), 'yticklabel',(ticky.*res))
export_fig cell_highres_2D46.png -m1 -transparent

%% Cell2D thick plot
res=0.1;
cell_thick= cell3D_2(res);
cell_thick_61=cell_thick(:,:,61);

%color = [0, 0.565, 0.482
    %0.569, 0.753, 0.478
    %0.996, 0.984, 0.545]

image(cell_thick_61, 'CDataMapping','scaled')
colormap('summer')
colorbar
axis equal
set(gcf, 'Position', [100 100 800 800]);
tickx = get(gca,'xtick');
ticky = get(gca,'ytick');
set(gca,'xtick', ticky, 'xticklabel',(ticky.*res), 'yticklabel',(ticky.*res))
set(gca,'FontSize',20)
export_fig cell_thick_2D61.png -m1 -transparent

%% FDTD
close all;
load('fdtd_33_Ex.mat');
cell= Ex_ampl(80:230,80:230,:);

max_col = max(cell,[],1);
max_cell = max(max_col,[],2);
cell_scaled= cell;

for i=1:size(cell, 1);
    for j=1:size(cell, 2);
        for k=1:size(cell, 3);
         cell_scaled(i,j,k)= cell_scaled(i,j,k)./max_cell(1,1,k);
         
         if cell_scaled(i,j,k) <0.83;
             cell_scaled(i,j,k)= 0.001;
         end
        end
    end
end


%alpha_map= cell./max(cell(:));
%alpha_scale=0.9*ones(size(alpha_map, 1), size(alpha_map, 2), size(alpha_map, 3));
%alpha_map_scaled= alpha_map - alpha_scale;

alpha_map_scaled= cell_scaled;
vol3d('cdata', cell_scaled, 'texture', '3D', 'Alpha', alpha_map_scaled);
camlight(300,-70); lighting gouraud; light('Color', 'w'); axis equal;
view(3); 
colormap('summer');
caxis([0.83, 1]);
colorbar;
set(gcf, 'Position', [100 100 800 800]);
%set(gca, 'XLim', [0 90], 'YLim', [0 90], 'ZLim', [0 390]) ;
tickx = get(gca,'xtick');
ticky = get(gca,'ytick');
tickz = get(gca,'ztick');
set(gca,'xticklabel',(tickx.*0.1), 'yticklabel',(ticky.*0.1), 'zticklabel',(ticky.*0.1));
%set(gca,'ytick', tickx, 'xticklabel',(tickx.*res), 'yticklabel',(tickx.*res), 'zticklabel',(tickz.*res));
set(gca,'FontSize',20);
export_fig cell_thick.png -m1 -transparent


%sliderDisplayIm(Ex_ampl(80:230,80:230,:));

