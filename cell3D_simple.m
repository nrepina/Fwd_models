function [object,cellsize_z,matrix_res]=cell3D_simple(matrix_res);
%parameters for cell size and resolution between object matrix elements
cellsize_x= 20; %um
cellsize_y= 20; %um
cellsize_z= 8; %um
%matrix_res= 0.05; %um

%initialize object matrix
%NA=0.65
%res=0.4um (0.532/2*NA)
%DoF=1.25um (0.532/NA^2)
x=0:matrix_res:(cellsize_x-matrix_res);
y=0:matrix_res:(cellsize_y-matrix_res);
z=0:matrix_res:(cellsize_z-matrix_res);
[X, Y, Z]=meshgrid(x,y,z);

object=1.33.*ones(length(x), length(y), length(z)); %object as refractive index; initialize in water, n=1.33
% %cell parameters
% r_cell_outer=5; %um
% r_cell_inner=4.5; %um %r_cell_outer - r_cell_inner is thickness of plasma membrane
% 
% %specifies location and refractive index of all spherical objects
% for i=1:length(x)
%     for j=1:length(y)
%         for k=1:length(z)
%             if (i-(length(x)/2))^2 + (j-(length(y)/2))^2 + (k-(length(z)/2))^2<(r_cell_outer/matrix_res)^2
%                 if (i-(length(x)/2))^2 + (j-(length(y)/2))^2 + (k-(length(z)/2))^2>(r_cell_inner/matrix_res)^2
%                     object(i,j,k)= 1.5; %refractive index of plasma membrane
%                 else
%                     object(i,j,k)= 1.3; %refractive index of cytoplasm
%                 end
%             end
%         end
%     end
% end
% figure;
% graph1= slice(X,Y,Z,object, [], 10, 5);
%     colormap(gray)
%     caxis([1,1.5])
%     set(graph1, 'EdgeAlpha',0.1)
%     colorbar
% %%

%change to ellipse

rad=[8,7.9,... %cell outer radius; cell inner radius (um)
    3.5,3.4,... %nucleus outer radius; nucleus outer radius
    0.4,... %organelle #1 radius
    0.4,... %organelle #2 radius
    0.4,... %organelle #3 radius
    0.4]; %organelle #4 radius

nr=[1.5,1.35,... %cell membrane refractive index (n) ; cytoplasm n
    1.5,1.35,... %nuclear membrane refractive index (n) ; inner nucleus n
    1.5,... %organelle #1 n, res comparison
    1.5,... %organelle #2 n, x res comparison
    1.5,... %organelle #3 n, y res comparison
    1.5]; %organelle #4 n, DoF

%units in um
cent=[mean(x),mean(y),mean(z);...
    mean(x),mean(y),mean(z);...
    
    mean(x)-0.8,mean(y)+.8,mean(z);...  
    mean(x)-0.8,mean(y)+.8,mean(z);...
    
    %organelles separated by 1.6um
    mean(x)+5,mean(y),mean(z);...
    mean(x)+6.6,mean(y),mean(z);...    
    mean(x)+5,mean(y)+1.6,mean(z);...
    
    %DoF 1.25um
    mean(x),mean(y)-5,mean(z)+2.5];

for ii=1:length(rad)
    object(sqrt((X-cent(ii,1)).^2 + (Y-cent(ii,2)).^2 + (Z-cent(ii,3)).^2)<rad(ii))=nr(ii);
end


%sliderDisplayIm(object)
end








