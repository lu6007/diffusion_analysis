% function estimate_simulation(cell_name)
% estimate the diffusion coefficient based on simulation result

% Copyright: Shaoying Lu and Yingxiao Wang 2012
function estimate_simulation(cell_name)
% Load the image and define the boundary to
% be a little smaller than the original boundary to
% avoid the edge effect during image processing.
data = diffusion_init_data(cell_name);
path = data.path;
dt = data.dt; display(dt);
boundary_file = strcat(path, 'boundary_small.data');
im = imread(strcat(path, 'output/1.tiff'));

% This part cannot be replaced by the get_polygon function
% because of the special treatment of the boundary dots.
if ~exist(boundary_file, 'file'),
    figure(2); imshow(im); caxis auto; axis xy;
    hold on; 
    tmp = load(strcat(path, data.boundary_file));
    boundary = tmp.boundary;
    if(iscell(boundary)),
        x1 = boundary{1}(:,1);
        y1 = boundary{1}(:,2);
    else
        x1 = boundary(:,1);
        y1 = boundary(:,2);
    end; clear boundary;

    plot(x1, y1, 'b.');
    title('Please mark cell boundary');
    [~,~, ~, x,y] = roipoly();
    boundary = [x,y];
    save(boundary_file,'boundary', '-ascii'); 
else
    boundary = load(boundary_file);
end;

% Estimate the diffusion constant
% Creat a mesh and Refine the mesh
cbound = [100, 2^16];
figure(1); imshow(im); caxis auto; hold on; axis xy;
set(gca, 'YDir','reverse','CLim',cbound);
title('Image overlaid with boundary');
plot(boundary(:,1), boundary(:,2),'r.');
mesh =create_mesh(boundary);
new_mesh = refine_mesh(mesh, 'num_refines', 3,'method', 1);
p_image = new_mesh.node; tri = new_mesh.tri; 
pdemesh(p_image,new_mesh.edge,new_mesh.tri);
clear mesh new_mesh;

% assemble matrices
p = scale_by_magnification(p_image, data.mag);
[K,M] = assemble_matrix(p,tri);

% read images as concentrations.
num_images = 2;
image_size = size(im); 
con_array = zeros([image_size,num_images]);
%x0 = [1 image_size(2)]; y0 = [1 image_size(1)];
%bw = roipoly(x0, y0, im, x1, y1); %figure; imagesc(bw);
nhood = [10,10];
for i=1:num_images,
    file_name = strcat(path, 'output/', num2str(i), '.tiff'); 
    image = imread(file_name);
    %temp = rgb2gray(image); % not needed for b/w images
    %con_array(:,:,i) = image;
    % smoothing is sometimes very important
    temp = image;
    %temp = image(:,:,1);
    con_array(:,:,i) = wiener2(temp,nhood); 
    %con_array(:,:,i) = region_wiener2(uint16(temp), nhood, bw);
    clear temp image;
end;

% convert concentration to vector
num_nodes = size(p_image,2);
display(num_nodes);
u = zeros(num_nodes,num_images);
p_ij = [p_image(2,:); p_image(1,:)];
for i= 1:num_images,
    temp = concentration_to_vector(con_array(:,:,i),p_ij);
    u(:,i) = double(temp); clear temp;
%     sum_Mu_i = sum(M*u(:,i));
%     display(sum_Mu_i);
    
    figure(4+i); pdesurf(p_image,tri,u(:,i));
    set(gca,'YDir','reverse','CLim',cbound);
    view(2); shading interp; colormap jet; colorbar;
    title(strcat('u\_',num2str(i)));
end;

% estimate diffusion coefficient
for i = 1:num_images-1,
    [diff_const, est_u, R] = ...
        get_diffusion_constant(u(:,i),u(:,i+1),M, K, dt);
    display(diff_const);
    display(R);
    Mdu = M*(u(:,i+1)-u(:,i));
    dtKu = -dt*K*u(:,i+1);
    %region = mark_region(Mdu,dtKu);
    relative_residual_norm = norm(u(:,i+1)-est_u,2)/...
        norm(u(:,i+1),2);
    display(relative_residual_norm);
    figure(10); plot(Mdu,dtKu,'+');
    set(gca, 'FontSize', 16);
    xlabel('Mdu'); ylabel('dtKu');
    title('Scattered plot of linear fitting');
%     figure; trimesh(tri(1:3,:)',p_image(1,:)',p_image(2,:)',region); 
%     view(2); colormap jet; shading interp;
    h = figure(11); pdesurf(p_image, tri,u(:,i+1)); 
    set(gca,'YDir','reverse','CLim',cbound);
    view(2); colormap jet; colorbar; 
    set(gca, 'FontSize', 16);
    title(strcat('u_', num2str(i+1)));
    saveas(h, strcat(path,'u_',num2str(i+1),'.png'));
    h = figure(12); pdesurf(p_image,tri,est_u);
    set(gca,'YDir','reverse','CLim',cbound);
    view(2); colormap jet; colorbar; 
    set(gca, 'FontSize', 16);
    title(strcat('est u_',num2str(i+1)));
    saveas(h, strcat(path,'est_u_',num2str(i+1),'.png'));
    h = figure(13); pdesurf(p_image, tri,est_u-u(:,i+1)); 
    set(gca,'YDir','reverse','CLim',[-1000 1000]);
    view(2); colormap jet; colorbar; 
    set(gca, 'FontSize', 16); 
    title(strcat('est u_',num2str(i+1),'-u_', num2str(i+1)));
    saveas(h, strcat(path,'cest_u_',num2str(i+1),...
        '_u_',num2str(i+1),'.png'));

end;
return;


% function region = mark_region(Mdu, dtKu)
% one_over_D = Mdu\dtKu;
% n = size(Mdu,1);
% region = zeros(size(Mdu));
% % shift = 0.3;
% for i = 1:n,
%     region(i) = dtKu(i)-one_over_D*Mdu(i);
% %     if abs(dtKu(i)-one_over_D*Mdu(i))>shift,
% %         region(i) = dtKu(i)-one_over_D*Mdu(i);
% %     end;
% end;
% return;

