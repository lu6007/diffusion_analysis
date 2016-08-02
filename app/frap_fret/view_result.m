% function view_result(data)
%
% Example:
% >> cell_name = 'photobleach_cell';
% >> data = sample_diffusion_init_data(cell_name);
% >> view_result(data);
%

% Copyright: Shaoying Lu and Yingxiao Wang 2013-2016
function view_result(data)
    pa = data.path;
    result = load(strcat(pa,'output/result.mat'));
    tri = result.tri; p_image = result.p_image;
    u = result.u;
    u_act = result.u_act;
    tri = tri(1:3,:)'; x = p_image(1,:); y = p_image(2,:);
    cbound = [0.0,0.85];
    xylimit= [100, 240, 100,200];
    load(strcat(pa, 'cmap_0208_2007.mat'));

    % % make movie
    % h = figure('Position', [1,150,1000, 500]);
    % for i = 1:num_frames,
    %     subplot(1,2,1, 'replace');
    %     trisurf(tri,x,y,u(:,i),'LineStyle', 'none'); 
    %     axis(xylimit);
    %     view(2); grid off; shading interp; caxis(cbound); 
    %     colormap(cmap); colorbar;
    %     title('FRET Image');
    %     subplot(1,2,2, 'replace');
    %     trisurf(tri,x,y,u_act(:,i),'LineStyle', 'none'); 
    %     axis(xylimit);
    %     view(2); grid off; shading interp; caxis(cbound); 
    %     colormap(cmap); colorbar;
    %     title('Subtract Diffusion');
    %     movie(i) = getframe(h);
    % end;
    % save movie.mat movie;
    % movie2avi(movie, 'lyn-src.avi', 'Compression', 'Cinepak', 'fps',3);



    % egf_pp1_lyn2: first = 12, last = 24
    % draw 3 images
    %first = data.image_index(1);
    %index = [2, 7, 12];
    ii = 1;
    for i = data.image_index,
        h = figure;
        figure(h); trisurf(tri,x,y,u(:,ii),'LineStyle','none');
        axis(xylimit); 
        view(2); grid off; shading interp; caxis(cbound); colormap(cmap); colorbar;
        saveas(h, sprintf('%soutput/u.%03d.jpg',pa,i));
        h = figure;
        figure(h); trisurf(tri,x,y,u_act(:,ii),'LineStyle','none');
        axis(xylimit); 
        view(2); grid off; shading interp; caxis(cbound); colormap(cmap); colorbar;
        saveas(h, sprintf('%soutput/u_act.%03d.jpg',pa, i));
        ii = ii+1;
     end;




    % % save grayscale images
    % % cbound = [-1.0, 4.0];
    % cbound = [0.0,2.0];
    % h= figure;
    % for i = 1:num_frames,
    %       figure(h); trisurf(tri,x,y,u(:,i),'LineStyle','none');
    %       view(2); grid off; caxis(cbound); colorbar; colormap gray;
    %       saveas(h, sprintf('u.%03d.tif',first+i-1), 'tif');
    %       figure(h); trisurf(tri,x,y,u_act(:,i),'LineStyle','none');
    %       view(2); grid off; caxis(cbound); colorbar; colormap gray;
    %       saveas(h, sprintf('u_act.%03d.tif',first+i-1),'tif');
    % end;

    % % define regions using roiply() function
    % % [x,y,regions(1).bw,regions(1).x, regions(1).y] = roipoly;
    % load regions.mat;
    % num_regions = size(regions,2);
    % col = ['r' 'g' 'b' 'y'];
    % % 
    % % plot regions
    % im = imread('u_act.010.png');
    % figure; imshow(im);
    % hold on;
    % for i = 1:num_regions;
    %     plot(regions(i).x, regions(i).y,col(i),'LineWidth',3);
    % end;
    % 
    % 
    % % compute average intensities
    % u_ai = zeros(num_frames,num_regions);
    % u_act_ai = zeros(num_frames, num_regions);
    % for i = 1:num_frames,
    %     im = imread(sprintf('u.%03d.tif',first+i-1));
    %     matrix = double(im(:,:,1));
    %     for j = 1:num_regions,
    %         bw = regions(j).bw;
    %         u_ai(i,j) = sum(sum(matrix.*bw))/sum(sum(bw));
    %     end;
    %     im = imread(sprintf('u_act.%03d.tif', first+i-1));
    %     matrix = double(im(:,:,1));
    %     for j = 1:num_regions,
    %         bw = regions(j).bw;
    %         u_act_ai(i,j) = sum(sum(matrix.*bw))/sum(sum(bw));
    %     end;
    % end;
    % save u_ai.mat u_ai;
    % save u_act_ai.mat u_act_ai;
    % 
    % time = [1:num_frames]*60;
    % load u_ai.mat; load u_act_ai.mat;
    % figure; hold on;
    % for j = 1:num_regions,
    %     plot(time, u_ai(:,j),col(j),'LineWidth',3);
    % end;
    % figure; hold on;
    % for j = 1:num_regions,
    %     plot(time, u_act_ai(:,j),col(j),'LineWidth',3);
    % end;
    % 

return

          
