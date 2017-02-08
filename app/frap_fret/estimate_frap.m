% Estimate the diffusion coefficient based on 
% FRAP image results.
% function estimate_frap(cell_name,starting_step)
%
% Example
% >> cell_name = 'mem17';
% >> data = sample_diffusion_init_data(cell_name);
% >> estimate_frap(data);
%
% Sub-functions
% function [h1, h2] = plot_Mdu_dtKu_residual(diff_const, ...
%    u1, u2, M, K, dt)
% function is_boundary = mark_mesh_boundary(p_image, edge)

% This application was used to estimate the diffusion coefficient of
% LeiLei's TRPC6 paper.
% The diffusion analysis procedure resembles that described in our previous 
% publication [ref: Lu S Wang Y 2008 PLoS Computational Biology].  
% Specifically, the intensity images were normalized by dividing 
% the median of five images before photobleaching to calculate 
% the ratio images which represent the relative concentration maps of 
% the fluorescence proteins. The concentration maps were then filtered 
% by a 100x100 (15 µm x 15µm) local adaptive filter restricted within 
% a cell body detected by segmentation [ref: Lu S Wang Y 2008 PLoS 
% Computational Biology]. As a result, the concentration maps before 
% photobleach had values close to 1 at all pixels within the cell, 
% while those after photobleach had values <1 near the photobleached 
% region. Seven filtered concentration maps right after photobleach were 
% used to estimate 6 values of diffusion coefficients. Since the 
% fluorescence proteins are expected to diffuse during photobleach, 
% a mobile subcellular region can be determined by selecting the 
% subcellular region where the concentration map right after photobleach 
% had values <0.98. Therefore, within the selected mobile subcellular 
% region, the movement of the fluorescence proteins can be modeled by 
% Fick’s Law of diffusion [ref: Lu S Wang Y 2008 PLoS Computational 
% Biology]. The concentration maps within the shape of the cells were 
% discretized using the finite element method with a triangular mesh, 
% to obtain a linear model of “the weighted discrete Laplacian of 
% concentration” and “the weighted change of concentration in time” 
% [ref: Lu S Wang Y 2008 PLoS Computational Biology]. By linear regression, 
% a diffusion coefficient D can be estimated with a residual r and the 
% correlation coefficient R which measures the quality of regression. 
% A value of R no less 0.7 (or R>=0.5 and r<=twice of the maximal residual 
% when R>=0.7) is considered a good fit. The corresponding diffusion 
% coefficients, D, were then entered in the statistical comparison. The 
% mean values of diffusion coefficients from different groups were 
% compared with one-tailed t-test of unequal variance. P-value<=0.05 was 
% considered significant difference. The fluorescence intensity images 
% were all background subtracted, filtered by a 3x3 median filter and 
% cropped before processing. It was also confirmed that the cell shape 
% did not change significantly within the 110 sec when the images before 
% or after photobleach were collected.

% Copyright: Shaoying Lu and Yingxiao Wang 2012
function estimate_frap(data, varargin)
parameter_name = {'remove_result','type'};
default_value = {0,'cell'};
[remove_result, type] = parse_parameter(parameter_name, default_value, varargin);

% data = diffusion_init_data(cell_name);
% if ~exist(data.output_path, 'dir'),
%     mkdir(data.output_path);
% end;
result_file = strcat(data.output_path, 'result.mat');
if remove_result,
    delete(result_file);
end;

mag = data.magnification;
nhood = [min(mag,100), min(mag,100)];
first_file = strcat(data.path, data.first_file);
image_index = data.image_index;
index_before = data.index_before;
index_after = data.index_after;
pattern = data.index_pattern;
correct_photobleaching = 0;
if isfield(data,'time_format'),
    time_format = data.time_format;
else
    time_format = 1;
end;

% Load the images, preprocess and filter the images
last_frame_before = max(index_before); % last frame before pdgf    
num_images = length(image_index);
if ~exist(result_file,'file'),
    % initialize the cropping rectangle and background
    im = imread(first_file);
    [temp,data ] = preprocess(im, data); 
    im = uint16(zeros([size(temp),1, num_images])); clear temp;
    bd = cell(num_images, 1);
    bw = cell(num_images, 1);
    total_int = zeros(num_images, 1);
    if time_format ==1,
        time = get_time(data, 'ref_time', data.photobleach_time/100/60);
    elseif time_format ==2,
        time = zeros(length(image_index),1);
    end;
    for i = 1:num_images;
        index_str = sprintf(pattern{2}, image_index(i));
        image_file = regexprep(first_file, pattern{1}, index_str);
        if time_format ==2,
            time(i) = get_time_2(image_file);
        end;
        temp = imread(image_file);
        im_t = preprocess(temp,data);
        % figure; imagesc(im{i}); hold on;
        min_area = 1000;
        switch type
            case 'cell',
                if isfield(data,'threshold') && data.threshold,
                    [bd_t, ~] = get_cell_edge(im_t, 'method', 1, 'threshold', data.threshold, ...
                        'smoothing_factor', 3, 'min_area', min_area, 'show_figure', 1);
                elseif isfield(data,'brightness_factor') && data.brightness_factor,
                    figure; imagesc(im_t);
                    [bd_t,~] = get_cell_edge(im_t, 'method',1, ...
                        'smoothing_factor', 3,'show_figure',1, 'brightness_factor', data.brightness_factor,...
                        'min_area', min_area, 'show_figure',1);
                else % there is no threshold
                    figure; imagesc(im_t);
                     [bd_t, ~] = get_cell_edge(im_t, 'method', 1, ...
                        'smoothing_factor', 3,'min_area', min_area,  'show_figure', 1);
                end;
            case 'solution',
                if i ==1,
                    figure; imagesc(im_t); title('Please Select the Region for the Solution');
                    [~, xi, yi] = roipoly();
                    %hold on; plot(yi,xi, 'r');
                    bd_t = [yi xi];
                else % i>1
                    % bw_t = bw{1};
                    bd_t = bd{1};
                end
        end
                
        im(:,:,1,i) = im_t;
        if size(bd_t,1) ==1, %Compatibility
            bd{i} = bd_t{1};
        else
            bd{i} = bd_t;
        end;
        bw_t = true(size(im_t));
        % The mask was calculated based on the boundary to 
        % avoid the problem that the detected mask does not
        % give a connected mask image.
        bw{i} = bd2im(bw_t, bd{i}(:,2), bd{i}(:,1)); 
        clear index_str image_file temp bd_t bw_t im_t;
    end;

    % filter the cropped images
    % convert images to concentration map
    % calculate the total intensity for photobleaching correction
    image_array = zeros(size(im));
    for i = 1:num_images,
        image_array(:,:,1,i) = region_wiener2(im(:,:,1,i),nhood,...
            bw{last_frame_before});
        total_int(i) = sum(sum(double(bw{i}).*double(im(:,:,1,i))));
    end;
    dimension = 4;
    image_0 = median(im(:,:,1,index_before), dimension);
%image_0 = im(:,:,1,last_frame_before);
    %
    con = zeros(size(image_array));
    for i = 1:num_images,
        temp = compute_ratio(im(:,:,1,i), image_0,'shift', 1.0e-4);
        con(:,:,1,i) = region_wiener2(temp,nhood,bw{last_frame_before}); 
        clear temp;
    end;
    save(result_file,'im', 'bd', 'image_0', 'image_array', 'con','total_int', 'time');
    clear im_array image_array;
else %if ~exist(result_file,'file'),
    res = load(result_file);
    con = res.con;
    %image_0 = res.image_0;
    im = res.im;
    bd = res.bd;
    total_int = res.total_int;
    if isfield(res, 'time'),
        time = res.time;
    else % for backward compatibility
        time = get_time(data, 'ref_time', data.photobleach_time/100/60);
        save(result_file, 'im', 'con', 'bd', 'time');
    end
end; %if ~exist(result_file,'file'),

% Confirm that the boundaries are ok.
color = ['r' 'w' 'g' 'y' 'k' 'b' 'g'];
ll = last_frame_before;
figure; imagesc(im(:,:,1,ll)); hold on;
title('Last image before PB with its (r) and other edges (w,g, etc) before PB');
for i = 0:length(index_before)-1,
    plot(bd{ll-i}(:,2), bd{ll-i}(:,1), color(i+1), 'LineWidth',2); 
end;
figure; imagesc(im(:,:,1,ll)); hold on;
title('Last image before PB with its (r) and other edges (w,g, etc) after PB');
for i = 0:length(index_after)-1,
   plot(bd{ll+i}(:,2), bd{ll+i}(:,1), color(i+1), 'LineWidth',2); 
end;
if ~correct_photobleaching,
    pb_factor = 1;
else
    % correct for photobleaching during imaging
    % pb_factor estimates photobleach
    % here we assume that all the frames before photobleaching are connected
    % and each frame had similar amount of photobleaching during imaging.
    pb_factor = mean(total_int(2:ll)./total_int(1:ll-1));
end;

%
% figure; imagesc(image_0); colorbar; title('image\_0');
% for i = last_frame_before: last_frame_before+4,
%     figure; imagesc(con(:,:,1, i)); 
%     caxis([0.5 1.1]); colorbar;
%     title(strcat('Concentration Map ', num2str(image_index(i))));
% end;
% 

%% estimate diffusion coefficients.
% mark a more restricted boundary by shrinking 10 pixels
% or, alternatively, manually select the mobile region
first_frame_after = min(index_after); % first frame after pdgf
temp = con(:,:,1, first_frame_after);

% num_pixels_erode = 4;
% bw_t = roipoly(temp, bd{ll}(:,2), bd{ll}(:,1));
% [num_rows num_cols] = size(bw_t);
% bw_t(1,:) = 0; bw_t(num_rows, :) = 0;
% bw_t(:,1) = 0; bw_t(:, num_cols) = 0;
% struct_ele = strel('disk',num_pixels_erode);
% bw_small = imerode(bw_t, struct_ele);
% figure; imshow(bw_small); hold on; plot(bd{ll}(:,2), bd{ll}(:,1), 'r');
% bd_t = get_cell_edge(bw_small, 'method', 1, 'smoothing_factor', 3,...
% 'show_figure',1);
% num_refines = 2;
% smaller_boundary = bd_t{1};

boundary_file = strcat(data.output_path,'smaller_boundary.mat');
title_str = 'Please selected smaller boundary';
[~, bd_t]= get_polygon(temp, boundary_file,title_str);
% swap two rows in smaller_boundary in concentration_to_vector
smaller_boundary = bd_t{1}; %[bd_t{1}(:,2) bd_t{1}(:,1)];
display('Got the polygon, and continue...');
clear bw_t bd_t temp;

%
cbound = [0.5 1.1];
num_refines = 2;
method = 2; 

% create mesh and assemble the matrices
mesh = create_mesh(smaller_boundary);
figure; imagesc(con(:,:,1, first_frame_after));
hold on; pdemesh(mesh.node,mesh.edge,mesh.tri); colorbar;
new_mesh = refine_mesh(mesh, 'num_refines', num_refines, 'method',1);
clear mesh; mesh = new_mesh; clear new_mesh;
p_image = mesh.node; edge = mesh.edge; tri = mesh.tri;
clear mesh;
p = scale_by_magnification(p_image,mag);
[K,M] = assemble_matrix(p,tri);
%is_boundary = mark_mesh_boundary(p_image, edge);
figure; imagesc(con(:,:,1, first_frame_after));
hold on; pdemesh(p_image,edge,tri); colorbar;

% convert concentration to solution vector
num_nodes = size(p_image,2);
display(sprintf('num_nodes = %d', num_nodes));
num_steps = length(index_after);
u = zeros(num_nodes,num_steps);
est_u = zeros(num_nodes,num_steps);
r = zeros(num_nodes, num_steps);
for i = 1:num_steps,
    temp = con(:,:,:,first_frame_after+i-1);
    u(:,i) = concentration_to_vector(temp,p_image, 'method', method);
%     figure; imagesc(temp); title(strcat('Concentration map #', num2str(i))); 
%     colorbar; hold on;
% %     pdemesh(p_image,edge,tri, u(:,i)'); 
%     figure; pdesurf([p_image(2,:); p_image(1,:)],tri,u(:,i));
%     view(90,90); colormap jet; caxis([0 1.1]); colorbar;
%     title(strcat('solution #', num2str(i))); 
    clear temp;
end;


% estimate diffusion constant
display(sprintf('%30s%20s%20s%20s',...
    'diffusion coef (um^2/sec)', 'residual_1', 'residual_2','R value'));
diff_const = zeros(num_steps-1,1);
for i= 1:num_steps-1,
    % dt in the unit of minutes
    dt = time(first_frame_after+i)-time(first_frame_after+i-1);
    u_i = u(:,i)*pb_factor;
    Mdu = M*(u(:,i+1)-u_i);
    dtKu = -dt*K*0.5*(u(:,i+1)+u_i);
    [diff_const(i) , est_u(:,i+1) , R] = ...
        get_diffusion_constant(u_i,u(:,i+1),M,K,dt);

    relative_residual_norm = norm(u(:,i+1)-est_u(:,i+1),2)/...
        norm(u(:,i+1),2);
    relative_residual_norm_2 = norm(u(:,i+1)-est_u(:,i+1),2)/...
        norm(u(:,i+1)-u_i);
    display(sprintf('%34.4f%30.4f%20.4f%20.4f', diff_const(i)/60, ...
        relative_residual_norm, relative_residual_norm_2, R));
end;

show_figure_result = 1;
if show_figure_result,
    %cbound = [0.5, 1.1];
    s = max(p_image)-min(p_image);
    for i = 1:2,
        h = figure; 
            pdesurf(p_image,tri,u(:,i));
        view(90,90); caxis(cbound); colorbar;
        set(gca, 'FontSize', 16);
        title(strcat('Concentration map ', num2str(i))); 
        colormap jet;
%         saveas(h,sprintf('%su_%d.fig',path,i+starting_step));
%         saveas(h,sprintf('%su_%d.png',path,i+starting_step));
    end;
    for i =2:2,
        h = figure;
        pdesurf(p_image,tri,est_u(:,i));
        view(90,90); colormap jet; caxis(cbound); colorbar;
        set(gca, 'FontSize', 16);
        title('Estimated concentration map'); 
%         saveas(h,sprintf('%sest_u_%d.fig',path,i+starting_step));
%         saveas(h,sprintf('%sest_u_%d.png',path,i+starting_step));
    end;

    error = zeros(size(u));
    for i = 1:1,
       %  error(:,i+1) = (u(:,i+1)-est_u(:,i+1))/abs(diff_const(i))/dt;
        error(:,i+1) = (u(:,i+1)-est_u(:,i+1))/norm(u(:,i+1)-u(:,i)*pb_factor);
    end;
%     for i = 2:3,
%         [R,P]=corrcoef(error(:,i), error(:,i+1))
%     end;

    cbound = [-0.005, 0.02];
    for i = 1:1,
        h = figure; 
        pdesurf(p_image,tri, abs(error(:,i+1)));
        view(90,90); colormap jet; caxis(cbound);colorbar;
        set(gca, 'FontSize', 16);
        title('abs(u - est u) '); 
%         saveas(h,sprintf('%sabs_diff_u_%d_est_u_%d.fig',path,i+1,i+1));
%         saveas(h,sprintf('%sabs_diff_u_%d_est_u_%d.png',path,i+1,i+1));
     end;

    for i = 1:1,
        [h1, h2] = plot_Mdu_dtKu_residual(diff_const(i), u(:,i)*pb_factor, u(:,i+1), M, K, dt);
%         saveas(h1, sprintf('%sMdu_dtKu_%d.fig', path, i+starting_step));
%         saveas(h1, sprintf('%sMdu_dtKu_%d.png', path, i+starting_step));
%         saveas(h2, sprintf('%sresidual_%d.fig', path, i+starting_step));
%         saveas(h2, sprintf('%sresidual_%d.png', path, i+starting_step));
    end;

end; % if show_figure_results

return;
function [h1, h2] = plot_Mdu_dtKu_residual(diff_const, ...
    u1, u2, M, K, dt)
    alpha = 0.5; 
    beta = 0.5;
    Mdu = M*(u2-u1);
    dtKu = -dt*K*(alpha*u1+beta*u2);
    %diff_const = diff_const/2;
    r0 = dtKu-Mdu/diff_const;
    h1 = figure;
    % r = r0.*is_boundary;
    % plot( dtKu+r/diff_const,Mdu, '+'); hold on; box off;
    plot(Mdu, dtKu, '+'); hold on; box off;
    set(gca, 'FontSize', 16);
    set(gca, 'LineWidth',2);
    xbound = [-4e-3 6e-3];
    set(gca, 'xlim', xbound, 'xtick', [-4e-3, 0, 4e-3]);
    % mdu1 = min(Mdu); mdu2 = max(Mdu); mdu12 = [mdu1, mdu2];
    % plot(mdu12, 1./diff_const*mdu12,'r-','LineWidth',3);
    plot(xbound , 1./diff_const*xbound, 'r-', 'LineWidth',3);
    title('dtKu vs. Mdu'); xlabel('Mdu'); ylabel('dtKu');
    %legend('Mdu', 'Linear Fitting', 'Location','SouthEast');
    %legend('boxoff');

    h2 = figure; 
    plot(Mdu, r0,'+'); hold on; box off;
    set(gca, 'FontSize', 16);
    set(gca, 'LineWidth', 2);
    xbound = [-4e-3 6e-3];
    set(gca, 'xlim', xbound, 'xtick', [-4e-3, 0, 4e-3]);
    plot(xbound , 0.*xbound, 'r-', 'LineWidth',3);
    title('residual vs. Mdu'); xlabel('Mdu'); ylabel('residual');
return;


function is_boundary = mark_mesh_boundary(p_image, edge)
    num_nodes = size(p_image,2);
    num_edges = size(edge,2);
    is_boundary = zeros(num_nodes,1);
    for i = 1: num_edges,
        is_boundary(edge(1,i)) = 1;
        is_boundary(edge(2,i)) = 1;
    end;
return;

