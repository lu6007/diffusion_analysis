% function create_concentration_map(cell_name,starting_step)

% Copyright: Shaoying Lu and Yingxiao Wang 2012
function create_concentration_map(cell_name,starting_step)

%[path, prefix, start_0, end_0,shift, ~] = set_parameter(cell_name);
    data = set_parameter(cell_name);
    path = data.path;
    prefix = data.prefix;
    start_0 = data.start_0;
    end_0 = data.end_0;
    shift = data.shift;

    num_steps = 4;
    num_images_per_layer =1;
    num_images = num_steps *num_images_per_layer;
    dimension = 4;
    shift = shift+starting_step*num_images_per_layer;
    nhood = [40,40];
    save_results = 1;

    %mark image_rectangle
    file_name = sprintf('%s%s%03d',path,prefix,end_0);
    im = imread(file_name);
    rectangle_file = strcat(path,'image_rectangle.data');
    image_rectangle = get_rectangle(im, rectangle_file,'format', '-ascii');

    % mark boundary
    image_cropped = imcrop(im, image_rectangle);
    %if need_mark_boundary,
    boundary_file = strcat(path, 'boundary.data');
    %boundary = mark_boundary(image_cropped, boundary_file);
    % updated 07/06/2014
    % after running this, need to change the file name from
    % 'boundary.data' to 'boundar.mat' and change 'format' to '-ascii'
    [~, temp] = get_polygon(image_cropped, boundary_file, 'format', '-ascii');
    boundary = temp{1};


    % load boundary data
    bw = roipoly(image_cropped, boundary(:,1), boundary(:,2));
    clear image_cropped;
    % read and crop the image files
    image_0_array = ...
        get_cropped_image_array(path,prefix,start_0,end_0);
    for i = 1:end_0-start_0+1
        image_0_array(:,:,1,i) = ...
            region_wiener2(image_0_array(:,:,1,i),nhood,bw);
    end
    image_0_filtered = median(image_0_array,dimension);
    clear image_0_array ;
    image_0_filtered = image_0_filtered + im2uint16(~bw);
    image_array = ...
        get_cropped_image_array(path,prefix,end_0+shift,...
        end_0+shift+num_images-1);
    for i = 1:num_images
        image_array(:,:,1,i) = ...
            region_wiener2(image_array(:,:,1,i),nhood,bw);
    end

    % convert images to concentration map
    con = zeros(size(image_array));
    for i = 1:num_images
        con(:,:,:,i) = compute_ratio(image_array(:,:,:,i), image_0_filtered, ...
            'shift',0);
    end
    %clear image_array;
    for i = 1:num_images
        con(:,:,1,i) = region_wiener2(con(:,:,1,i),nhood,bw);
    end

    % take median in time dimension
    num_steps = floor(num_images/num_images_per_layer);
    display(num_steps);
    if num_steps <2
        error('num_steps must be no less then 2');
    end
    size_image= size(con);
    med_con = zeros([size_image(1:3),num_steps]);
    for i = 1:num_steps
        the_start = (i-1)*num_images_per_layer+1;
        the_end = i*num_images_per_layer;
        med_con(:,:,:,i) = median(con(:,:,:,the_start:the_end),dimension);
    end

    if save_results
        file_name = sprintf('%sresults_%d.mat',path,starting_step);
        save(file_name, 'image_0_filtered', 'image_array','med_con','con');
    end
return;


