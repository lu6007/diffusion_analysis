% function data = sample_diffusion_init_data(cell_name)
% Initialize cell data

% Copyright: Shaoying Lu and Yingxiao Wang 2014
function data = sample_diffusion_init_data(cell_name)
    % root = 'D:/sof/data/diffusion_sample/';
    root = '/Users/Yiwen/Desktop/for_yiwen/diffusion_data/diffusion/';
    data.cell_name = cell_name;
    switch cell_name
        case 'test',
            data.path = strcat(root,'simulation/test/');
            data.diff_const = [5.0 20.0 30.0 1.0]; % mu m^2/s from outside to inside
            data.first_file = 'cell_after_photobleach.PNG';
            image_before = imread(strcat(data.path, 'cell_before_photobleach.PNG'));
            rectangle_file = strcat(data.path, 'image_rectangle.data');
            data.rectangle = get_rectangle(image_before, rectangle_file, 'format', '-ascii');
            temp = imread(strcat(data.path,'cell_after_photobleach.PNG'));
            data.image_0 = imcrop(temp, data.rectangle); clear temp;
            data.num_bits = 8;
            data.mag = 100; %magnification
            % dt is not needed for computer simulation,
            % but it is needed for compute recovery curve and
            % estimate the diffusion coefficient from the simulation
            data.dt = 0.25;
        case 'mem17', % the best Lyn-Src cell
            data.path = strcat(root, 'frap/mem17/');
            data.output_path = strcat(data.path, 'output/');
            data.first_file = 'mem171.016';
            data.index_pattern = {'016', '%03d'};
            data.image_index = (26:36);
            data.index_before = 1:5; % Before 26-30, After 31-36
            data.index_after = 6:11;
            data.photobleach_time = 7287;
            data.dt = 10.0/60; %min
            data.subtract_background = 1;
            data.crop_image = 1;
            data.median_filter = 1;
            data.magnification = 100;

        case 'egf_pp1_lyn1',
            data.path = strcat(root, 'fret/EGF_pp1_lyn/');
            data.output_path = strcat(data.path, 'output/');
            data.first_cfp_file = 'PP1_1_CROP_SUB1.009';
            data.index_pattern = {'009', '%03d'};
            data.channel_pattern= {'SUB1', 'SUB2'}; %cfp, fret
            data.image_index = 9:21;
            data.starting_step = 0;
            data.magification = 40;
            data.dt = 30; %sec
            data.diff_coef = 0.1126*2.34*2.34;
            data.subtract_background = 1;
            data.crop_image = 0;
            data.median_filter = 1;
            data.ratio_bound = [0.3, 0.52]; % for lyn
            data.intensity_bound = [400,1600];

        case 'square_5_circle',
            data.path = strcat(root,'simulation/square_5_circle/');
            data.diff_const = 1.0; % mu m^2/s
            data.first_file = 'square_5_circle.TIFF';
            data.boundary_file = 'boundary.mat';
            data = simulation_get_boundary(cell_name, data);
            data.num_bits = 16;
            data.mag = 100; %magnification
            data.rectangle = [1 1 512 512];
            % dt is not needed for computer simulation,
            % but it is needed for compute recovery curve and
            % estimate the diffusion coefficient from the simulation
            data.dt = 0.0156;

        case 'layered_diffusion',
            data.path = strcat(root,'simulation/layered_diffusion/');
            data.diff_const = [5.0 20.0 30.0 1.0]; % mu m^2/s from outside to inside
            data.first_file = 'cell_after_photobleach.PNG';
            image_before = imread(strcat(data.path, 'cell_before_photobleach.PNG'));
            rectangle_file = strcat(data.path, 'image_rectangle.data');
            data.rectangle = get_rectangle(image_before, rectangle_file, 'format', '-ascii');
            temp = imread(strcat(data.path,'cell_after_photobleach.PNG'));
            data.image_0 = imcrop(temp, data.rectangle); clear temp;
            data.num_bits = 8;
            data.mag = 100; %magnification
            % dt is not needed for computer simulation,
            % but it is needed for compute recovery curve and
            % estimate the diffusion coefficient from the simulation
            data.dt = 0.25;

        case 'spot_diffusion',
            data.path = strcat(root,'simulation/spot_diffusion/');
            data.diff_const = [29.0 5.0]; % mu m^2/s
            data.first_file = 'cell_after_photobleach.PNG';
            image_before = imread(strcat(data.path, 'cell_before_photobleach.PNG'));
            rectangle_file = strcat(data.path, 'image_rectangle.data');
            data.rectangle = get_rectangle(image_before, rectangle_file, 'format', '-ascii');
            temp = imread(strcat(data.path,'cell_after_photobleach.PNG'));
            data.image_0 = imcrop(temp, data.rectangle); clear temp;
            data.num_bits = 8;
            data.mag = 100; %magnification
            % additional input for estimate_simulation
            data.dt = 0.25;

        case 'tensor_diffusion',
            data.path = strcat(root,'simulation/tensor_diffusion/');
            data.diff_const = [5.0 100.0 5.0];
            % mu m^2/s
            % first number  - white part;
            % second number - along filament;
            % third number  - perpendicular to filament.
            data.first_file = 'cell_after_photobleach.PNG';
            image_before = imread(strcat(data.path, 'cell_before_photobleach.PNG'));
            rectangle_file = strcat(data.path, 'image_rectangle.data');
            data.rectangle = get_rectangle(image_before, rectangle_file, 'format', '-ascii');
            temp = imread(strcat(data.path,'cell_after_photobleach.PNG'));
            data.image_0 = imcrop(temp, data.rectangle); clear temp;
            data.num_bits = 8;
            data.mag = 100; %magnification
            data.dt = 1; %sec

        case 'tensor_cross_2',
            data.path = strcat(root,'simulation/tensor_cross_2/');
            data.diff_const = [5.0 100.0 5.0];
            % mu m^2/s
            % first number  - white part;
            % second number - along filament;
            % third number  - perpendicular to filament.
            data.first_file = 'cell_after_photobleach.PNG';
            image_before = imread(strcat(data.path, 'cell_before_photobleach.PNG'));
            rectangle_file = strcat(data.path, 'image_rectangle.data');
            data.rectangle = get_rectangle(image_before, rectangle_file, 'format', '-ascii');
            temp = imread(strcat(data.path,'cell_after_photobleach.PNG'));
            data.image_0 = imcrop(temp, data.rectangle); clear temp;
            data.num_bits = 8;
            data.mag = 100; %magnification
            data.dt = 1; %sec

        case 'tensor_cross',
            data.path = strcat(root,'simulation/tensor_cross/');
            data.output_path = strcat(data.path, 'output/');
            data.diff_const = [5.0 100.0 5.0];
            % mu m^2/s
            % first number  - white part;
            % second number - along filament;
            % third number  - perpendicular to filament.
            data.first_file = 'cell_after_photobleach.PNG';
            image_before = imread(strcat(data.path, 'cell_before_photobleach.PNG'));
            rectangle_file = strcat(data.path, 'image_rectangle.data');
            data.rectangle = get_rectangle(image_before, rectangle_file, 'format', '-ascii');
            temp = imread(strcat(data.path,'cell_after_photobleach.PNG'));
            data.image_0 = imcrop(temp, data.rectangle); clear temp;
            data.num_bits = 8;
            data.mag = 100; %magnification
            data.dt = 1; %sec

        case 'photobleach_cell',
            data.path = strcat(root, 'simulation/photobleach_cell/');
            data.diff_const = 29.0;
            data.first_file = 'cell_after_photobleach.PNG';
            image_before = imread(strcat(data.path, 'cell_before_photobleach.PNG'));
            rectangle_file = strcat(data.path, 'image_rectangle.data');
            data.rectangle = get_rectangle(image_before, rectangle_file, 'format', '-ascii');
            temp = imread(strcat(data.path,'cell_after_photobleach.PNG'));
            data.image_0 = imcrop(temp, data.rectangle); clear temp;
            data.num_bits = 8;
            data.mag = 100;
            data.dt = 4.8828e-004;

        case 'photobleach_cell_2',
            % photobleached cell with variable diffusion coefficient
            % diffusion coefficient = 29 when y_image >= 200 pixel or y >= 31.25 um
            %                          1 when y_image <  200 pixel or y <  31.25 um
            data.path = strcat(root, 'simulation/photobleach_cell_2/');
            data.diff_const = [29.0 1.0];
            % criteria for diff_cont(1)
            data.diff_vector_statement = '(tri_centroid(2,:)''>=200)';
            data.first_file = 'cell_after_photobleach.PNG';
            image_before = imread(strcat(data.path, 'cell_before_photobleach.PNG'));
            rectangle_file = strcat(data.path, 'image_rectangle.data');
            data.rectangle = get_rectangle(image_before, rectangle_file, 'format', '-ascii');
            temp = imread(strcat(data.path,'cell_after_photobleach.PNG'));
            data.image_0 = imcrop(temp, data.rectangle); clear temp;
            data.num_bits = 8;
            data.mag = 100;
            data.dt = 9.7656e-004;

        case 'photobleach_cell_2b',
            % photobleached cell with variable diffusion coefficient
            % diffusion coefficient = 29 when y_image >= 200 pixel or y >= 31.25 um
            %                          1 when y_image <  200 pixel or y <  31.25 um
            data.path = strcat(root, 'simulation/photobleach_cell_2b/');
            data.diff_const = [29.0 15.0];
            % criteria for diff_cont(1)
            data.diff_vector_statement = '(tri_centroid(2,:)''>=200)';
            data.first_file = 'cell_after_photobleach.PNG';
            image_before = imread(strcat(data.path, 'cell_before_photobleach.PNG'));
            rectangle_file = strcat(data.path, 'image_rectangle.data');
            data.rectangle = get_rectangle(image_before, rectangle_file, 'format', '-ascii');
            temp = imread(strcat(data.path,'cell_after_photobleach.PNG'));
            data.image_0 = imcrop(temp, data.rectangle); clear temp;
            data.num_bits = 8;
            data.mag = 100;
            data.dt = 4.8828e-004;

        case 'photobleach_cell_3',
            % photobleached cell with variable diffusion coefficient
            % diffusion coefficient = 29 when image distance to (200,200) <= 40 pixel;
            %                          1 when image distance to (200,200) >  40 pixles
            % or D = 29 when distance((x,y)-(31.25, 31.25))<=6.25 um
            %    D =  1 when distance > 6.25 um
            data.path = strcat(root, 'simulation/photobleach_cell_3/');
            data.diff_const = [29.0 1.0];
            % criterial for diff_cont(1)
            data.diff_vector_statement = ...
                '((tri_centroid(1,:)''-200).^2+(tri_centroid(2,:)''-200).^2<=1600)';
            data.first_file = 'cell_after_photobleach.PNG';
            image_before = imread(strcat(data.path, 'cell_before_photobleach.PNG'));
            rectangle_file = strcat(data.path, 'image_rectangle.data');
            data.rectangle = get_rectangle(image_before, rectangle_file, 'format', '-ascii');
            temp = imread(strcat(data.path,'cell_after_photobleach.PNG'));
            data.image_0 = imcrop(temp, data.rectangle); clear temp;
            data.num_bits = 8;
            data.mag = 100;
            data.dt = 9.7656e-004;

        case 'photobleach_cell_3b',
            % photobleached cell with variable diffusion coefficient
            % diffusion coefficient = 29 when image distance to (200,200) <= 40 pixel;
            %                          1 when image distance to (200,200) >  40 pixles
            % or D = 29 when distance((x,y)-(31.25, 31.25))<=6.25 um
            %    D =  1 when distance > 6.25 um
            data.path = strcat(root, 'simulation/photobleach_cell_3b/');
            data.diff_const = [29.0 15.0];

            % criterial for diff_cont(1)
            data.diff_vector_statement = ...
                '((tri_centroid(1,:)''-200).^2+(tri_centroid(2,:)''-200).^2<=1600)';
            data.first_file = 'cell_after_photobleach.PNG';
            image_before = imread(strcat(data.path, 'cell_before_photobleach.PNG'));
            rectangle_file = strcat(data.path, 'image_rectangle.data');
            data.rectangle = get_rectangle(image_before, rectangle_file, 'format', '-ascii');
            temp = imread(strcat(data.path,'cell_after_photobleach.PNG'));
            data.image_0 = imcrop(temp, data.rectangle); clear temp;
            data.num_bits = 8;
            data.mag = 100;
            data.dt = 4.8828e-004; % required for estimate_simulation()
    end;


    if ~isfield(data, 'path'),
        data = init_data_0703_2014(cell_name);
    end;

    % convert from file name to file prefix here
    if ~isfield(data, 'prefix') && isfield(data,'first_file'),
        dot_i = regexpi(data.first_file, '\.');
        data.prefix = data.first_file(1:dot_i-2);
    end;

return;
