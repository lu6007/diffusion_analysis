% function data = diffusion_init_data(cell_name)
% Initialize cell data

% Copyright: Shaoying Lu and Yingxiao Wang 2014
function data = diffusion_init_data(cell_name)
% root='/Users/shirleywu/Desktop/data/';
root = 'D:/sof/data/';
data.cell_name = cell_name;
switch cell_name
    case 'mem17' % the best Lyn-Src cell 
        data.path = strcat(root, 'diffusion_sample/frap/mem17/');
        %data.path = strcat('/Users/shirleywu/Desktop/Research UCSD/Programs/diffusion_1.1/data/frap/mem17/');
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
    case 'egf_pp1_lyn1'
        data.path = strcat(root, 'diffusion_sample/fret/EGF_pp1_lyn/');
        data.output_path = strcat(data.path, 'output/');
        data.first_cfp_file = 'PP1_1_CROP_SUB1.009';
        data.index_pattern = {'009', '%03d'};
        data.channel_pattern= {'SUB1', 'SUB2'}; %cfp, fret
        data.image_index = 9:21;
        data.starting_step = 0;
        data.magification = 40;
        data.dt = 30; %sec
        % um^2/sec, updated 07/05/2014
        data.diff_coef = 0.1126*2.34*2.34;   
        data.subtract_background = 1;
        data.crop_image = 0;
        data.median_filter = 1;
        data.ratio_bound = [0.3, 0.52]; % for lyn
        data.intensity_bound = [400,1600];
    case '03_14_trpc6_gel2'
        data.path = strcat(root, '03_14_FRAP_TRPC6/gel-2/');
        data.output_path = strcat(data.path, 'output/');
        data.first_file = 'G-21.001';
        data.index_pattern = {'001', '%03d'};
        data.image_index = 8:(8+11);
        data.index_before = 1:5;
        data.index_after = 6:12;
        % photobleach between frames 12 and 13
        data.photobleach_time = 10506;
        data.threshold = 0.015;
        data.subtract_background = 1;
        data.crop_image = 1;
        data.median_filter = 1;
        data.magnification = 100;        
    case 'square_5_circle'
        %data.path = strcat(root,'03_09_2012_for_ru/2/');
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
    case 'layered_diffusion'
        data.path = strcat(root,'diffusion_sample/simulation/layered_diffusion/');
        data.diff_const = [5.0 20.0 30.0 1.0]; % mu m^2/s from outside to inside
        %data.diff_const = [5.0 20.0 30.0 10.0 1.0]; % changed to 5 layers 03/19/15
        data.first_file = 'cell_after_photobleach.PNG';
        %
        image_before = imread(strcat(data.path, 'cell_before_photobleach.PNG'));
        rectangle_file = strcat(data.path, 'image_rectangle.data');
        data.rectangle = get_rectangle(image_before, rectangle_file, 'format', '-ascii');
        temp = imread(strcat(data.path,'cell_after_photobleach.PNG'));
        data.image_0 = imcrop(temp, data.rectangle); clear temp;
        %
%       data = simulation_get_boundary(cell_name, data);
        data.num_bits = 8;
        data.mag = 100; %magnification
        % dt is not needed for computer simulation,
        % but it is needed for compute recovery curve and
        % estimate the diffusion coefficient from the simulation
        data.dt = 0.002;
    case 'spot_diffusion'
        % computer_simulation input
        data.path = strcat(root,'diffusion_sample/simulation/spot_diffusion/');
        data.diff_const = [29.0 5.0]; % mu m^2/s
        data.first_file = 'cell_after_photobleach.PNG';
        %
        image_before = imread(strcat(data.path, 'cell_before_photobleach.PNG'));
        rectangle_file = strcat(data.path, 'image_rectangle.data');
        data.rectangle = get_rectangle(image_before, rectangle_file, 'format', '-ascii');
        temp = imread(strcat(data.path,'cell_after_photobleach.PNG'));
        data.image_0 = imcrop(temp, data.rectangle); clear temp;
        %
%       data = simulation_get_boundary(cell_name, data);
        data.num_bits = 8;
        data.mag = 100; %magnification
        % additional input for estimate_simulation
        data.dt = 9.7656e-04;
    case 'tensor_diffusion'
        % computer_simulation input
        data.path = strcat(root,'diffusion_sample/simulation/tensor_diffusion/');
        data.diff_const = [5.0 100.0 5.0]; % mu m^2/s %% first number - 
        %white part; second number - along filament; third number - perpendicular 
        %to filament.
        %data.diff_const ={ [5.0 0; 0 5.0], [60.0 0; 0 0.0]}; % mu m^2/s
        data.first_file = 'cell_after_photobleach.PNG';
        %
        image_before = imread(strcat(data.path, 'cell_before_photobleach.PNG'));
        rectangle_file = strcat(data.path, 'image_rectangle.data');
        data.rectangle = get_rectangle(image_before, rectangle_file, 'format', '-ascii');
        temp = imread(strcat(data.path,'cell_after_photobleach.PNG'));
        data.image_0 = imcrop(temp, data.rectangle); clear temp;
        %
%       data = simulation_get_boundary(cell_name, data);
        data.num_bits = 8;
        data.mag = 100; %magnification
        % additional input for estimate_simulation
        data.dt = 9.7656e-04; %sec
        
    case 'tensor_cross'
        % computer_simulation input
        data.path = strcat(root,'diffusion_sample/simulation/tensor_cross/');
        data.output_path = strcat(data.path, 'output/');
        data.diff_const = [5.0 100.0 5.0]; % mu m^2/s %% first number - 
        %white part; second number - along filament; third number - perpendicular 
        %to filament.
        %data.diff_const ={ [5.0 0; 0 5.0], [60.0 0; 0 0.0]}; % mu m^2/s
        data.first_file = 'cell_after_photobleach.PNG';
        %
        image_before = imread(strcat(data.path, 'cell_before_photobleach.PNG'));
        rectangle_file = strcat(data.path, 'image_rectangle.data');
        data.rectangle = get_rectangle(image_before, rectangle_file, 'format', '-ascii');
        temp = imread(strcat(data.path,'cell_after_photobleach.PNG'));
        data.image_0 = imcrop(temp, data.rectangle); clear temp;
        %
%       data = simulation_get_boundary(cell_name, data);
        data.num_bits = 8;
        data.mag = 100; %magnification
        % additional input for estimate_simulation
        data.dt = 9.7656e-04; %sec

    case 'photobleach_cell'
        %data.path = 'C:/sof/diffusion_1.1/data/simulation/';
        data.path = strcat(root, 'diffusion_sample/simulation/photobleach_cell/');
        data.diff_const = 29.0;
        data.first_file = 'cell_after_photobleach.PNG';
        %
        image_before = imread(strcat(data.path, 'cell_before_photobleach.PNG'));
        rectangle_file = strcat(data.path, 'image_rectangle.data');
        data.rectangle = get_rectangle(image_before, rectangle_file, 'format', '-ascii');
        temp = imread(strcat(data.path,'cell_after_photobleach.PNG'));
        data.image_0 = imcrop(temp, data.rectangle); clear temp;
        %
%       data = simulation_get_boundary(cell_name, data);
        data.num_bits = 8;
        data.mag = 100;
        %
        data.dt = 4.8828e-004; 
    case 'photobleach_cell_2' 
        % photobleached cell with variable diffusion coefficient
        % diffusion coefficient = 29 when y_image>=200 pixel or y>=31.25 um
        %                         1 when y_image<200 pixel or y<31.25 um
        data.path = strcat(root, 'diffusion_sample/simulation/photobleach_cell_2/');
        data.diff_const = [29.0 1.0];
        % criteria for diff_cont(1)
        data.diff_vector_statement = '(tri_centroid(2,:)''>=200)';
        data.first_file = 'cell_after_photobleach.PNG';
        %
        image_before = imread(strcat(data.path, 'cell_before_photobleach.PNG'));
        rectangle_file = strcat(data.path, 'image_rectangle.data');
        data.rectangle = get_rectangle(image_before, rectangle_file, 'format', '-ascii');
        temp = imread(strcat(data.path,'cell_after_photobleach.PNG'));
        data.image_0 = imcrop(temp, data.rectangle); clear temp;
        %
%       data = simulation_get_boundary(cell_name, data);
        data.num_bits = 8;
        data.mag = 100;
        %
        data.dt = 9.7656e-004;
    case 'photobleach_cell_2b' 
        % photobleached cell with variable diffusion coefficient
        % diffusion coefficient = 29 when y_image>=200 pixel or y>=31.25 um
        %                         1 when y_image<200 pixel or y<31.25 um
        data.path = strcat(root, 'diffusion_sample/simulation/photobleach_cell_2b/');
        data.diff_const = [29.0 15.0];
        % criteria for diff_cont(1)
        data.diff_vector_statement = '(tri_centroid(2,:)''>=200)';
        data.first_file = 'cell_after_photobleach.PNG';
        %
        image_before = imread(strcat(data.path, 'cell_before_photobleach.PNG'));
        rectangle_file = strcat(data.path, 'image_rectangle.data');
        data.rectangle = get_rectangle(image_before, rectangle_file, 'format', '-ascii');
        temp = imread(strcat(data.path,'cell_after_photobleach.PNG'));
        data.image_0 = imcrop(temp, data.rectangle); clear temp;
        %
%       data = simulation_get_boundary(cell_name, data);
        data.num_bits = 8;
        data.mag = 100;
        %
        data.dt = 4.8828e-004;
    case 'photobleach_cell_3'
        % photobleached cell with variable diffusion coefficient
        % diffusion coefficient = 29 when image distance to (200,200)<=40 pixel;
        % 1 when image distance to (200,200)> 40 pixles
        % or D = 29 when distance((x,y)-(31.25, 31.25))<=6.25 um
        %    D = 1 when distance > 6.25 um
        data.path = strcat(root, 'diffusion_sample/simulation/photobleach_cell_3/');
        data.diff_const = [29.0 1.0];
        % criterial for diff_cont(1)
        data.diff_vector_statement = ...
            '((tri_centroid(1,:)''-200).^2+(tri_centroid(2,:)''-200).^2<=1600)';
        data.first_file = 'cell_after_photobleach.PNG';
        %
        image_before = imread(strcat(data.path, 'cell_before_photobleach.PNG'));
        rectangle_file = strcat(data.path, 'image_rectangle.data');
        data.rectangle = get_rectangle(image_before, rectangle_file, 'format', '-ascii');
        temp = imread(strcat(data.path,'cell_after_photobleach.PNG'));
        data.image_0 = imcrop(temp, data.rectangle); clear temp;
        %
%       data = simulation_get_boundary(cell_name, data);
        data.num_bits = 8;
        data.mag = 100;
        %
        data.dt = 9.7656e-004;
    case 'photobleach_cell_3b'
        % photobleached cell with variable diffusion coefficient
        % diffusion coefficient = 29 when image distance to (200,200)<=40 pixel;
        % 1 when image distance to (200,200)> 40 pixles
        % or D = 29 when distance((x,y)-(31.25, 31.25))<=6.25 um
        %    D = 15 when distance > 6.25 um
        data.path = strcat(root, 'diffusion_sample/simulation/photobleach_cell_3b/');
        data.diff_const = [29.0 15.0];
        % criterial for diff_cont(1)
        data.diff_vector_statement = ...
            '((tri_centroid(1,:)''-200).^2+(tri_centroid(2,:)''-200).^2<=1600)';
        data.first_file = 'cell_after_photobleach.PNG';
        %
        image_before = imread(strcat(data.path, 'cell_before_photobleach.PNG'));
        rectangle_file = strcat(data.path, 'image_rectangle.data');
        data.rectangle = get_rectangle(image_before, rectangle_file, 'format', '-ascii');
        temp = imread(strcat(data.path,'cell_after_photobleach.PNG'));
        data.image_0 = imcrop(temp, data.rectangle); clear temp;
        %
%       data = simulation_get_boundary(cell_name, data);
        data.num_bits = 8;
        data.mag = 100;
        %
        data.dt = 4.8828e-004; % required for estimate_simulation()
        
    % 0702/2014 More data from Qin Peng
    case 'h3k9_wt2' 
        data.path = strcat(root, 'qin_frap/h3k9_0610/WT/WT MEF 2 15s/');
        data.output_path = strcat(data.path, 'output/');
        data.first_file = 'H3K924.009';
        data.index_pattern = {'009', '%03d'};
        data.image_index = 9:17;
        data.index_before = 1:3;
        data.index_after = 4:length(data.image_index);
        % photobleach between frames 11 and 12
        data.photobleach_time = 50020;
        data.threshold = 0;
        data.subtract_background = 1;
        data.crop_image = 1;
        data.median_filter = 1;
        data.magnification = 100;
        data.scale_image = 16;
    case 'h3k9_wt3' %not good
        data.path = strcat(root, 'qin_frap/h3k9_0610/WT/WT MEF 3 20s/');
        data.output_path = strcat(data.path, 'output/');
        data.first_file = 'H3K924.007';
        data.index_pattern = {'007', '%03d'};
        data.image_index = 7:15;
        data.index_before = 1:3;
        data.index_after = 4:length(data.image_index);
        % photobleach between frames 9 and 10
        data.photobleach_time = 42907;
        data.threshold = 0;
        data.subtract_background = 1;
        data.crop_image = 1;
        data.median_filter = 1;
        data.magnification = 100;
        data.scale_image = 16;
    case 'h3k9_wt5'
        data.path = strcat(root, 'qin_frap/h3k9_0610/WT/WT MEF 5 25s/');
        data.output_path = strcat(data.path, 'output/');
        data.first_file = 'H3K924.007';
        data.index_pattern = {'007', '%03d'};
        data.image_index = 7:15;
        data.index_before = 1:3;
        data.index_after = 4:length(data.image_index);
        % photobleach between frames 9 and 10
        data.photobleach_time = 39650;
        data.threshold = 0;
        data.subtract_background = 1;
        data.crop_image = 1;
        data.median_filter = 1;
        data.magnification = 100;
        data.brightness_factor = 1.0;
    case 'h3k9_wt6'
        data.path = strcat(root, 'qin_frap/h3k9_0610/WT/WT MEF 6 25s/');
        data.output_path = strcat(data.path, 'output/');
        data.first_file = 'H3K924.010';
        data.index_pattern = {'010', '%03d'};
        data.image_index = 10:18;
        data.index_before = 1:3;
        data.index_after = 4:length(data.image_index);
        data.threshold = 0;
        % photobleach between frames 12 and 13
        data.photobleach_time = 61517;
        data.subtract_background = 1;
        data.crop_image = 1;
        data.median_filter = 1;
        data.magnification = 100;
        data.scale_image = 16;
    case 'h3k9_wt7'
        data.path = strcat(root, 'qin_frap/h3k9_0610/WT/WT MEF 7 25s/');
        data.output_path = strcat(data.path, 'output/');
        data.first_file = 'H3K924.007';
        data.index_pattern = {'007', '%03d'};
        data.image_index = 7:15;
        data.index_before = 1:3;
        data.index_after = 4:length(data.image_index);
        data.threshold = 0;
        % photobleach between frames 9 and 10
        data.photobleach_time = 29306;
        data.subtract_background = 1;
        data.crop_image = 1;
        data.median_filter = 1;
        data.magnification = 100;
        data.scale_image = 16;
    case 'h3k9_wt8'
        data.path = strcat(root, 'qin_frap/h3k9_0610/WT/WT MEF 8 25s/');
        data.output_path = strcat(data.path, 'output/');
        data.first_file = 'H3K924.012';
        data.index_pattern = {'012', '%03d'};
        data.image_index = 12:20;
        data.index_before = 1:3;
        data.index_after = 4:length(data.image_index);
        data.threshold = 0;
        % photobleach between frames 9 and 10
        data.photobleach_time = 35378;
        data.subtract_background = 1;
        data.crop_image = 1;
        data.median_filter = 1;
        data.magnification = 100;
        data.scale_image = 16;
    case 'h3k9_wt9'
        data.path = strcat(root, 'qin_frap/h3k9_0610/WT/WT MEF 9 25s/');
        data.output_path = strcat(data.path, 'output/');
        data.first_file = 'H3K924.009';
        data.index_pattern = {'009', '%03d'};
        data.image_index = 9:17;
        data.index_before = 1:3;
        data.index_after = 4:length(data.image_index);
        data.threshold = 0;
        % photobleach between frames 9 and 10
        data.photobleach_time = 27363;
        data.subtract_background = 1;
        data.crop_image = 1;
        data.median_filter = 1;
        data.magnification = 100;
        data.scale_image = 16;
    case 'h3k9_wt10'
        data.path = strcat(root, 'qin_frap/h3k9_0610/WT/WT MEF 10 25s/');
        data.output_path = strcat(data.path, 'output/');
        data.first_file = 'H3K924.009';
        data.index_pattern = {'009', '%03d'};
        data.image_index = 9:17;
        data.index_before = 1:3;
        data.index_after = 4:length(data.image_index);
        data.threshold = 0;
        % photobleach between frames 11 and 12
        data.photobleach_time = 0;
        data.subtract_background = 1;
        data.crop_image = 0;
        data.median_filter = 1;
        data.magnification = 100;
        data.scale_image = 16;
    case 'h3k9_wt11'
        data.path = strcat(root, 'qin_frap/h3k9_0610/WT/WT MEF 11 25s/');
        data.output_path = strcat(data.path, 'output/');
        data.first_file = 'H3K924.008';
        data.index_pattern = {'008', '%03d'};
        data.image_index = 8:17;
        data.index_before = 1:3;
        data.index_after = 4:length(data.image_index);
        data.threshold = 0;
        % photobleach between frames 10 and 11
        data.photobleach_time = 0;
        data.subtract_background = 1;
        data.crop_image = 0;
        data.median_filter = 1;
        data.magnification = 100;
        data.scale_image = 16;
    case 'h3k9_wt12'
        data.path = strcat(root, 'qin_frap/h3k9_0610/WT/WT MEF 12 25s/');
        data.output_path = strcat(data.path, 'output/');
        data.first_file = 'H3K924.008';
        data.index_pattern = {'008', '%03d'};
        data.image_index = 8:17;
        data.index_before = 1:3;
        data.index_after = 4:length(data.image_index);
        data.threshold = 0;
        % photobleach between frames 10 and 11
        data.photobleach_time = 0;
        data.subtract_background = 1;
        data.crop_image = 0;
        data.median_filter = 1;
        data.magnification = 100;
        data.scale_image = 16;
     case 'h3k9_wt13'
        data.path = strcat(root, 'qin_frap/h3k9_0610/WT/WT MEF 13 25s/');
        data.output_path = strcat(data.path, 'output/');
        data.first_file = 'H3K924.008';
        data.index_pattern = {'008', '%03d'};
        data.image_index = 8:17;
        data.index_before = 1:3;
        data.index_after = 4:length(data.image_index);
        data.threshold = 0;
        % photobleach between frames 10 and 11
        data.photobleach_time = 0;
        data.subtract_background = 1;
        data.crop_image = 0;
        data.median_filter = 1;
        data.magnification = 100;
        data.scale_image = 16;
     case 'h3k9_s10a1'
        data.path = strcat(root, 'qin_frap/h3k9_0610/S10A/1 25s/');
        data.output_path = strcat(data.path, 'output/');
        data.first_file = 'H3K924.008';
        data.index_pattern = {'008', '%03d'};
        data.image_index = 9:18;
        data.index_before = 1:3;
        data.index_after = 4:length(data.image_index);
        data.threshold = 0;
        % photobleach between frames 11 and 12
        data.photobleach_time = 0;
        data.subtract_background = 1;
        data.crop_image = 0;
        data.median_filter = 1;
        data.magnification = 100;
        data.scale_image = 16;
     case 'h3k9_s10a2'
        data.path = strcat(root, 'qin_frap/h3k9_0610/S10A/2 25s/');
        data.output_path = strcat(data.path, 'output/');
        data.first_file = 'H3K924.006';
        data.index_pattern = {'006', '%03d'};
        data.image_index = 6:15;
        data.index_before = 1:3;
        data.index_after = 4:length(data.image_index);
        data.threshold = 0;
        % photobleach between frames 8 and 9
        data.photobleach_time = 0;
        data.subtract_background = 1;
        data.crop_image = 0;
        data.median_filter = 1;
        data.magnification = 100;
        data.scale_image = 16;
     case 'h3k9_s10a3'
        data.path = strcat(root, 'qin_frap/h3k9_0610/S10A/3 25s/');
        data.output_path = strcat(data.path, 'output/');
        data.first_file = 'H3K924.007';
        data.index_pattern = {'007', '%03d'};
        data.image_index = 7:16;
        data.index_before = 1:3;
        data.index_after = 4:length(data.image_index);
        data.threshold = 0;
        % photobleach between frames 9 and 10
        data.photobleach_time = 0;
        data.subtract_background = 1;
        data.crop_image = 0;
        data.median_filter = 1;
        data.magnification = 100;
        data.scale_image = 16;
     case 'h3k9_s10a4'
        data.path = strcat(root, 'qin_frap/h3k9_0610/S10A/4 25s/');
        data.output_path = strcat(data.path, 'output/');
        data.first_file = 'H3K924.007';
        data.index_pattern = {'007', '%03d'};
        data.image_index = 7:16;
        data.index_before = 1:3;
        data.index_after = 4:length(data.image_index);
        data.threshold = 0;
        % photobleach between frames 9 and 10
        data.photobleach_time = 0;
        data.subtract_background = 1;
        data.crop_image = 0;
        data.median_filter = 1;
        data.magnification = 100;
        data.scale_image = 16;
     case 'h3k9_s10a5'
        data.path = strcat(root, 'qin_frap/h3k9_0610/S10A/5 25s/');
        data.output_path = strcat(data.path, 'output/');
        data.first_file = 'H3K924.009';
        data.index_pattern = {'009', '%03d'};
        data.image_index = 9:18;
        data.index_before = 1:3;
        data.index_after = 4:length(data.image_index);
        data.threshold = 0;
        % photobleach between frames 11 and 12
        data.photobleach_time = 0;
        data.subtract_background = 1;
        data.crop_image = 0;
        data.median_filter = 1;
        data.magnification = 100;
        data.scale_image = 16;
     case 'h3k9_s10a6'
        data.path = strcat(root, 'qin_frap/h3k9_0610/S10A/6 25s/');
        data.output_path = strcat(data.path, 'output/');
        data.first_file = 'H3K924.007';
        data.index_pattern = {'007', '%03d'};
        data.image_index = 7:16;
        data.index_before = 1:3;
        data.index_after = 4:length(data.image_index);
        data.threshold = 0;
        % photobleach between frames 9 and 10
        data.photobleach_time = 0;
        data.subtract_background = 1;
        data.crop_image = 1;
        data.median_filter = 1;
        data.magnification = 100;
        data.scale_image = 16;
     case 'h3k9_s10a7'
        data.path = strcat(root, 'qin_frap/h3k9_0610/S10A/7 25s/');
        data.output_path = strcat(data.path, 'output/');
        data.first_file = 'H3K924.009';
        data.index_pattern = {'009', '%03d'};
        data.image_index = 9:18;
        data.index_before = 1:3;
        data.index_after = 4:length(data.image_index);
        data.threshold = 0;
        % photobleach between frames 11 and 12
        data.photobleach_time = 0;
        data.subtract_background = 1;
        data.crop_image = 1;
        data.median_filter = 1;
        data.magnification = 100;
        data.scale_image = 16;
     case 'h3k9_s10a8'
        data.path = strcat(root, 'qin_frap/h3k9_0610/S10A/8 25s/');
        data.output_path = strcat(data.path, 'output/');
        data.first_file = 'H3K924.010';
        data.index_pattern = {'010', '%03d'};
        data.image_index = 10:19;
        data.index_before = 1:3;
        data.index_after = 4:length(data.image_index);
        data.threshold = 0;
        % photobleach between frames 12 and 13
        data.photobleach_time = 0;
        data.subtract_background = 1;
        data.crop_image = 1;
        data.median_filter = 1;
        data.magnification = 100;
        data.scale_image = 1;
     case 'h3k9_s10a9' % cell looks weird, not good. 
        data.path = strcat(root, 'qin_frap/h3k9_0610/S10A/9 25s/');
        data.output_path = strcat(data.path, 'output/');
        data.first_file = 'H3K924.009';
        data.index_pattern = {'009', '%03d'};
        data.image_index = 9:18;
        data.index_before = 1:3;
        data.index_after = 4:length(data.image_index);
        data.threshold = 0;
        % photobleach between frames 11 and 12
        data.photobleach_time = 0;
        data.subtract_background = 1;
        data.crop_image = 1;
        data.median_filter = 1;
        data.magnification = 100;
        data.scale_image = 16;
     case 'h3k9_at1' 
        data.path = strcat(root, 'qin_frap/h3k9_0610/AT/1 25s/');
        data.output_path = strcat(data.path, 'output/');
        data.first_file = 'H3K924.007';
        data.index_pattern = {'007', '%03d'};
        data.image_index = 7:16;
        data.index_before = 1:3;
        data.index_after = 4:length(data.image_index);
        data.threshold = 0;
        % photobleach between frames 9 and 10
        data.photobleach_time = 0;
        data.subtract_background = 1;
        data.crop_image = 1;
        data.median_filter = 1;
        data.magnification = 100;
        data.scale_image = 16;
     case 'h3k9_at2' 
        data.path = strcat(root, 'qin_frap/h3k9_0610/AT/2 25s/');
        data.output_path = strcat(data.path, 'output/');
        data.first_file = 'H3K924.007';
        data.index_pattern = {'007', '%03d'};
        data.image_index = 7:16;
        data.index_before = 1:3;
        data.index_after = 4:length(data.image_index);
        data.threshold = 0;
        % photobleach between frames 9 and 10
        data.photobleach_time = 0;
        data.subtract_background = 1;
        data.crop_image = 0;
        data.median_filter = 1;
        data.magnification = 100;
        data.scale_image = 1;
     case 'h3k9_at3' 
        data.path = strcat(root, 'qin_frap/h3k9_0610/AT/3 25s/');
        data.output_path = strcat(data.path, 'output/');
        data.first_file = 'H3K924.007';
        data.index_pattern = {'007', '%03d'};
        data.image_index = 7:16;
        data.index_before = 1:3;
        data.index_after = 4:length(data.image_index);
        data.threshold = 0;
        % photobleach between frames 9 and 10
        data.photobleach_time = 0;
        data.subtract_background = 1;
        data.crop_image = 0;
        data.median_filter = 1;
        data.magnification = 100;
        data.scale_image = 16;
     case 'h3k9_at4' 
        data.path = strcat(root, 'qin_frap/h3k9_0610/AT/4 25s/');
        data.output_path = strcat(data.path, 'output/');
        data.first_file = 'H3K924.007';
        data.index_pattern = {'007', '%03d'};
        data.image_index = 7:16;
        data.index_before = 1:3;
        data.index_after = 4:length(data.image_index);
        data.threshold = 0;
        % photobleach between frames 9 and 10
        data.photobleach_time = 0;
        data.subtract_background = 1;
        data.crop_image = 0;
        data.median_filter = 1;
        data.magnification = 100;
        data.scale_image = 16;
     case 'h3k9_at5' 
        data.path = strcat(root, 'qin_frap/h3k9_0610/AT/5 25s/');
        data.output_path = strcat(data.path, 'output/');
        data.first_file = 'H3K924.008';
        data.index_pattern = {'008', '%03d'};
        data.image_index = 8:17;
        data.index_before = 1:3;
        data.index_after = 4:length(data.image_index);
        data.threshold = 0;
        % photobleach between frames 10 and 11
        data.photobleach_time = 0;
        data.subtract_background = 1;
        data.crop_image = 0;
        data.median_filter = 1;
        data.magnification = 100;
        data.scale_image = 16;
     case 'h3k9_at6' 
        data.path = strcat(root, 'qin_frap/h3k9_0610/AT/6 25s/');
        data.output_path = strcat(data.path, 'output/');
        data.first_file = 'H3K924.008';
        data.index_pattern = {'008', '%03d'};
        data.image_index = 8:17;
        data.index_before = 1:3;
        data.index_after = 4:length(data.image_index);
        data.threshold = 0;
        % photobleach between frames 10 and 11
        data.photobleach_time = 0;
        data.subtract_background = 1;
        data.crop_image = 0;
        data.median_filter = 1;
        data.magnification = 100;
        data.scale_image = 16;
     case 'h3k9_at7' 
        data.path = strcat(root, 'qin_frap/h3k9_0610/AT/7 25s/');
        data.output_path = strcat(data.path, 'output/');
        data.first_file = 'H3K924.007';
        data.index_pattern = {'007', '%03d'};
        data.image_index = 7:16;
        data.index_before = 1:3;
        data.index_after = 4:length(data.image_index);
        data.threshold = 0;
        % photobleach between frames 9 and 10
        data.photobleach_time = 0;
        data.subtract_background = 1;
        data.crop_image = 0;
        data.median_filter = 1;
        data.magnification = 100;
        data.scale_image = 16;
        
        %%% 08/20/2014 Bo's FRAP of membrane Dil
    case '0818_down_sh13'
        data.path = strcat(root, 'bo_frap\0818_frap\shear_start\sh13\');
        data.output_path = strcat(data.path, 'output/');
        data.first_file = 'FRAP1.011';
        data.index_pattern = {'011', '%03d'};
        data.image_index = 11:20;
        % photbleach between frames 13 and 14
        data.index_before = 1:3;
        data.index_after = 4:length(data.image_index);
        data.threshold = 0;
        data.photobleach_time = 0;
        data.subtract_background = 1;
        data.crop_image = 0;
        data.median_filter = 1;
        data.magnification = 100;
        data.scale_image = 1;
    case '0818_down_sh17'
        data.path = strcat(root, 'bo_frap\0818_frap\shear_start\sh17\');
        data.output_path = strcat(data.path, 'output/');
        data.first_file = 'FRAP1.024';
        data.index_pattern = {'024', '%03d'};
        data.image_index = 24:33;
        % photbleach between frames 26 and 27
        data.index_before = 1:3;
        data.index_after = 4:length(data.image_index);
        data.threshold = 0;
        data.photobleach_time = 0;
        data.subtract_background = 1;
        data.crop_image = 0;
        data.median_filter = 1;
        data.magnification = 100;
        data.scale_image = 1;
    case '0818_down_sh18'
        data.path = strcat(root, 'bo_frap\0818_frap\shear_start\sh18\');
        data.output_path = strcat(data.path, 'output/');
        data.first_file = 'FRAP1.012';
        data.index_pattern = {'012', '%03d'};
        data.image_index = 12:21;
        % photbleach between frames 14 and 15
        data.index_before = 1:3;
        data.index_after = 4:length(data.image_index);
        data.threshold = 0;
        data.photobleach_time = 0;
        data.subtract_background = 1;
        data.crop_image = 1;
        data.median_filter = 1;
        data.magnification = 100;
        data.scale_image = 1;
    case '0818_down_sh19'
        data.path = strcat(root, 'bo_frap\0818_frap\shear_start\sh19\');
        data.output_path = strcat(data.path, 'output/');
        data.first_file = 'FRAP1.010';
        data.index_pattern = {'010', '%03d'};
        data.image_index = (10:19);
        % photbleach between frames 12 and 14
        data.index_before = 1:3;
        data.index_after = 4:length(data.image_index);
        data.threshold = 0;
        data.photobleach_time = 0;
        data.subtract_background = 1;
        data.crop_image = 0;
        data.median_filter = 1;
        data.magnification = 100;
        data.scale_image = 1;
        data.brightness_factor = 0.4;
    case '0818_up_sh22'
        data.path = strcat(root, 'bo_frap\0818_frap\shear_upstream\sh22\');
        data.output_path = strcat(data.path, 'output/');
        data.first_file = 'FRAP1.017';
        data.index_pattern = {'017', '%03d'};
        data.image_index = 17:26;
        % photbleach between frames 12 and 14
        data.index_before = 1:3;
        data.index_after = 4:length(data.image_index);
        data.threshold = 0;
        data.photobleach_time = 0;
        data.subtract_background = 1;
        data.crop_image = 0;
        data.median_filter = 1;
        data.magnification = 100;
        data.scale_image = 1;
        data.brightness_factor =0.7;
    case '0818_up_sh24'
        data.path = strcat(root, 'bo_frap\0818_frap\shear_upstream\sh24\');
        data.output_path = strcat(data.path, 'output/');
        data.first_file = 'FRAP1.016';
        data.index_pattern = {'016', '%03d'};
        data.image_index = [16:18, 20:26];
        % photbleach between frames 18 and 20
        data.index_before = 1:3;
        data.index_after = 4:length(data.image_index);
        data.threshold = 0;
        data.photobleach_time = 0;
        data.subtract_background = 1;
        data.crop_image = 0;
        data.median_filter = 1;
        data.magnification = 100;
        data.scale_image = 1;
        data.brightness_factor =1;
    case '0818_up_sh27'
        data.path = strcat(root, 'bo_frap\0818_frap\shear_upstream\sh27\');
        data.output_path = strcat(data.path, 'output/');
        data.first_file = 'FRAP1.016';
        data.index_pattern = {'016', '%03d'};
        data.image_index = (16:25);
        % photbleach between frames 18 and 19
        data.index_before = 1:3;
        data.index_after = 4:length(data.image_index);
        data.threshold = 0;
        data.photobleach_time = 0;
        data.subtract_background = 1;
        data.crop_image = 0;
        data.median_filter = 1;
        data.magnification = 100;
        data.scale_image = 1;
        data.brightness_factor =1;
    case '0818_up_sh29'
        data.path = strcat(root, 'bo_frap\0818_frap\shear_upstream\sh29\');
        data.output_path = strcat(data.path, 'output/');
        data.first_file = 'FRAP1.013';
        data.index_pattern = {'013', '%03d'};
        data.image_index = (13:22);
        % photbleach between frames 15 and 16
        data.index_before = 1:3;
        data.index_after = 4:length(data.image_index);
        data.threshold = 0;
        data.photobleach_time = 0;
        data.subtract_background = 1;
        data.crop_image = 1;
        data.median_filter = 1;
        data.magnification = 100;
        data.scale_image = 1;
        data.brightness_factor =1;
    case '0818_up_sh30'
        data.path = strcat(root, 'bo_frap\0818_frap\shear_upstream\sh30\');
        data.output_path = strcat(data.path, 'output/');
        data.first_file = 'FRAP1.013';
        data.index_pattern = {'013', '%03d'};
        data.image_index = (13:22);
        % photbleach between frames 15 and 16
        data.index_before = 1:3;
        data.index_after = 4:length(data.image_index);
        data.threshold = 0;
        data.photobleach_time = 0;
        data.subtract_background = 1;
        data.crop_image = 0;
        data.median_filter = 1;
        data.magnification = 100;
        data.scale_image = 1;
        data.brightness_factor =1;
    case '0818_control_sh03'
        data.path = strcat(root, 'bo_frap\0818_frap\noshear\sh03\');
        data.output_path = strcat(data.path, 'output/');
        data.first_file = 'FRAP1.008';
        data.index_pattern = {'008', '%03d'};
        data.image_index = (8:17);
        % photbleach between frames 10 and 11
        data.index_before = 1:3;
        data.index_after = 4:length(data.image_index);
        data.threshold = 0;
        data.photobleach_time = 0;
        data.subtract_background = 1;
        data.crop_image = 1;
        data.median_filter = 1;
        data.magnification = 100;
        data.scale_image = 1;
        data.brightness_factor =1;
    case '0813_control21'
        data.path = strcat(root, 'bo_frap\0813_2014\down2\1-24\');
        data.output_path = strcat(data.path, 'output/');
        data.first_file = 'FRAP1.009';
        data.index_pattern = {'009', '%03d'};
        data.image_index = (12:21);
        % photbleach between frames 10 and 11
        data.index_before = 1:3;
        data.index_after = 4:length(data.image_index);
        data.threshold = 0;
        data.photobleach_time = 0;
        data.subtract_background = 1;
        data.crop_image = 0;
        data.median_filter = 1;
        data.magnification = 100;
        data.scale_image = 1;
        data.brightness_factor =1;
    case '0813_control22'
        data.path = strcat(root, 'bo_frap\0813_2014\down2\50-71\');
        data.output_path = strcat(data.path, 'output/');
        data.first_file = 'FRAP1.061';
        data.index_pattern = {'061', '%03d'};
        data.image_index = [61:63, 65:71];
        % photbleach between frames 10 and 11
        data.index_before = 1:3;
        data.index_after = 4:length(data.image_index);
        data.threshold = 0;
        data.photobleach_time = 0;
        data.subtract_background = 1;
        data.crop_image = 0;
        data.median_filter = 1;
        data.magnification = 100;
        data.scale_image = 1;
        data.brightness_factor =1;
    case '0826_up01'
        data.path = strcat(root, 'bo_frap\0826_frap\shear_upstream\sh01\');
        data.output_path = strcat(data.path, 'output/');
        data.first_file = 'FRAP1.057';
        data.index_pattern = {'057', '%03d'};
        data.image_index = [57:59, 60:66];
        % photbleach between frames 10 and 11
        data.index_before = 1:3;
        data.index_after = 4:length(data.image_index);
        data.threshold = 0;
        data.photobleach_time = 0;
        data.subtract_background = 1;
        data.crop_image = 1;
        data.median_filter = 1;
        data.magnification = 100;
        data.scale_image = 1;
        data.brightness_factor =0.4;
    case '0826_up02' % not good
        data.path = strcat(root, 'bo_frap\0826_frap\shear_upstream\sh02\');
        data.output_path = strcat(data.path, 'output/');
        data.first_file = 'FRAP1.013';
        data.index_pattern = {'013', '%03d'};
        data.image_index = [13:15, 16:22];
        % photbleach between frames 15 and 16
        data.index_before = 1:3;
        data.index_after = 4:length(data.image_index);
        data.threshold = 0;
        data.photobleach_time = 0;
        data.subtract_background = 1;
        data.crop_image = 1;
        data.median_filter = 1;
        data.magnification = 100;
        data.scale_image = 1;
        data.brightness_factor =0.5;
    case '0826_up05'
        data.path = strcat(root, 'bo_frap\0826_frap\shear_upstream\sh05\');
        data.output_path = strcat(data.path, 'output/');
        data.first_file = 'FRAP1.018';
        data.index_pattern = {'018', '%03d'};
        data.image_index = [16:18, 19:25];
        % photbleach between frames 18 and 19
        data.index_before = 1:3;
        data.index_after = 4:length(data.image_index);
        data.threshold = 0;
        data.photobleach_time = 0;
        data.subtract_background = 1;
        data.crop_image = 1;
        data.median_filter = 1;
        data.magnification = 100;
        data.scale_image = 1;
        data.brightness_factor =0.5;
    case '0826_up08'
        data.path = strcat(root, 'bo_frap\0826_frap\shear_upstream\sh08\');
        data.output_path = strcat(data.path, 'output/');
        data.first_file = 'FRAP1.012';
        data.index_pattern = {'012', '%03d'};
        data.image_index = [13:15, 16:22];
        % photbleach between frames 15 and 16
        data.index_before = 1:3;
        data.index_after = 4:length(data.image_index);
        data.threshold = 0;
        data.photobleach_time = 0;
        data.subtract_background = 1;
        data.crop_image = 1;
        data.median_filter = 1;
        data.magnification = 100;
        data.scale_image = 1;
        data.brightness_factor =0.2;
    case '0826_control09' % not good, no good photobleach
        data.path = strcat(root, 'bo_frap\0826_frap\noshear\st09\');
        data.output_path = strcat(data.path, 'output/');
        data.first_file = 'FRAP1.008';
        data.index_pattern = {'008', '%03d'};
        data.image_index = [8:10, 12:18];
        % photbleach between frames 10 and 11
        data.index_before = 1:3;
        data.index_after = 4:length(data.image_index);
        data.threshold = 0;
        data.photobleach_time = 0;
        data.subtract_background = 1;
        data.crop_image = 1;
        data.median_filter = 1;
        data.magnification = 100;
        data.scale_image = 1;
        data.brightness_factor =0.8;
    case '0826_control14' % not good, no good photobleach
        data.path = strcat(root, 'bo_frap\0826_frap\noshear\st14\');
        data.output_path = strcat(data.path, 'output/');
        data.first_file = 'FRAP1.013';
        data.index_pattern = {'013', '%03d'};
        data.image_index = [13:16, 17:23];
        % photbleach between frames 16 and 17
        data.index_before = 1:4;
        data.index_after = 5:length(data.image_index);
        data.threshold = 0;
        data.photobleach_time = 0;
        data.subtract_background = 1;
        data.crop_image = 1;
        data.median_filter = 1;
        data.magnification = 100;
        data.scale_image = 1;
        data.brightness_factor =0.8;
        %%%% 09/27/2014 Jie PB data %%%
    case '0916_cell3' % Before PB 0-4; After PB 
        % Not good. 
        data.path = strcat(root, 'Jie_Protocell/0916_FRAP/cell3/');
        data.output_path = strcat(data.path, 'output/');
        data.first_file = '20140916_1_FRAP Series104_t000_ch00.tif';
        data.index_pattern = {'t000', 't%03d'};
        data.image_index = [0:4, 5:11];
        % photbleach between frames 4 and 5
        data.index_before = 1:4;
        data.index_after = 6:length(data.image_index);
        %data.threshold = 0;
        data.photobleach_time = 0;
        data.subtract_background = 1;
        data.crop_image = 0;
        data.median_filter = 1;
        data.magnification = 300; %17.203pixels/micron
        data.scale_image = 128;
        data.brightness_factor =1.0;
        data.time_format = 2;
    case '0916_cell4' % Before PB 0-4; After PB 
        % Not good. 
        data.path = strcat(root, 'Jie_Protocell/0916_FRAP/cell4/');
        data.output_path = strcat(data.path, 'output/');
        data.first_file = '20140916_1_FRAP Series110_t000_ch00.tif';
        data.index_pattern = {'t000', 't%03d'};
        data.image_index = [0:4, 5:11];
        % photbleach between frames 4 and 5
        data.index_before = 1:5;
        data.index_after = 6:length(data.image_index);
        %data.threshold = 0;
        data.photobleach_time = 0;
        data.subtract_background = 1;
        data.crop_image = 0;
        data.median_filter = 1;
        data.magnification = 300; %17.203pixels/micron
        data.scale_image = 128;
        data.brightness_factor =1.0;
        data.time_format = 2;   
    case '0916_30nmcore3' % Before PB 0-4; After PB 
        data.path = strcat('/Users/shirleywu/Desktop/data/Jie/0916_FRAP/30nmcore3/');
        data.output_path = strcat(data.path, 'output/');
        data.first_file = '20140916_1_FRAP Series40_t000_ch00.tif';
        data.index_pattern = {'t000', 't%03d'};
        data.image_index = [0:4, 5:11];
        % photbleach between frames 4 and 5
        data.index_before = 1:5;
        data.index_after = 6:length(data.image_index);
        %data.threshold = 0;
        data.photobleach_time = 0;
        data.subtract_background = 1;
        data.crop_image = 0;
        data.median_filter = 1;
        data.magnification = 300; %17.203pixels/micron
        data.scale_image = 128;
        data.brightness_factor =1.0;
        data.time_format = 2;
    case '0916_30nmcore4' % Before PB 0-4; After PB 
        data.path = strcat('/Users/shirleywu/Desktop/0916_FRAP/30nmcore4/');
        data.output_path = strcat(data.path, 'output/');
        data.first_file = '20140916_1_FRAP Series46_t000_ch00.tif';
        data.index_pattern = {'t000', 't%03d'};
        data.image_index = [0:4, 5:11];
        % photbleach between frames 4 and 5
        data.index_before = 1:5;
        data.index_after = 6:length(data.image_index);
        %data.threshold = 0;
        data.photobleach_time = 0;
        data.subtract_background = 1;
        data.crop_image = 0;
        data.median_filter = 1;
        data.magnification = 300; %17.203pixels/micron
        data.scale_image = 128;
        data.brightness_factor =1.0;
        data.time_format = 2;
    case '0916_30nmcore5' % Before PB 0-4; After PB 
        data.path = strcat('/Users/shirleywu/Desktop/0916_FRAP/30nmcore5/');
        data.output_path = strcat(data.path, 'output/');
        data.first_file = '20140916_1_FRAP Series52_t000_ch00.tif';
        data.index_pattern = {'t000', 't%03d'};
        data.image_index = [0:4, 5:11];
        % photbleach between frames 4 and 5
        data.index_before = 1:5;
        data.index_after = 6:length(data.image_index);
        %data.threshold = 0;
        data.photobleach_time = 0;
        data.subtract_background = 1;
        data.crop_image = 0;
        data.median_filter = 1;
        data.magnification = 300; %17.203pixels/micron
        data.scale_image = 128;
        data.brightness_factor =1.0;
        data.time_format = 2;   
    case '0916_silicacell1' % Before PB 0-4; After PB 
        data.path = strcat('/Users/shirleywu/Desktop/0916_FRAP/silicacell1/');
        data.output_path = strcat(data.path, 'output/');
        data.first_file = '20140916_1_FRAP Series59_t000_ch00.tif';
        data.index_pattern = {'t000', 't%03d'};
        data.image_index = [0:4, 5:11];
        % photbleach between frames 4 and 5
        data.index_before = 1:5;
        data.index_after = 6:length(data.image_index);
        %data.threshold = 0;
        data.photobleach_time = 0;
        data.subtract_background = 1;
        data.crop_image = 0;
        data.median_filter = 1;
        data.magnification = 300; %17.203pixels/micron
        data.scale_image = 128;
        data.brightness_factor =1.0;
        data.time_format = 2;  
   case '0916_silicacell2' % Before PB 0-4; After PB 
        data.path = strcat('/Users/shirleywu/Desktop/0916_FRAP/silicacell2/');
        data.output_path = strcat(data.path, 'output/');
        data.first_file = '20140916_1_FRAP Series64_t000_ch00.tif';
        data.index_pattern = {'t000', 't%03d'};
        data.image_index = [0:4, 5:11];
        % photbleach between frames 4 and 5
        data.index_before = 1:5;
        data.index_after = 6:length(data.image_index);
        %data.threshold = 0;
        data.photobleach_time = 0;
        data.subtract_background = 1;
        data.crop_image = 0;
        data.median_filter = 1;
        data.magnification = 300; %17.203pixels/micron
        data.scale_image = 128;
        data.brightness_factor =1.0;
        data.time_format = 2;       
    case '0916_silicacell3' % Before PB 0-4; After PB 
        data.path = strcat('/Users/shirleywu/Desktop/0916_FRAP/silicacell3/');
        data.output_path = strcat(data.path, 'output/');
        data.first_file = '20140916_1_FRAP Series71_t000_ch00.tif';
        data.index_pattern = {'t000', 't%03d'};
        data.image_index = [0:4, 5:11];
        % photbleach between frames 4 and 5
        data.index_before = 1:5;
        data.index_after = 6:length(data.image_index);
        %data.threshold = 0;
        data.photobleach_time = 0;
        data.subtract_background = 1;
        data.crop_image = 0;
        data.median_filter = 1;
        data.magnification = 300; %17.203pixels/micron
        data.scale_image = 128;
        data.brightness_factor =1.0;
        data.time_format = 2; 
    case '0916_silicacell4' % Before PB 0-4; After PB 
        data.path = strcat('/Users/shirleywu/Desktop/0916_FRAP/silicacell4/');
        data.output_path = strcat(data.path, 'output/');
        data.first_file = '20140916_1_FRAP Series76_t000_ch00.tif';
        data.index_pattern = {'t000', 't%03d'};
        data.image_index = [0:4, 5:11];
        % photbleach between frames 4 and 5
        data.index_before = 1:5;
        data.index_after = 6:length(data.image_index);
        %data.threshold = 0;
        data.photobleach_time = 0;
        data.subtract_background = 1;
        data.crop_image = 0;
        data.median_filter = 1;
        data.magnification = 300; %17.203pixels/micron
        data.scale_image = 128;
        data.brightness_factor =1.0;
        data.time_format = 2; 
    case '0916_silicacell5' % Before PB 0-4; After PB 
        data.path = strcat('/Users/shirleywu/Desktop/0916_FRAP/silicacell5/');
        data.output_path = strcat(data.path, 'output/');
        data.first_file = '20140916_1_FRAP Series81_t000_ch00.tif';
        data.index_pattern = {'t000', 't%03d'};
        data.image_index = [0:4, 5:11];
        % photbleach between frames 4 and 5
        data.index_before = 1:5;
        data.index_after = 6:length(data.image_index);
        %data.threshold = 0;
        data.photobleach_time = 0;
        data.subtract_background = 1;
        data.crop_image = 0;
        data.median_filter = 1;
        data.magnification = 300; %17.203pixels/micron
        data.scale_image = 128;
        data.brightness_factor =1.0;
        data.time_format = 2; 
    case '1004_livecell10' % Before PB 0-4; After PB 
        data.path = strcat('/Users/shirleywu/Desktop/data/Jie/1004_FRAP/live cells/live cell 10/');
        data.output_path = strcat(data.path, 'output/');
        data.first_file = '4T1SrcFRAP_FRAP Series64_t000_ch00.tif';
        data.index_pattern = {'t000', 't%03d'};
        data.image_index = [0:4, 5:11];
        % photbleach between frames 4 and 5
        data.index_before = 1:5;
        data.index_after = 6:length(data.image_index);
        %data.threshold = 0;
        data.photobleach_time = 0;
        data.subtract_background = 1;
        data.crop_image = 0;
        data.median_filter = 1;
        data.magnification = 300; %17.203pixels/micron
        data.scale_image = 128;
        data.brightness_factor =1.0;
        data.time_format = 2;
   case '1004_livecell17' % Before PB 0-2; After PB 
        data.path = strcat('/Users/shirleywu/Desktop/data/Jie/1004_FRAP/live cells/live cell 17/');
        data.output_path = strcat(data.path, 'output/');
        data.first_file = '4T1SrcFRAP_FRAP Series88_t000_ch00.tif';
        data.index_pattern = {'t000', 't%03d'};
        data.image_index = [0:2, 3:9];
        % photbleach between frames 4 and 5
        data.index_before = 1:3;
        data.index_after = 4:length(data.image_index);
        %data.threshold = 0;
        data.photobleach_time = 0;
        data.subtract_background = 1;
        data.crop_image = 0;
        data.median_filter = 1;
        data.magnification = 300; %17.203pixels/micron
        data.scale_image = 128;
        data.brightness_factor =1.0;
        data.time_format = 2;
    case '1004_livecell19' % Before PB 0-2; After PB 
        data.path = strcat('/Users/shirleywu/Desktop/data/Jie/1004_FRAP/live cells/live cell 19/');
        data.output_path = strcat(data.path, 'output/');
        data.first_file = 'FRAP_FRAP Series182_t000_ch00.tif';
        data.index_pattern = {'t000', 't%03d'};
        data.image_index = [0:2, 3:9];
        % photbleach between frames 4 and 5
        data.index_before = 1:3;
        data.index_after = 4:length(data.image_index);
        %data.threshold = 0;
        data.photobleach_time = 0;
        data.subtract_background = 1;
        data.crop_image = 0;
        data.median_filter = 1;
        data.magnification = 300; %17.203pixels/micron
        data.scale_image = 128;
        data.brightness_factor =1.0;
        data.time_format = 2;
    case '1004_livecell20' % Before PB 0-2; After PB 
        data.path = strcat('/Users/shirleywu/Desktop/data/Jie/1004_FRAP/live cells/live cell 20/');
        data.output_path = strcat(data.path, 'output/');
        data.first_file = '4T1SrcFRAP_FRAP Series105_t000_ch00.tif';
        data.index_pattern = {'t000', 't%03d'};
        data.image_index = [0:2, 3:9];
        % photbleach between frames 4 and 5
        data.index_before = 1:3;
        data.index_after = 4:length(data.image_index);
        %data.threshold = 0;
        data.photobleach_time = 0;
        data.subtract_background = 1;
        data.crop_image = 0;
        data.median_filter = 1;
        data.magnification = 300; %17.203pixels/micron
        data.scale_image = 128;
        data.brightness_factor =1.0;
        data.time_format = 2;
   case '1004_livecell26' % Before PB 0-2; After PB 
        data.path = strcat('/Users/shirleywu/Desktop/data/Jie/1004_FRAP/live cells/live cell 26/');
        data.output_path = strcat(data.path, 'output/');
        data.first_file = 'FRAP_FRAP Series245_t000_ch00.tif';
        data.index_pattern = {'t000', 't%03d'};
        data.image_index = [0:2, 3:9];
        % photbleach between frames 4 and 5
        data.index_before = 1:3;
        data.index_after = 4:length(data.image_index);
        %data.threshold = 0;
        data.photobleach_time = 0;
        data.subtract_background = 1;
        data.crop_image = 0;
        data.median_filter = 1;
        data.magnification = 300; %17.203pixels/micron
        data.scale_image = 128;
        data.brightness_factor =1.0;
        data.time_format = 2;
    case '1004_livecell27' % Before PB 0-2; After PB 
        data.path = strcat('/Users/shirleywu/Desktop/data/Jie/1004_FRAP/live cells/live cell 27/');
        data.output_path = strcat(data.path, 'output/');
        data.first_file = '4T1SrcFRAP_FRAP Series142_t000_ch00.tif';
        data.index_pattern = {'t000', 't%03d'};
        data.image_index = [0:2, 3:9];
        % photbleach between frames 4 and 5
        data.index_before = 1:3;
        data.index_after = 4:length(data.image_index);
        %data.threshold = 0;
        data.photobleach_time = 0;
        data.subtract_background = 1;
        data.crop_image = 1;
        data.median_filter = 1;
        data.magnification = 300; %17.203pixels/micron
        data.scale_image = 128;
        data.brightness_factor =1.0;
        data.time_format = 2;
    case '1004_livecell29' % Before PB 0-2; After PB 
        data.path = strcat('/Users/shirleywu/Desktop/data/Jie/1004_FRAP/live cells/live cell 29/');
        data.output_path = strcat(data.path, 'output/');
        data.first_file = '4T1SrcFRAP_FRAP Series154_t000_ch00.tif';
        data.index_pattern = {'t000', 't%03d'};
        data.image_index = [0:2, 3:9];
        % photbleach between frames 4 and 5
        data.index_before = 1:3;
        data.index_after = 4:length(data.image_index);
        %data.threshold = 0;
        data.photobleach_time = 0;
        data.subtract_background = 1;
        data.crop_image = 0;
        data.median_filter = 1;
        data.magnification = 300; %17.203pixels/micron
        data.scale_image = 128;
        data.brightness_factor =1.0;
        data.time_format = 2;
    case '1004_livecell31' % Before PB 0-2; After PB 
        data.path = strcat('/Users/shirleywu/Desktop/data/Jie/1004_FRAP/live cells/live cell 31/');
        data.output_path = strcat(data.path, 'output/');
        data.first_file = 'FRAP_FRAP Series275_t000_ch00.tif';
        data.index_pattern = {'t000', 't%03d'};
        data.image_index = [0:2, 3:9];
        % photbleach between frames 4 and 5
        data.index_before = 1:3;
        data.index_after = 4:length(data.image_index);
        %data.threshold = 0;
        data.photobleach_time = 0;
        data.subtract_background = 1;
        data.crop_image = 0;
        data.median_filter = 1;
        data.magnification = 300; %17.203pixels/micron
        data.scale_image = 128;
        data.brightness_factor =1.0;
        data.time_format = 2;
    case '1004_silicacell00' % Before PB 0-2; After PB 
        data.path = strcat('/Users/shirleywu/Desktop/data/Jie/1004_FRAP/silica cells/sc 00/');
        data.output_path = strcat(data.path, 'output/');
        data.first_file = 'silica cell_FRAP Series04_t000_ch00.tif';
        data.index_pattern = {'t000', 't%03d'};
        data.image_index = [0:2, 3:9];
        % photbleach between frames 4 and 5
        data.index_before = 1:3;
        data.index_after = 4:length(data.image_index);
        %data.threshold = 0;
        data.photobleach_time = 0;
        data.subtract_background = 1;
        data.crop_image = 1;
        data.median_filter = 1;
        data.magnification = 300; %17.203pixels/micron
        data.scale_image = 128;
        data.brightness_factor =1.0;
        data.time_format = 2;
   case '1004_silicacell01' % Before PB 0-2; After PB 
        data.path = strcat('/Users/shirleywu/Desktop/data/Jie/1004_FRAP/silica cells/sc 01/');
        data.output_path = strcat(data.path, 'output/');
        data.first_file = 'silica cell_FRAP Series10_t000_ch00.tif';
        data.index_pattern = {'t000', 't%03d'};
        data.image_index = [0:2, 3:9];
        % photbleach between frames 4 and 5
        data.index_before = 1:3;
        data.index_after = 4:length(data.image_index);
        %data.threshold = 0;
        data.photobleach_time = 0;
        data.subtract_background = 1;
        data.crop_image = 0;
        data.median_filter = 1;
        data.magnification = 300; %17.203pixels/micron
        data.scale_image = 128;
        data.brightness_factor =1.0;
        data.time_format = 2;
    case '1004_silicacell03' % Before PB 0-2; After PB 
        data.path = strcat('/Users/shirleywu/Desktop/data/Jie/1004_FRAP/silica cells/sc 03/');
        data.output_path = strcat(data.path, 'output/');
        data.first_file = 'silica cell_FRAP Series22_t000_ch00.tif';
        data.index_pattern = {'t000', 't%03d'};
        data.image_index = [0:2, 3:9];
        % photbleach between frames 4 and 5
        data.index_before = 1:3;
        data.index_after = 4:length(data.image_index);
        %data.threshold = 0;
        data.photobleach_time = 0;
        data.subtract_background = 1;
        data.crop_image = 0;
        data.median_filter = 1;
        data.magnification = 300; %17.203pixels/micron
        data.scale_image = 128;
        data.brightness_factor =1.0;
        data.time_format = 2;
    case '1004_silicacell06' % Before PB 0-2; After PB 
        data.path = strcat('/Users/shirleywu/Desktop/data/Jie/1004_FRAP/silica cells/sc 06/');
        data.output_path = strcat(data.path, 'output/');
        data.first_file = 'silica cell_FRAP Series38_t000_ch00.tif';
        data.index_pattern = {'t000', 't%03d'};
        data.image_index = [0:2, 3:9];
        % photbleach between frames 4 and 5
        data.index_before = 1:3;
        data.index_after = 4:length(data.image_index);
        %data.threshold = 0;
        data.photobleach_time = 0;
        data.subtract_background = 1;
        data.crop_image = 0;
        data.median_filter = 1;
        data.magnification = 300; %17.203pixels/micron
        data.scale_image = 128;
        data.brightness_factor =1.0;
        data.time_format = 2;
    case '1004_silicacell07' % Before PB 0-2; After PB 
        data.path = strcat('/Users/shirleywu/Desktop/data/Jie/1004_FRAP/silica cells/sc 07/');
        data.output_path = strcat(data.path, 'output/');
        data.first_file = 'silica cell_FRAP Series44_t000_ch00.tif';
        data.index_pattern = {'t000', 't%03d'};
        data.image_index = [0:2, 3:9];
        % photbleach between frames 4 and 5
        data.index_before = 1:3;
        data.index_after = 4:length(data.image_index);
        %data.threshold = 0;
        data.photobleach_time = 0;
        data.subtract_background = 1;
        data.crop_image = 0;
        data.median_filter = 1;
        data.magnification = 300; %17.203pixels/micron
        data.scale_image = 128;
        data.brightness_factor =1.0;
        data.time_format = 2;
    case '1004_silicacell18' % Before PB 0-2; After PB 
        data.path = strcat('/Users/shirleywu/Desktop/data/Jie/1004_FRAP/silica cells/sc 18/');
        data.output_path = strcat(data.path, 'output/');
        data.first_file = 'silica cell_FRAP Series107_t000_ch00.tif';
        data.index_pattern = {'t000', 't%03d'};
        data.image_index = [0:2, 3:9];
        % photbleach between frames 4 and 5
        data.index_before = 1:3;
        data.index_after = 4:length(data.image_index);
        %data.threshold = 0;
        data.photobleach_time = 0;
        data.subtract_background = 1;
        data.crop_image = 1;
        data.median_filter = 1;
        data.magnification = 300; %17.203pixels/micron
        data.scale_image = 128;
        data.brightness_factor =1.0;
        data.time_format = 2;
    case '1004_silicacell19' % Before PB 0-2; After PB 
        data.path = strcat('/Users/shirleywu/Desktop/data/Jie/1004_FRAP/silica cells/sc 19/');
        data.output_path = strcat(data.path, 'output/');
        data.first_file = 'silica cell_FRAP Series113_t000_ch00.tif';
        data.index_pattern = {'t000', 't%03d'};
        data.image_index = [0:2, 3:9];
        % photbleach between frames 4 and 5
        data.index_before = 1:3;
        data.index_after = 4:length(data.image_index);
        %data.threshold = 0;
        data.photobleach_time = 0;
        data.subtract_background = 1;
        data.crop_image = 0;
        data.median_filter = 1;
        data.magnification = 300; %17.203pixels/micron
        data.scale_image = 128;
        data.brightness_factor =1.0;
        data.time_format = 2;
    case '1004_silicacell20' % Before PB 0-2; After PB 
        data.path = strcat('/Users/shirleywu/Desktop/data/Jie/1004_FRAP/silica cells/sc 20/');
        data.output_path = strcat(data.path, 'output/');
        data.first_file = 'silica cell_FRAP Series120_t000_ch00.tif';
        data.index_pattern = {'t000', 't%03d'};
        data.image_index = [0:2, 3:9];
        % photbleach between frames 4 and 5
        data.index_before = 1:3;
        data.index_after = 4:length(data.image_index);
        %data.threshold = 0;
        data.photobleach_time = 0;
        data.subtract_background = 1;
        data.crop_image = 0;
        data.median_filter = 1;
        data.magnification = 300; %17.203pixels/micron
        data.scale_image = 128;
        data.brightness_factor =1.0;
        data.time_format = 2;
    case '1004_core00' % Before PB 0-2; After PB 
        data.path = strcat('/Users/shirleywu/Desktop/data/Jie/1004_FRAP/30nm synthetic cores/core 00/');
        data.output_path = strcat(data.path, 'output/');
        data.first_file = '30nmcores_FRAP Series02_t000_ch00.tif';
        data.index_pattern = {'t000', 't%03d'};
        data.image_index = [0:2, 3:9];
        % photbleach between frames 4 and 5
        data.index_before = 1:3;
        data.index_after = 4:length(data.image_index);
        %data.threshold = 0;
        data.photobleach_time = 0;
        data.subtract_background = 1;
        data.crop_image = 0;
        data.median_filter = 1;
        data.magnification = 300; %17.203pixels/micron
        data.scale_image = 128;
        data.brightness_factor =1.0;
        data.time_format = 2;
    case '1004_core01' % Before PB 0-2; After PB 
        data.path = strcat('/Users/shirleywu/Desktop/data/Jie/1004_FRAP/30nm synthetic cores/core 01/');
        data.output_path = strcat(data.path, 'output/');
        data.first_file = '30nmcores_FRAP Series08_t000_ch00.tif';
        data.index_pattern = {'t000', 't%03d'};
        data.image_index = [0:2, 3:9];
        % photbleach between frames 4 and 5
        data.index_before = 1:3;
        data.index_after = 4:length(data.image_index);
        %data.threshold = 0;
        data.photobleach_time = 0;
        data.subtract_background = 1;
        data.crop_image = 0;
        data.median_filter = 1;
        data.magnification = 300; %17.203pixels/micron
        data.scale_image = 128;
        data.brightness_factor =1.0;
        data.time_format = 2;
    case '1004_core03' % Before PB 0-2; After PB 
        data.path = strcat('/Users/shirleywu/Desktop/data/Jie/1004_FRAP/30nm synthetic cores/core 03/');
        data.output_path = strcat(data.path, 'output/');
        data.first_file = '30nmcores_FRAP Series20_t000_ch00.tif';
        data.index_pattern = {'t000', 't%03d'};
        data.image_index = [0:2, 3:9];
        % photbleach between frames 4 and 5
        data.index_before = 1:3;
        data.index_after = 4:length(data.image_index);
        %data.threshold = 0;
        data.photobleach_time = 0;
        data.subtract_background = 1;
        data.crop_image = 0;
        data.median_filter = 1;
        data.magnification = 300; %17.203pixels/micron
        data.scale_image = 128;
        data.brightness_factor =1.0;
        data.time_format = 2;
    case '1004_core04' % Before PB 0-2; After PB 
        data.path = strcat('/Users/shirleywu/Desktop/data/Jie/1004_FRAP/30nm synthetic cores/core 04/');
        data.output_path = strcat(data.path, 'output/');
        data.first_file = '30nmcores_FRAP Series03_t000_ch00.tif';
        data.index_pattern = {'t000', 't%03d'};
        data.image_index = [0:2, 3:9];
        % photbleach between frames 4 and 5
        data.index_before = 1:3;
        data.index_after = 4:length(data.image_index);
        %data.threshold = 0;
        data.photobleach_time = 0;
        data.subtract_background = 1;
        data.crop_image = 0;
        data.median_filter = 1;
        data.magnification = 300; %17.203pixels/micron
        data.scale_image = 128;
        data.brightness_factor =1.0;
        data.time_format = 2;
    case '1004_core05' % Before PB 0-2; After PB 
        data.path = strcat('/Users/shirleywu/Desktop/data/Jie/1004_FRAP/30nm synthetic cores/core 05/');
        data.output_path = strcat(data.path, 'output/');
        data.first_file = '30nmcores_FRAP Series16_t000_ch00.tif';
        data.index_pattern = {'t000', 't%03d'};
        data.image_index = [0:2, 3:9];
        % photbleach between frames 4 and 5
        data.index_before = 1:3;
        data.index_after = 4:length(data.image_index);
        %data.threshold = 0;
        data.photobleach_time = 0;
        data.subtract_background = 1;
        data.crop_image = 0;
        data.median_filter = 1;
        data.magnification = 300; %17.203pixels/micron
        data.scale_image = 128;
        data.brightness_factor =1.0;
        data.time_format = 2;
    case '1004_core06' % Before PB 0-2; After PB 
        data.path = strcat('/Users/shirleywu/Desktop/data/Jie/1004_FRAP/30nm synthetic cores/core 06/');
        data.output_path = strcat(data.path, 'output/');
        data.first_file = '30nmcores_FRAP Series22_t000_ch00.tif';
        data.index_pattern = {'t000', 't%03d'};
        data.image_index = [0:2, 3:9];
        % photbleach between frames 4 and 5
        data.index_before = 1:3;
        data.index_after = 4:length(data.image_index);
        %data.threshold = 0;
        data.photobleach_time = 0;
        data.subtract_background = 1;
        data.crop_image = 0;
        data.median_filter = 1;
        data.magnification = 300; %17.203pixels/micron
        data.scale_image = 128;
        data.brightness_factor =1.0;
        data.time_format = 2;
    case '1004_core07' % Before PB 0-2; After PB 
        data.path = strcat('/Users/shirleywu/Desktop/data/Jie/1004_FRAP/30nm synthetic cores/core 07/');
        data.output_path = strcat(data.path, 'output/');
        data.first_file = '30nmcores_FRAP Series28_t000_ch00.tif';
        data.index_pattern = {'t000', 't%03d'};
        data.image_index = [0:2, 3:9];
        % photbleach between frames 4 and 5
        data.index_before = 1:3;
        data.index_after = 4:length(data.image_index);
        %data.threshold = 0;
        data.photobleach_time = 0;
        data.subtract_background = 1;
        data.crop_image = 0;
        data.median_filter = 1;
        data.magnification = 300; %17.203pixels/micron
        data.scale_image = 128;
        data.brightness_factor =1.0;
        data.time_format = 2;
    case '1004_core09' % Before PB 0-2; After PB 
        data.path = strcat('/Users/shirleywu/Desktop/data/Jie/1004_FRAP/30nm synthetic cores/core 09/');
        data.output_path = strcat(data.path, 'output/');
        data.first_file = '30nmcores_FRAP Series59_t000_ch00.tif';
        data.index_pattern = {'t000', 't%03d'};
        data.image_index = [0:2, 3:9];
        % photbleach between frames 4 and 5
        data.index_before = 1:3;
        data.index_after = 4:length(data.image_index);
        %data.threshold = 0;
        data.photobleach_time = 0;
        data.subtract_background = 1;
        data.crop_image = 0;
        data.median_filter = 1;
        data.magnification = 300; %17.203pixels/micron
        data.scale_image = 128;
        data.brightness_factor =1.0;
        data.time_format = 2;
    case '1101_cell 04' % Before PB 0-4; After PB 
        data.path = strcat(root,'Jie/1101_FRAP/cell 04/');
        data.output_path = strcat(data.path, 'output/');
        data.first_file = 'FRAPlive_cell04_t000_ch00.tif';
        data.index_pattern = {'t000', 't%03d'};
        data.image_index = [0:4, 5:11];
        % photbleach between frames 4 and 5
        data.index_before = 1:5;
        data.index_after = 6:length(data.image_index);
        %data.threshold = 0;
        data.photobleach_time = 0;
        data.subtract_background = 1;
        data.crop_image = 0;
        data.median_filter = 1;
        data.magnification = 150; %8.602pixels/micron
        data.scale_image = 128;
        data.brightness_factor =1.0;
        data.time_format = 2;
    case '1101_cell 06' % Before PB 0-4; After PB 
        data.path = strcat(root,'Jie/1101_FRAP/cell 06/');
        data.output_path = strcat(data.path, 'output/');
        data.first_file = 'FRAPlive_cell06_t000_ch00.tif';
        data.index_pattern = {'t000', 't%03d'};
        data.image_index = [0:4, 5:11];
        % photbleach between frames 4 and 5
        data.index_before = 1:5;
        data.index_after = 6:length(data.image_index);
        %data.threshold = 0;
        data.photobleach_time = 0;
        data.subtract_background = 1;
        data.crop_image = 0;
        data.median_filter = 1;
        data.magnification = 150; %8.602pixels/micron
        data.scale_image = 128;
        data.brightness_factor =1.0;
        data.time_format = 2;
    case '1101_cell 07' % Before PB 0-4; After PB 
        data.path = strcat(root,'Jie/1101_FRAP/cell 07/');
        data.output_path = strcat(data.path, 'output/');
        data.first_file = 'FRAPlive_cell07_t000_ch00.tif';
        data.index_pattern = {'t000', 't%03d'};
        data.image_index = [0:4, 5:11];
        % photbleach between frames 4 and 5
        data.index_before = 1:5;
        data.index_after = 6:length(data.image_index);
        %data.threshold = 0;
        data.photobleach_time = 0;
        data.subtract_background = 1;
        data.crop_image = 0;
        data.median_filter = 1;
        data.magnification = 150; %8.602pixels/micron
        data.scale_image = 128;
        data.brightness_factor =1.0;
        data.time_format = 2;
    case '1101_cell 09' % Before PB 0-4; After PB 
        data.path = strcat(root,'Jie/1101_FRAP/cell 09/');
        data.output_path = strcat(data.path, 'output/');
        data.first_file = 'FRAPlive_cell09_t000_ch00.tif';
        data.index_pattern = {'t000', 't%03d'};
        data.image_index = [0:4, 5:11];
        % photbleach between frames 4 and 5
        data.index_before = 1:5;
        data.index_after = 6:length(data.image_index);
        %data.threshold = 0;
        data.photobleach_time = 0;
        data.subtract_background = 1;
        data.crop_image = 0;
        data.median_filter = 1;
        data.magnification = 150; %8.602pixels/micron
        data.scale_image = 128;
        data.brightness_factor =1.0;
        data.time_format = 2;
    case '1101_cell 15' % Before PB 0-4; After PB 
        data.path = strcat(root,'Jie/1101_FRAP/cell 15/');
        data.output_path = strcat(data.path, 'output/');
        data.first_file = 'FRAPlive_cell15_t000_ch00.tif';
        data.index_pattern = {'t000', 't%03d'};
        data.image_index = [0:4, 5:11];
        % photbleach between frames 4 and 5
        data.index_before = 1:5;
        data.index_after = 6:length(data.image_index);
        %data.threshold = 0;
        data.photobleach_time = 0;
        data.subtract_background = 1;
        data.crop_image = 0;
        data.median_filter = 1;
        data.magnification = 150; %8.602pixels/micron
        data.scale_image = 128;
        data.brightness_factor =1.0;
        data.time_format = 2;
    case '1101_cell 21' % Before PB 0-4; After PB 
        data.path = strcat(root,'Jie/1101_FRAP/cell 21/');
        data.output_path = strcat(data.path, 'output/');
        data.first_file = 'FRAPlive_cell21_t000_ch00.tif';
        data.index_pattern = {'t000', 't%03d'};
        data.image_index = [0:4, 5:11];
        % photbleach between frames 4 and 5
        data.index_before = 1:5;
        data.index_after = 6:length(data.image_index);
        %data.threshold = 0;
        data.photobleach_time = 0;
        data.subtract_background = 1;
        data.crop_image = 0;
        data.median_filter = 1;
        data.magnification = 150; %8.602pixels/micron
        data.scale_image = 128;
        data.brightness_factor =1.0;
        data.time_format = 2;
    case '1101_cell 23' % Before PB 0-4; After PB 7
        data.path = strcat(root,'Jie/1101_FRAP/cell 23/');
        data.output_path = strcat(data.path, 'output/');
        data.first_file = 'FRAPlive_cell23_t000_ch00.tif';
        data.index_pattern = {'t000', 't%03d'};
        data.image_index = [0:4, 5:11];
        % photbleach between frames 4 and 5
        data.index_before = 1:5;
        data.index_after = 6:length(data.image_index);
        %data.threshold = 0;
        data.photobleach_time = 0;
        data.subtract_background = 1;
        data.crop_image = 0;
        data.median_filter = 1;
        data.magnification = 150; %8.602pixels/micron
        data.scale_image = 128;
        data.brightness_factor =1.0;
        data.time_format = 2;
    case '1113_silica01' % Before PB 0-4; After PB 7
        data.path = strcat(root,'Jie/2014 11 13 in vitro in silica imaging/30nm01/');
        data.output_path = strcat(data.path, 'output/');
        data.first_file = 'invitroinsilica_30nm01_t000_ch00.tif';
        data.index_pattern = {'t000', 't%03d'};
        data.image_index = [0:4, 5:11];
        % photbleach between frames 4 and 5
        data.index_before = 1:5;
        data.index_after = 6:length(data.image_index);
        %data.threshold = 0;
        data.photobleach_time = 0;
        data.subtract_background = 1;
        data.crop_image = 0;
        data.median_filter = 1;
        data.magnification = 150; %8.602pixels/micron
        data.scale_image = 128;
        data.brightness_factor =1.0;
        data.time_format = 2;
    case '1113_silica02' % Before PB 0-4; After PB 7
        data.path = strcat(root,'Jie/2014 11 13 in vitro in silica imaging/30nm02/');
        data.output_path = strcat(data.path, 'output/');
        data.first_file = 'invitroinsilica_30nm02_t000_ch00.tif';
        data.index_pattern = {'t000', 't%03d'};
        data.image_index = [0:4, 5:11];
        % photbleach between frames 4 and 5
        data.index_before = 1:5;
        data.index_after = 6:length(data.image_index);
        %data.threshold = 0;
        data.photobleach_time = 0;
        data.subtract_background = 1;
        data.crop_image = 0;
        data.median_filter = 1;
        data.magnification = 150; %8.602pixels/micron
        data.scale_image = 128;
        data.brightness_factor =1.0;
        data.time_format = 2;
    case '1113_invitro01' % Before PB 0-4; After PB 7
        data.path = strcat(root,'Jie/2014 11 13 in vitro in silica imaging/invitro01/');
        data.output_path = strcat(data.path, 'output/');
        data.first_file = 'invitroinsilica_invitro01_t000_ch00.tif';
        data.index_pattern = {'t000', 't%03d'};
        data.image_index = [0:4, 5:11];
        % photbleach between frames 4 and 5
        data.index_before = 1:5;
        data.index_after = 6:length(data.image_index);
        %data.threshold = 0;
        data.photobleach_time = 0;
        data.subtract_background = 0; % for protein solution, no subtraction of bg.
        data.crop_image = 0;
        data.median_filter = 1;  
        data.magnification = 150; %8.602pixels/micron
        data.scale_image = 1;% for protein solution, make the multiplication factor 1 instead of 128.
        data.brightness_factor =1.0;
        data.time_format = 2;
    case '1113_invitro03' % Before PB 0-4; After PB 7
        data.path = strcat(root,'Jie/2014 11 13 in vitro in silica imaging/invitro03/');
        data.output_path = strcat(data.path, 'output/');
        data.first_file = 'invitroinsilica_invitro03_t000_ch00.tif';
        data.index_pattern = {'t000', 't%03d'};
        data.image_index = [0:4, 5:11];
        % photbleach between frames 4 and 5
        data.index_before = 1:5;
        data.index_after = 6:length(data.image_index);
        %data.threshold = 0;
        data.photobleach_time = 0;
        data.subtract_background = 0;
        data.crop_image = 0;
        data.median_filter = 1;
        data.magnification = 150; %8.602pixels/micron
        data.scale_image = 1;
        data.brightness_factor =1.0;
        data.time_format = 2;
    case '1113_invitro04' % Before PB 0-4; After PB 7
        data.path = strcat(root,'Jie/2014 11 13 in vitro in silica imaging/invitro04/');
        data.output_path = strcat(data.path, 'output/');
        data.first_file = 'invitroinsilica_invitro04_t000_ch00.tif';
        data.index_pattern = {'t000', 't%03d'};
        data.image_index = [0:4, 5:11];
        % photbleach between frames 4 and 5
        data.index_before = 1:5;
        data.index_after = 6:length(data.image_index);
        %data.threshold = 0;
        data.photobleach_time = 0;
        data.subtract_background = 0;
        data.crop_image = 0;
        data.median_filter = 1;
        data.magnification = 150; %8.602pixels/micron
        data.scale_image = 1;
        data.brightness_factor =1.0;
        data.time_format = 2;
    case '1113_invitro05' % Before PB 0-4; After PB 7
        data.path = strcat(root,'Jie/2014 11 13 in vitro in silica imaging/invitro05/');
        data.output_path = strcat(data.path, 'output/');
        data.first_file = 'invitroinsilica_invitro05_t000_ch00.tif';
        data.index_pattern = {'t000', 't%03d'};
        data.image_index = [0:4, 5:11];
        % photbleach between frames 4 and 5
        data.index_before = 1:5;
        data.index_after = 6:length(data.image_index);
        %data.threshold = 0;
        data.photobleach_time = 0;
        data.subtract_background = 0;
        data.crop_image = 0;
        data.median_filter = 1;
        data.magnification = 150; %8.602pixels/micron
        data.scale_image = 1;
        data.brightness_factor =1.0;
        data.time_format = 2;
    case '1113_invitro06' % Before PB 0-4; After PB 7
        data.path = strcat(root,'Jie/2014 11 13 in vitro in silica imaging/invitro06/');
        data.output_path = strcat(data.path, 'output/');
        data.first_file = 'invitroinsilica_invitro06_t000_ch00.tif';
        data.index_pattern = {'t000', 't%03d'};
        data.image_index = [0:4, 5:11];
        % photbleach between frames 4 and 5
        data.index_before = 1:5;
        data.index_after = 6:length(data.image_index);
        %data.threshold = 0;
        data.photobleach_time = 0;
        data.subtract_background = 0;
        data.crop_image = 0;
        data.median_filter = 1;
        data.magnification = 150; %8.602pixels/micron
        data.scale_image = 1;
        data.brightness_factor =1.0;
        data.time_format = 2;
end

if ~isfield(data, 'path')
    data = init_data_0703_2014(cell_name);
end

% convert from file name to file prefix here
if ~isfield(data, 'prefix') && isfield(data,'first_file')
    dot_i = regexpi(data.first_file, '\.');
    data.prefix = data.first_file(1:dot_i-2);
end

return;





