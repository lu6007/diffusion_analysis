% function data = diffusion_init_data(cell_name)
% Initialize cell data

% Copyright: Shaoying Lu and Yingxiao Wang 2014
function data = diffusion_init_data_1023(cell_name)
% root='/Users/shirleywu/Desktop/data/';
% root = 'D:/sof/data/';
% root = '/Volumes/KathyWD2TB/data/2016/rongxue/';
root = '/Users/kathylu/Documents/data/2018/';
data.cell_name = cell_name;
switch cell_name
    case 'yf5_328_1018'
        data.path = strcat(root, 'rongxue_1023/yf328_1018/5/');
        data.output_path = strcat(data.path, 'output/');
        data.first_file = '5_w1Imaging RFP_s1_t1.TIF';
        data.index_pattern = {'t1', 't%d'};
        data.image_index = (1:10);
        data.index_before = 1:4; 
        data.index_after = 5:10;
        data.photobleach_time = 0; 
        data.dt = 10.0/60; %min
        data.subtract_background = 1;
        data.crop_image = 1;
        data.median_filter = 1;
        data.magnification = 60;
        data.time_format = 2; 
        data.brightness_factor = 1.4;
    case 'yf4_328_1018'
        data.path = strcat(root, 'rongxue_1023/yf328_1018/4/');
        data.output_path = strcat(data.path, 'output/');
        data.first_file = '4_w1Imaging RFP_s1_t1.TIF';
        data.index_pattern = {'t1', 't%d'};
        data.image_index = (1:10);
        data.index_before = 1:4; 
        data.index_after = 5:10;
        data.photobleach_time = 0; 
        data.dt = 10.0/60; %min
        data.subtract_background = 1;
        data.crop_image = 1;
        data.median_filter = 1;
        data.magnification = 60;
        data.time_format = 2; 
        data.brightness_factor = 1.2;
    case 'yf3_328_1018'
        data.path = strcat(root, 'rongxue_1023/yf328_1018/3/');
        data.output_path = strcat(data.path, 'output/');
        data.first_file = '3_w1Imaging RFP_s1_t1.TIF';
        data.index_pattern = {'t1', 't%d'};
        data.image_index = (1:10);
        data.index_before = 1:4; 
        data.index_after = 5:10;
        data.photobleach_time = 0; 
        data.dt = 10.0/60; %min
        data.subtract_background = 1;
        data.crop_image = 1;
        data.median_filter = 1;
        data.magnification = 60;
        data.time_format = 2; 
        data.brightness_factor = 0.95;
    case 'yf2_328_1018' % not good
        data.path = strcat(root, 'rongxue_1023/yf328_1018/2/');
        data.output_path = strcat(data.path, 'output/');
        data.first_file = '2_w1Imaging RFP_s1_t1.TIF';
        data.index_pattern = {'t1', 't%d'};
        data.image_index = (1:10);
        data.index_before = 1:4; 
        data.index_after = 5:10;
        data.photobleach_time = 0; 
        data.dt = 10.0/60; %min
        data.subtract_background = 1;
        data.crop_image = 1;
        data.median_filter = 1;
        data.magnification = 60;
        data.time_format = 2; 
        data.brightness_factor = 1.5;
    case 'yf1_328_1018' 
        data.path = strcat(root, 'rongxue_1023/yf328_1018/1/');
        data.output_path = strcat(data.path, 'output/');
        data.first_file = '1_w1Imaging RFP_s1_t1.TIF';
        data.index_pattern = {'t1', 't%d'};
        data.image_index = (1:10);
        data.index_before = 1:4; 
        data.index_after = 5:10;
        data.photobleach_time = 0; 
        data.dt = 10.0/60; %min
        data.subtract_background = 1;
        data.crop_image = 1;
        data.median_filter = 1;
        data.magnification = 60;
        data.time_format = 2; 
        data.brightness_factor = 1.2;
    case 'yf7_1023' 
        data.path = strcat(root, 'rongxue_1023/yf_1023/7/');
        data.output_path = strcat(data.path, 'output/');
        data.first_file = '7_w1Imaging RFP_s1_t1.TIF';
        data.index_pattern = {'t1', 't%d'};
        data.image_index = (1:10);
        data.index_before = 1:4; 
        data.index_after = 5:10;
        data.photobleach_time = 0; 
        data.dt = 10.0/60; %min
        data.subtract_background = 1;
        data.crop_image = 1;
        data.median_filter = 1;
        data.magnification = 60;
        data.time_format = 2; 
        data.brightness_factor = 1.7;
    case 'yf6_1023' 
        data.path = strcat(root, 'rongxue_1023/yf_1023/6/');
        data.output_path = strcat(data.path, 'output/');
        data.first_file = '6_w1Imaging RFP_s1_t1.TIF';
        data.index_pattern = {'t1', 't%d'};
        data.image_index = (1:10);
        data.index_before = 1:4; 
        data.index_after = 5:10;
        data.photobleach_time = 0; 
        data.dt = 10.0/60; %min
        data.subtract_background = 1;
        data.crop_image = 1;
        data.median_filter = 1;
        data.magnification = 60;
        data.time_format = 2; 
        data.brightness_factor = 1.5;
    case 'yf5_1023' 
        data.path = strcat(root, 'rongxue_1023/yf_1023/5/');
        data.output_path = strcat(data.path, 'output/');
        data.first_file = '5_w1Imaging RFP_s1_t1.TIF';
        data.index_pattern = {'t1', 't%d'};
        data.image_index = (1:10);
        data.index_before = 1:4; 
        data.index_after = 5:10;
        data.photobleach_time = 0; 
        data.dt = 10.0/60; %min
        data.subtract_background = 1;
        data.crop_image = 1;
        data.median_filter = 1;
        data.magnification = 60;
        data.time_format = 2; 
        data.brightness_factor = 1.35;
    case 'yf4_1023' 
        data.path = strcat(root, 'rongxue_1023/yf_1023/4/');
        data.output_path = strcat(data.path, 'output/');
        data.first_file = '4_w1Imaging RFP_s1_t1.TIF';
        data.index_pattern = {'t1', 't%d'};
        data.image_index = (1:10);
        data.index_before = 1:4; 
        data.index_after = 5:10;
        data.photobleach_time = 0; 
        data.dt = 10.0/60; %min
        data.subtract_background = 1;
        data.crop_image = 1;
        data.median_filter = 1;
        data.magnification = 60;
        data.time_format = 2; 
        data.brightness_factor = 1.4;
    case 'yf3_1023' 
        data.path = strcat(root, 'rongxue_1023/yf_1023/3/');
        data.output_path = strcat(data.path, 'output/');
        data.first_file = '3_w1Imaging RFP_s1_t1.TIF';
        data.index_pattern = {'t1', 't%d'};
        data.image_index = (1:10);
        data.index_before = 1:4; 
        data.index_after = 5:10;
        data.photobleach_time = 0; 
        data.dt = 10.0/60; %min
        data.subtract_background = 1;
        data.crop_image = 1;
        data.median_filter = 1;
        data.magnification = 60;
        data.time_format = 2; 
        data.brightness_factor = 1.2;
    case 'yf2_1023' % YF2 cell
        data.path = strcat(root, 'rongxue_1023/yf_1023/2/');
        data.output_path = strcat(data.path, 'output/');
        data.first_file = '2_w1Imaging RFP_s1_t1.TIF';
        data.index_pattern = {'t1', 't%d'};
        data.image_index = (1:10);
        data.index_before = 1:4; 
        data.index_after = 5:10;
        data.photobleach_time = 0; 
        data.dt = 10.0/60; %min
        data.subtract_background = 1;
        data.crop_image = 1;
        data.median_filter = 1;
        data.magnification = 60;
        data.time_format = 2; 
        data.brightness_factor = 1.4;
    case 'yf1_1023' % YF1 cell
        data.path = strcat(root, 'rongxue_1023/yf_1023/1/');
        data.output_path = strcat(data.path, 'output/');
        data.first_file = '1_w1Imaging RFP_s1_t1.TIF';
        data.index_pattern = {'t1', 't%d'};
        data.image_index = (1:10);
        data.index_before = 1:4; 
        data.index_after = 5:10;
        data.photobleach_time = 0; 
        data.dt = 10.0/60; %min
        data.subtract_background = 1;
        data.crop_image = 1;
        data.median_filter = 1;
        data.magnification = 60;
        data.time_format = 2; 
        data.brightness_factor = 2.0;
    case 'yf1_1018' % YF1 cell
        data.path = strcat(root, 'rongxue_1023/yf_1018/1/');
        data.output_path = strcat(data.path, 'output/');
        data.first_file = '1_w1Imaging RFP_s1_t1.TIF';
        data.index_pattern = {'t1', 't%d'};
        data.image_index = (1:10);
        data.index_before = 1:4; 
        data.index_after = 5:10;
        data.photobleach_time = 0; 
        data.dt = 10.0/60; %min
        data.subtract_background = 1;
        data.crop_image = 1;
        data.median_filter = 1;
        data.magnification = 60;
        data.time_format = 2; 
    case 'yf2_1018' % YF2 cell
        data.path = strcat(root, 'rongxue_1023/yf_1018/2/');
        data.output_path = strcat(data.path, 'output/');
        data.first_file = '2_w1Imaging RFP_s1_t1.TIF';
        data.index_pattern = {'t1', 't%d'};
        data.image_index = (1:10);
        data.index_before = 1:4; 
        data.index_after = 5:10;
        data.photobleach_time = 0; 
        data.dt = 10.0/60; %min
        data.subtract_background = 1;
        data.crop_image = 1;
        data.median_filter = 1;
        data.magnification = 60;
        data.time_format = 2; 
        data.brightness_factor = 1.5;
    case 'yf3_1018' 
        data.path = strcat(root, 'rongxue_1023/yf_1018/3/');
        data.output_path = strcat(data.path, 'output/');
        data.first_file = '3_w1Imaging RFP_s1_t1.TIF';
        data.index_pattern = {'t1', 't%d'};
        data.image_index = (1:10);
        data.index_before = 1:4; 
        data.index_after = 5:10;
        data.photobleach_time = 0;
        data.dt = 10.0/60; %min
        data.subtract_background = 1;
        data.crop_image = 1;
        data.median_filter = 1;
        data.magnification = 60;
        data.time_format = 2; 
        data.brightness_factor = 1.5;
    case 'yf4_1018' 
        data.path = strcat(root, 'rongxue_1023/yf_1018/4/');
        data.output_path = strcat(data.path, 'output/');
        data.first_file = '4_w1Imaging RFP_s1_t1.TIF';
        data.index_pattern = {'t1', 't%d'};
        data.image_index = (1:10);
        data.index_before = 1:4; 
        data.index_after = 5:10;
        data.photobleach_time = 0;
        data.dt = 10.0/60; %min
        data.subtract_background = 1;
        data.crop_image = 1;
        data.median_filter = 1;
        data.magnification = 60;
        data.time_format = 2; 
        data.brightness_factor = 1.5;
    case 'yf5_1018' 
        data.path = strcat(root, 'rongxue_1023/yf_1018/5/');
        data.output_path = strcat(data.path, 'output/');
        data.first_file = '5_w1Imaging RFP_s1_t1.TIF';
        data.index_pattern = {'t1', 't%d'};
        data.image_index = (1:10);
        data.index_before = 1:4; 
        data.index_after = 5:10;
        data.photobleach_time = 0;
        data.dt = 10.0/60; %min
        data.subtract_background = 1;
        data.crop_image = 1;
        data.median_filter = 1;
        data.magnification = 60;
        data.time_format = 2; 
        data.brightness_factor = 1.5;
    case 'yf6_1018' 
        data.path = strcat(root, 'rongxue_1023/yf_1018/6/');
        data.output_path = strcat(data.path, 'output/');
        data.first_file = '6_w1Imaging RFP_s1_t1.TIF';
        data.index_pattern = {'t1', 't%d'};
        data.image_index = (1:10);
        data.index_before = 1:4; 
        data.index_after = 5:10;
        data.photobleach_time = 0;
        data.dt = 10.0/60; %min
        data.subtract_background = 1;
        data.crop_image = 1;
        data.median_filter = 1;
        data.magnification = 60;
        data.time_format = 2; 
        data.brightness_factor = 1.5;
    case 'wt3_1018' 
        data.path = strcat(root, 'rongxue_1023/wt_1018/3/');
        data.output_path = strcat(data.path, 'output/');
        data.first_file = '3_w1Imaging RFP_s1_t1.TIF';
        data.index_pattern = {'t1', 't%d'};
        data.image_index = (1:10);
        data.index_before = 1:3; 
        data.index_after = 4:9;
        data.photobleach_time = 0;
        data.dt = 10.0/60; %min
        data.subtract_background = 1;
        data.crop_image = 1;
        data.median_filter = 1;
        data.magnification = 60;
        data.time_format = 2; 
        data.brightness_factor = 1.25;
    case 'wt4_1018' 
        data.path = strcat(root, 'rongxue_1023/wt_1018/4/');
        data.output_path = strcat(data.path, 'output/');
        data.first_file = '4_w1Imaging RFP_s1_t1.TIF';
        data.index_pattern = {'t1', 't%d'};
        data.image_index = (1:10);
        data.index_before = 1:3; 
        data.index_after = 4:9;
        data.photobleach_time = 0;
        data.dt = 10.0/60; %min
        data.subtract_background = 1;
        data.crop_image = 1;
        data.median_filter = 1;
        data.magnification = 60;
        data.time_format = 2; 
        data.brightness_factor = 1.2;
    case 'wt5_1018' 
        data.path = strcat(root, 'rongxue_1023/wt_1018/5/');
        data.output_path = strcat(data.path, 'output/');
        data.first_file = '5_w1Imaging RFP_s1_t1.TIF';
        data.index_pattern = {'t1', 't%d'};
        data.image_index = (1:10);
        data.index_before = 1:3; 
        data.index_after = 4:9;
        data.photobleach_time = 0;
        data.dt = 10.0/60; %min
        data.subtract_background = 1;
        data.crop_image = 1;
        data.median_filter = 1;
        data.magnification = 60;
        data.time_format = 2; 
        data.brightness_factor = 1.2;
    case 'wt6_1018' 
        data.path = strcat(root, 'rongxue_1023/wt_1018/6/');
        data.output_path = strcat(data.path, 'output/');
        data.first_file = '6_w1Imaging RFP_s1_t1.TIF';
        data.index_pattern = {'t1', 't%d'};
        data.image_index = (1:10);
        data.index_before = 1:3; 
        data.index_after = 4:9;
        data.photobleach_time = 0;
        data.dt = 10.0/60; %min
        data.subtract_background = 1;
        data.crop_image = 1;
        data.median_filter = 1;
        data.magnification = 60;
        data.time_format = 2; 
        data.brightness_factor = 1.2;
    case 'wt7_1018' 
        data.path = strcat(root, 'rongxue_1023/wt_1018/7/');
        data.output_path = strcat(data.path, 'output/');
        data.first_file = '7_w1Imaging RFP_s1_t1.TIF';
        data.index_pattern = {'t1', 't%d'};
        data.image_index = (1:10);
        data.index_before = 1:4; 
        data.index_after = 5:10;
        data.photobleach_time = 0;
        data.dt = 10.0/60; %min
        data.subtract_background = 1;
        data.crop_image = 1;
        data.median_filter = 1;
        data.magnification = 60;
        data.time_format = 2; 
        data.brightness_factor = 1.2;
    case 'wt8_1018' 
        data.path = strcat(root, 'rongxue_1023/wt_1018/8/');
        data.output_path = strcat(data.path, 'output/');
        data.first_file = '8_w1Imaging RFP_s1_t1.TIF';
        data.index_pattern = {'t1', 't%d'};
        data.image_index = (1:10);
        data.index_before = 1:4; 
        data.index_after = 5:10;
        data.photobleach_time = 0;
        data.dt = 10.0/60; %min
        data.subtract_background = 1;
        data.crop_image = 1;
        data.median_filter = 1;
        data.magnification = 60;
        data.time_format = 2; 
        data.brightness_factor = 1.2;
    case 'wt9_1018' % not used, cell moved too much, wt10 had detection problem
        data.path = strcat(root, 'rongxue_1023/wt_1018/9/');
        data.output_path = strcat(data.path, 'output/');
        data.first_file = '9_w1Imaging RFP_s1_t1.TIF';
        data.index_pattern = {'t1', 't%d'};
        data.image_index = (1:10);
        data.index_before = 1:4; 
        data.index_after = 5:10;
        data.photobleach_time = 0;
        data.dt = 10.0/60; %min
        data.subtract_background = 1;
        data.crop_image = 1;
        data.median_filter = 1;
        data.magnification = 60;
        data.time_format = 2; 
        data.brightness_factor = 1.2;
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





