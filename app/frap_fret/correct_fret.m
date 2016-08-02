% Correct for biosensor diffusion and enhance the signal.
% function correct_fret(data);
% 
% Example:
% >> cell_name = 'mem17';
% >> data = sample_diffusion_init_data(cell_name);
% >> correct_fret(data);

% Copyright: Shaoying Lu and Yingxiao Wang 2012-2016
function correct_fret(data)
    cell_name = data.cell_name; 
    path = data.path;
    first_cfp_file = data.first_cfp_file;
    magnification = data.magification;
    diff_coef = data.diff_coef;
    channel_pattern = data.channel_pattern;
    image_index = data.image_index;
    %

    % Mark rectangle and boundary
    fret_file = regexprep(first_cfp_file, channel_pattern{1}, channel_pattern{2});
    temp = imread(strcat(path, data.first_cfp_file));
    cfp_im = preprocess(temp, data); clear temp;
    temp = imread(strcat(path, fret_file));
    fret_im = preprocess(temp, data); clear temp;
    ratio = compute_ratio(cfp_im, fret_im);
    im_1 = get_imd_image(ratio, fret_im, 'ratio_bound', data.ratio_bound,...
        'intensity_bound', data.intensity_bound);
    clear index  ratio yfp;
    %
    rectangle_file = strcat(path, 'image_rectangle.data');
    image_rectangle = get_rectangle(im_1,rectangle_file,'format', '-ascii');
    im_1 = imcrop(im_1, image_rectangle);
    boundary_file = strcat(path, 'boundary.mat');
    [~, temp] = get_polygon(im_1, boundary_file, 'format', '-mat');
    boundary = temp{1};
    %boundary = mark_boundary(im_1, boundary_file);
    % 
    index = sprintf(data.index_pattern{2}, max(image_index));
    fret_file = regexprep(fret_file, data.index_pattern{1}, index);
    temp = imread(strcat(path, fret_file));
    last_fret = imcrop(temp, image_rectangle);
    figure; imshow(last_fret); caxis auto; hold on;
    plot(boundary(:,1), boundary(:,2),'--','LineWidth',2);
    title('Check the cropping rectangle on the last image');
    clear index fret_file temp last_fret;

    % create mesh
    num_refines = 2;
    mesh = create_mesh(boundary);
    new_mesh = refine_mesh(mesh, 'num_refines', num_refines);
    p_image = new_mesh.node;  tri = new_mesh.tri;
    figure; pdemesh(mesh.node, mesh.edge, mesh.tri);
    title('Initial mesh');
    figure; pdemesh(new_mesh.node, new_mesh.edge, new_mesh.tri);
    title('Refined mesh');
    clear mesh new_mesh;

    % assemble matrices 
    p = scale_by_magnification(p_image,magnification);
    [K_diff_coef,M] = assemble_matrix(p,tri,'diff_coef', diff_coef);
    num_steps =length(data.image_index);
    %num_steps = 4;
    num_nodes = size(p_image,2);
    u = zeros(num_nodes, num_steps);
    est_u = zeros(size(u));
    act = zeros(size(u));
    u_act = zeros(size(u));
    for i = 1:num_steps-1,
        fprintf('i = %d\n', i); 
        index = sprintf(data.index_pattern{2}, image_index(i));
        cfp_file = regexprep(first_cfp_file, data.index_pattern{1}, index);
        fret_file = regexprep(cfp_file, channel_pattern{1}, channel_pattern{2});
        temp = imread(strcat(path, cfp_file));
        cfp_im = preprocess( temp, data); clear temp;
        temp = imread(strcat(path, fret_file));
        fret_im = preprocess( temp, data); clear temp;
        temp = compute_ratio(cfp_im, fret_im);
        ratio = imcrop(temp, image_rectangle); clear temp;
        u(:,i) = concentration_to_vector(ratio, p_image, 'method', 2);
        if strcmp(cell_name, 'egf_pp1_lyn1'),
            dt = special_treat_egf_pp1_lyn1(image_index(i));
        end;
        est_u(:,i+1) = simulate_diffusion(u(:,i), dt, K_diff_coef,M,'method', 2);
        clear index cfp_file yfp_file ratio;
        clear cfp_im fret_im;
    end;
    index = sprintf(data.index_pattern{2}, image_index(num_steps));
    cfp_file = regexprep(first_cfp_file, data.index_pattern{1}, index);
    fret_file = regexprep(cfp_file, channel_pattern{1}, channel_pattern{2});
    temp = imread(strcat(path, cfp_file));
    cfp_im = preprocess( temp, data); clear temp;
    temp = imread(strcat(path, fret_file));
    fret_im = preprocess( temp, data); clear temp;
    temp = compute_ratio(cfp_im, fret_im);
    ratio = imcrop(temp, image_rectangle); clear temp;
    u(:,num_steps) = concentration_to_vector(ratio, p_image, 'method', 2);

    u_act(:,1) = u(:,1);
    for i = 2: num_steps,
        act(:,i) = u(:,i)-est_u(:,i);
        u_act(:,i) = u_act(:,i-1)+act(:,i);
    %     i
    %     sum_u = sum(M*u(:,i))
    %     sum_u_act = sum(M*u_act(:,i))
    end;
    sum_Mu = sum(M*u(:,num_steps));
    display(sum_Mu);
    %
    output_dir = strcat(path, 'output/');
    if ~isdir(output_dir),
        mkdir(output_dir);
    end;
    save(strcat(output_dir,'result.mat'), 'u', 'act', 'u_act',...
        'tri', 'p_image', 'p', 'M', 'K_diff_coef', 'image_index', 'dt');
return;


function dt = special_treat_egf_pp1_lyn1(index)
% get dt in seconds
    dt = 30;
    if index == 12,
        dt = 16.19;
    elseif index == 16,
        dt = 16.89;
    elseif index ==17,
        dt = 42.19;
    elseif index >=18,
            dt = 60;
    end;
return;

