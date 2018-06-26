% function computer_simulation(data, varargin)
% test the diffusion_simulation function
%     parameter_name = {'load_file', 'save_file','dt', 'estimate_diffusion', 'mark_subregion',...
%         'simulation_time'}; % simulation_time in seconds
%     default_value = { false, true, 0.25, false, false, 10};
%
% Sub-functions
% function mark_mesh_subregion(mesh)
%
% Example:
% >> cell_name = 'phtobleach_cell';
% >> data = sample_diffusion_init_data(cell_name);
% >> computer_simulation(data);
%
% The variable cell_name can take the values of 'photobleach_cell',
% 'photobleach_cell_2', 'square_5_circle', 'layered_diffusion',
% 'tensor_diffusion', etc

% Copyright: Shaoying Lu and Yingxiao Wang 2012-2016
function data = computer_simulation(data, varargin)
    parameter_name = {'load_file', 'save_file','dt', 'estimate_diffusion', 'mark_subregion',...
        'simulation_time'}; % simulation_time in seconds
    default_value = { false, true, data.dt, false, false, 10};
    [load_file, save_file, dt, estimate_diffusion, mark_subregion, simulation_time] = ...
        parse_parameter(parameter_name, default_value, varargin);
    data.dt = dt;

    % Initialize data.
    cell_name =data.cell_name;
    pa = data.path;
    diff_const = data.diff_const;
    num_bits = data.num_bits;
    output_dir =strcat(pa, 'output/');
    if ~exist(output_dir, 'dir')
        mkdir(output_dir);
    end

    % Initialize the first image after photobleach
    image_0 = imread(strcat(pa,data.first_file));
    temp = imcrop(image_0, data.rectangle); clear image_0;
    if num_bits ==8 % convert from uint8 to uint16
        image_0 = uint16(double(temp))*2^(16-num_bits);
        %num_bits = 16;
    else % num_bits == 16,
        image_0 = temp;
    end; clear temp;
    data.image_0 = image_0;

    % Calculate the diffusion map image data.diff_map
    % Obtain the boundary of the cell and the photobleach region (same in
    % all the cell diffusion cases)
    data.boundary = simulation_get_boundary(cell_name, data, 'load_file', load_file);
    data.diffusion_map_file = strcat(pa, 'diffusion_map.tiff'); % need data.boundary{1};
    data = init_diffusion_map(cell_name, data, 'load_file', load_file, 'save_file', save_file);
    % obtaine the boundary of subdomains with different diffusion
    % coefficients
    boundary = simulation_get_boundary(cell_name, data, 'option',2,'load_file', true);

    % Create the initial mesh and refine the mesh
    % p_image: node points with x and y values
    % edge and tri defines the mesh.
    % dl is probably the boundary condition
    rr = data.rectangle;

    if ~iscell(boundary)
        % create mesh with 1 boundary
        mesh = create_mesh(boundary, 'method', 1);
    else % create mesh with multi boundaries
        mesh = create_mesh(boundary, 'method', 2);
        temp = boundary{1}; clear boundary;
        boundary = temp; clear temp;
    end

    c_max = 2^16;
    color_bound = [0 c_max-1];
    axis_bound = [1 rr(3) 1 rr(4) 1 c_max];
    graph_pos = [200 200 rr(3) rr(4)];

    % draw image and plot initial mesh
    figure(1);
    set(gcf, 'color', 'w'); hold on;
    image(image_0);
    view(2); axis ij;
    colormap(gray);
    set(gca, 'FontSize', 12, 'FontWeight','bold',...
             'Box', 'off', 'LineWidth', 1.5);
    title('Image Right after Photobleach');


    if isfield(data, 'diff_map')
        figure(2);
        set(gcf, 'color', 'w');
        imagesc(data.diff_map); hold on;
        pdemesh(mesh.node,mesh.edge,mesh.tri);
        view(2); axis ij;
        set(gca, 'FontSize', 12, 'FontWeight','bold',...
                 'Box', 'off', 'LineWidth', 1.5);
        title('Initial Triangular Mesh');
    end

    % refine mesh
    refine_improve = 2;
    new_mesh = refine_mesh(mesh, 'method',refine_improve, 'num_refines',1);
    % p_image contains the (x,y) positions of the nodes in the image
    % matrix, not the row-index and col-index in the matrix.
    p_image = new_mesh.node; edge = new_mesh.edge; tri = new_mesh.tri;
    dl = new_mesh.dl;
    clear new_mesh;
    mag = data.mag;
    p = scale_by_magnification(p_image,mag);

    % plot mesh
    h = figure(3);
    set(gcf, 'color', 'w'); hold on;
    image(image_0);
    pdemesh(p_image,edge,tri);
    view(2); axis ij; colormap(gray);
    set(gca, 'FontSize', 12, 'FontWeight','bold',...
             'Box', 'off', 'LineWidth', 1.5);
    title('Refined Mesh');
    if mark_subregion
        mark_mesh_subregion(mesh,'handle', h);
    end

    % Compute initial solution
    lw = 0.5; % LineWidth
    p_ij = [p_image(2,:); p_image(1,:)]; % switch from (x,y) to (i,j)
    v_0 = concentration_to_vector(image_0(:,:,1),p_ij);

    figure(4);
    set(gcf, 'color', 'w'); hold on;
    x1 = boundary(:,1); y1 = boundary(:,2);
    plot3(x1, y1, c_max*ones(size(x1)),'r', 'LineWidth', lw);
    pdesurf(p_image,tri,v_0);
    set(gca, 'FontSize', 12, 'FontWeight', 'Bold', 'LineWidth',1.5, ...
             'YDir', 'reverse');
    view(2); shading interp; caxis(color_bound); colormap(gray);
    colorbar;
    title('Initial Concentration before Smoothing');

    % Calculate the diffusion vector
    data.tri_centroid = (p_image(:,tri(1,:)) + p_image(:,tri(2,:)) + p_image(:,tri(3,:))) / 3.0;
    data = simulation_get_diffusion_vector(cell_name, data, tri(4,:));
    diff_coef = data.diff_coef;
    diff_tag  = data.diff_tag;

    % Visualize the diffusion vector diff_coef
    if isfield(data, 'diff_map')
        figure; set(gcf, 'color', 'w');
        num_nodes = length(mesh.node);
        pdesurf(p_image,tri,(diff_coef)'); hold on;
        pdemesh(mesh.node,mesh.edge,mesh.tri, (max(max(diff_coef))+1)*ones(1,num_nodes));
        view(2); axis ij;
        set(gca, 'FontSize', 12, 'FontWeight','bold',...
                 'Box', 'off', 'LineWidth', 1.5);
        title('Diffusion Coefficients Overlaid with Mesh');
        clear mesh;
    end

    if isfield(data, 'diff_map')
        figure; set(gcf, 'color', 'w');
        num_nodes = length(p_image);
        pdesurf(p_image,tri,(diff_coef)'); hold on;
        pdemesh(p_image,edge,tri, (max(max(diff_coef))+1)*ones(1,num_nodes));
        view(2); axis ij;
        set(gca, 'FontSize', 12, 'FontWeight','bold',...
                 'Box', 'off', 'LineWidth', 1.5);
        title('Diffusion Coefficients Overlaid with Refined Mesh');
        clear mesh;
    end

    % assemble linear system
    % let the initial solution diffuse for t0 sec
    % so that it smoothes out at boundary;
    % [K M] = assemble_matrix(p, tri);
    [K_diff_coef, M] = assemble_matrix(p, tri, 'diff_coef', diff_coef');
    t0 = 2.9/mean(diff_const);
    backward_eular= 2;
    u_0 = simulate_diffusion(v_0, t0, K_diff_coef,M,'method', backward_eular);

    % adjust time stepping
    % until the solution converges
    % using the backward Eular method
    % For better format, use str_1 = sprintf(...); display(str_1);
    display(dt);

    % simulate diffusion for <= 30 seconds
    % until the fluorescence recovery stops.
    % Compute fluorescence recovery curve.
    num_steps = floor(simulation_time/dt); % default = 40;
    show_image = 40;
    u = zeros(length(u_0),num_steps+1);
    for i=1:num_steps+1
        if i ==1
            u(:,i) = u_0;
        else
            u(:,i) = simulate_diffusion(u(:,i-1), dt, ...
                K_diff_coef,M,'method', backward_eular);
        end

        % Save images
        file_name = strcat(output_dir, num2str(i), '.tiff');
        if ((i<=6) || i ==floor(i/show_image)*show_image)
            figure;
            i
            set(gcf, 'position', graph_pos, 'color', 'w');
            pdesurf(p_image,tri,u(:,i)); hold on;

            if ~exist(file_name, 'file')
                cmap = gray;
            else
                cmap = jet;
            end

            view(2); shading interp; colormap(cmap);
            axis(axis_bound); axis off; caxis([0 2^16]);
            set(gca, 'position', [0 0 1 1],'YDir','reverse');
            pause(1);

            if ~exist(file_name, 'file')

                % The print statement is needed for the getframe function
                print -dbitmap
                ff = getframe(gca);

                % Convert the output of getframe from uint8 to uint16
                im_i = 2^8*uint16(double(ff.cdata(:,:,1)));
                imwrite(im_i,file_name,'tiff','compression','none');
                clear ff im_i;

            end
        end
        clear file_name;
    end

    output_file = strcat(output_dir,'simulation_result.mat');
    output_file2 = strcat(output_dir,'mass_stiffness_matrix.mat');

    if save_file
        save(output_file, 'u','v_0','M','p','p_image','tri','edge');
        save(output_file2, 'M','K_diff_coef');
    end

    figure;
    set(gcf, 'position', graph_pos, 'color', 'w');
    pdesurf(p_image,tri,u(:,6)-u(:,5)); hold on; %12/3/14
    view(2); shading interp;
    colormap(jet);
    axis(axis_bound); axis off; caxis([-500 500]);
    set(gca, 'position', [0 0 1 1],'YDir','reverse');
    clear time;

    if estimate_diffusion
        % % estimate the diffusion constant
        % % using the last quarter of solution profile.
        % % The Crank-Nicolson Scheme is used for accuracy
        % i_0 = floor(0.75*num_steps);
        % drive = -dt*K*0.5*(u(:,i_0)+u(:,i_0+1));
        % Mdu = M*(u(:,i_0+1)-u(:,i_0));
        % for i= i_0+1:num_steps-1,
        %     drive = [drive; -dt*K*0.5*(u(:,i)+u(:,i+1))];
        %     Mdu = [Mdu; M*(u(:,i+1)-u(:,i))];
        % end;
        % est_diff_const = drive\Mdu

        % Estimate the diffusion constant using backward Eular,
        % then refine mesh and test the estimation of diff_const
        [K, M] = assemble_matrix(p, tri);
        [est_diff_const, ~, ~] = get_diffusion_constant(u(:,1),u(:,2), M, K, dt);
        display(est_diff_const);

        [p1_image, ~, t1, u1] = refinemesh(dl,p_image,edge,tri,u(:,1:2));
        p1 = scale_by_magnification(p1_image,mag);

        [K1, M1] = assemble_matrix(p1, t1);
        [est_diff_const, ~, ~] = get_diffusion_constant(u1(:,1), u1(:,2),...
                                                        M1, K1, dt);
        display(est_diff_const);
    end

    % output for pltmg
    if save_file

        vxvyu2u3 = [p' u(:,2:3)];
        num_nodes = size(p,2);   % 11153
        num_tris = size(tri,2);  % 21888
        itnodes = tri';          % num_tris x 4
        num_edges = size(edge,2);
        ibndary = edge([1 2 7],:)'; % num_edges x 3
        % The 4th column of itnodes marks the region of photobleach
        % it is not important in the examples.
        % In the layered problem the 4th columns is the labeling of the
        % subregions.
        if isfield(data, 'diff_map') && ~isempty(diff_tag)
            itnodes(:,4) = diff_tag;
        end

        save(strcat(output_dir,cell_name, '.data'), 'num_nodes','vxvyu2u3',...
            'num_tris', 'itnodes', 'num_edges', 'ibndary','-ascii');
        clear vxvyu2u3 ibndary itnodes;

    end

return;

% function mark_mesh_subregion(mesh)
% Mark the triangles with the subregion it belongs on the current figure
% Mark the edges with the subregion it belongs.

% Copyright: Shaoying Lu and Yingxiao Wang 2016
function mark_mesh_subregion(mesh, varargin)

    parameter_name = {'handle'};
    default_value = {gcf};
    handle = parse_parameter(parameter_name, default_value, varargin);

    figure(handle); hold on;
    ee = mesh.edge([1, 2, 7],:)';
    pp = mesh.node';
    mm = 0.5*(pp(ee(:,1),:)+pp(ee(:,2),:));
    for i = 0:10
        ii = (ee(:, 3)==i);
        text(mm(ii,1), mm(ii,2), num2str(i),'color','k');
        clear ii;
    end

    mm1 = 0.75*pp(ee(:,1),:)+0.25*pp(ee(:, 2),:);
    mm2 = 0.25*pp(ee(:,1),:)+0.75*pp(ee(:, 2),:);
    text(mm1(:,1), mm1(:,2), 'o');
    text(mm2(:,1), mm2(:,2), '+');
    clear mm mm1 mm2;

    % mark triangles
    tt = mesh.tri';
    for i = 1:8
        ii = (tt(:, 4) == i);
        mm = 0.3333*(pp(tt(ii,1),:)+pp(tt(ii,2),:)+pp(tt(ii,3),:));
        text(mm(:,1), mm(:,2), num2str(i), 'color', 'r');
        clear ii mm;
    end

return
