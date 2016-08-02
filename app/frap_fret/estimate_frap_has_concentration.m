% function estimate_frap_has_concentration(cell_name,starting_step)

% Copyright: Shaoying Lu and Yingxiao Wang 2012-2016
function estimate_frap_has_concentration(cell_name,starting_step)
    % Initialization
    data = set_parameter(cell_name);
    path = data.path; 
    shift = data.shift;
    dt = data.dt;
    mag = data.magnification;
    display(path); display(shift); display(dt); display(mag);

    num_steps = 2;
    num_images_per_layer = 1;
    num_refines = 2;
    show_figure_result =1;

    % load boundary and med_con
    file_name = sprintf('%sresults_%d.mat',path,starting_step);
    res = load(file_name);
    med_con = res.med_con;
    con = res.con;

    % mark a more restricted boundary
    boundary_file = strcat(path, 'smaller_boundary.data');
    % in the new version, need to change mark_boundary to get_polygon
    smaller_boundary = mark_boundary(con(:,:,:,1), boundary_file);
    clear con;

    % create mesh and assemble the matrices
    mesh = create_mesh(smaller_boundary);
    new_mesh = refine_mesh(mesh, 'num_refines', num_refines, 'method',1);
    p_image = new_mesh.node; edge = new_mesh.edge;
    tri = new_mesh.tri;
    clear mesh new_mesh;
    p = scale_by_magnification(p_image,mag);
    [K,M] = assemble_matrix(p,tri);
    is_boundary = mark_mesh_boundary(p_image, edge);

    % convert concentration to solution vector
    num_nodes = size(p_image,2);
    display(num_nodes);
    u = zeros(num_nodes,num_steps);
    est_u = zeros(num_nodes,num_steps);
    r = zeros(num_nodes, num_steps);
    for i = 1:num_steps,
        u(:,i) = concentration_to_vector(med_con(:,:,:,i),p_image);
    end;

    % estimate diffusion constant
    display(sprintf('%20s%20s%20s%20s',...
        'diffusion coef', 'residual_1', 'residual_2','R value'));
    dt = dt*num_images_per_layer;
    diff_const = zeros(num_steps-1,1);
    for i= 1:num_steps-1,
        %Mdu = M*(u(:,i+1)-u(:,i));
        %dtKu = -dt*K*0.5*(u(:,i+1)+u(:,i));
        [diff_const(i) , est_u(:,i+1) , R, r(:,i+1)] = ...
            get_diffusion_constant(u(:,i),u(:,i+1),M,K,dt,is_boundary);
                
        relative_residual_norm = norm(u(:,i+1)-est_u(:,i+1),2)/...
            norm(u(:,i+1),2);
        relative_residual_norm_2 = norm(u(:,i+1)-est_u(:,i+1),2)/...
            norm(u(:,i+1)-u(:,i));
        display(sprintf('%20f%20f%20f%20f', diff_const(i), ...
            relative_residual_norm, relative_residual_norm_2, R));
        %du = u(:,i+1)-u(:,i);
    end;

    if show_figure_result,
        cbound = [0.45, 1.0];
        s = max(p_image)-min(p_image);
        width = s(1)*30; 
        height = s(2)*20;
        for i = 1:2,
            h = figure('position', [0, 150, width, height]); pdesurf(p_image,tri,u(:,i));
            view(90,90); colormap jet; caxis(cbound); colorbar;
            set(gca, 'FontSize', 16);
            title('Concentration map'); 
            saveas(h,sprintf('%su_%d.fig',path,i+starting_step));
            saveas(h,sprintf('%su_%d.png',path,i+starting_step));
        end;
        for i =2:2,
            h = figure('position', [0, 150, width, height]); 
            pdesurf(p_image,tri,est_u(:,i));
            view(90,90); colormap jet; caxis(cbound); colorbar;
            set(gca, 'FontSize', 16);
            title('Estimated concentration map'); 
            saveas(h,sprintf('%sest_u_%d.fig',path,i+starting_step));
            saveas(h,sprintf('%sest_u_%d.png',path,i+starting_step));
        end;

        error = zeros(size(u));
        for i = 1:1,
           %  error(:,i+1) = (u(:,i+1)-est_u(:,i+1))/abs(diff_const(i))/dt;
            error(:,i+1) = (u(:,i+1)-est_u(:,i+1))/norm(u(:,i+1)-u(:,i));
        end;
    %     for i = 2:3,
    %         [R,P]=corrcoef(error(:,i), error(:,i+1))
    %     end;

        cbound = [-0.005, 0.02];
        for i = 1:1,
            h = figure('position', [0,150, width, height]); 
            pdesurf(p_image,tri, abs(error(:,i+1)));
            view(90,90); colormap jet; caxis(cbound);colorbar;
            set(gca, 'FontSize', 16);
            title('abs(u - est u) '); 
            saveas(h,sprintf('%sabs_diff_u_%d_est_u_%d.fig',path,i+1,i+1));
            saveas(h,sprintf('%sabs_diff_u_%d_est_u_%d.png',path,i+1,i+1));
         end;

        for i = 1:1,
            [h1, h2] = plot_Mdu_dtKu_residual(diff_const(i), u(:,i), u(:,i+1), M, K, dt);
            saveas(h1, sprintf('%sMdu_dtKu_%d.fig', path, i+starting_step));
            saveas(h1, sprintf('%sMdu_dtKu_%d.png', path, i+starting_step));
            saveas(h2, sprintf('%sresidual_%d.fig', path, i+starting_step));
            saveas(h2, sprintf('%sresidual_%d.png', path, i+starting_step));
        end;

    end;
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
    plot(Mdu, r0/dt,'+'); hold on; box off;
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

