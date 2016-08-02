% function compute_recovery_curve(data)
%
% Example:
% >> cell_name = 'photobleach_cell';
% >> data = sample_diffusion_init_data(cell_name);
% >> compute_recovery_curve(data);
%

% Copyright: Shaoying Lu and Yingxiao Wang 2013-2016
function compute_recovery_curve(data)
    pa = data.path;
    dt = data.dt;
    file = strcat(pa, 'output/simulation_result.mat');
    res =load(file);
    u = res.u;
    M = res.M;
    v_0 = res.v_0;
    num_steps = size(u,2);
    r = zeros(num_steps,1);
    t = (1:num_steps)'*dt;
    dd = mark_photobleach_region(v_0,'bound', 10000,'file',strcat(pa,'pb_region.mat'));
    r_0 = dd'*M*(ones(size(u(:,1))).*dd);
    for i = 1:num_steps,
        r(i) = dd'*M*(u(:,i).*dd)/r_0;
    end;
    time = [-2; -1; 0; t];
    recovery_curve = [65535; 65535; 65535; r];
    figure;
    plot(time, recovery_curve, 'LineWidth',1.5);
    axis([-4, 10, 10000,70000]);
    set(gca, 'XTick', [-4, 0, 4, 8], 'YTick', [1e4, 2e4,3e4,4e4,5e4,6e4,7e4]);
    ylabel('FI (AU)'); xlabel('Time (Sec)');
    title('Photobleach Recovery Curve');
return;
    
