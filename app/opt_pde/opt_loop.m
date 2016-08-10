% outer loop for optimization 
% uses solution from previous steps as initial guess for later iterations
% Created by Yiwen Shi, 07/2016, UCSD

% tensor_diffusion: 99 steps, [4.9987, 93.1667, 5.1059]
% spot_diffusion: 8 steps, [28.9978, 5, 5, 5, 5]
% layered_diffusion: 5 steps, [5, 19.9974]


function d = opt_loop(data, tol, cell_name)

    load(data);

    if strcmp(cell_name, 'tensor_cross_2')
        tri(4, tri(4, :) == 3) = 2;
        tri(4, tri(4, :) == 4) = 3;
        tri(4, tri(4, :) == 5) = 3;
        tri(4, tri(4, :) == 6) = 4;
        num_para = size(unique(tri(4, :)), 2);
    end
    
    d0 = ones(num_para, 1);
    d = ones(num_para, 1) * 10;

    step = 1;
    while norm(d0 - d) > tol
        d0 = d;
        d'
        [u,d, s_hist, err_hist, linerr_hist] = opt3(u0, u1, u2, d0, dt, num_nodes, tri, p, cell_name);
        fprintf('Outer Step %d \n', step);
        step = step + 1;
    end