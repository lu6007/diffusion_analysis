% function fh = diffusion_function()
%     fh.draw_surface_with_mesh = @draw_surface_with_mesh;
%     fh.adjust_time_step = @adjust_time_step; 
%
% function draw_surface_with_mesh(surf, mesh)
%
% function [G, time_stepping] = ...
%    adjust_time_stepping(diff_const,dt,u_0,M,K)
% adjust time stepping 
% until the solution converges
% using the backward Eular method

% Copyright: Shaoying Lu and Yingxiao Wang, Email: shaoying.lu@gmail.com
function fh = diffusion_function()
    fh.draw_surface_with_mesh = @draw_surface_with_mesh;
    fh.adjust_time_step = @adjust_time_step; 
return

function draw_surface_with_mesh(surf, mesh)
    set(gcf, 'color', 'w');
    %
    u = surf.u; 
    pdesurf(surf.p,surf.tri,u'); hold on;
    %
    num_node = length(mesh.node);
    pdemesh(mesh.node, mesh.edge, mesh.tri, (max(max(u))+1)*ones(1,num_node));
    view(2); axis ij; colormap jet; 
    set(gca, 'FontSize', 12, 'FontWeight','bold',...
             'Box', 'off', 'LineWidth', 1.5);
    colorbar; 
return

function [G, time_step] =  adjust_time_step(dt,u_0,M,K)
    solution_error = 5.0e-2*dt;
    G = (M+dt*K);
    while 1
        u_1 = G\(M*u_0);
        dtt = dt/2.0;
        H = (M+dtt*K);
        v_1 = H\(M*(H\(M*u_0)));
        err = norm((v_1-u_1)/size(v_1,1));
        if err>solution_error
            dt = dtt;
            clear G;
            G = H; clear H;
            continue;
        else
            time_step = dt;
            break;
        end
    end
return