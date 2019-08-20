function fh = diffusion_function()
    fh.draw_surface_with_mesh = @draw_surface_with_mesh;
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