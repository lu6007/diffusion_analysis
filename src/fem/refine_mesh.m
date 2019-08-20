% function mesh = refine_mesh(mesh, num_refines)
% parameter_name = {'method', 'num_refines'};
% default_value = {1, 1};
% refine the triangular mesh based on the boundary and
% num_refinement provided.
% method 1 -- refine and improve mesh quality, output 1 mesh.
% method 2 -- refine only, do not jigglemesh.

% Copyright: Shaoying Lu and Yingxiao Wang 2012
function new_mesh = refine_mesh(mesh, varargin)
para_name = {'method', 'num_refines'};
default_v = {1, 1};
[method, num_refines] = parse_parameter(para_name, default_v, varargin);
switch method
    case 1 
        % method 1 -- refine and improve mesh quality, output 1 mesh.
        for i=1:num_refines
           [temp, new_mesh.edge, new_mesh.tri] =...
               refinemesh(mesh.dl, mesh.node,mesh.edge,mesh.tri);
        %    q=pdetriq(p_image,tri); 
        %    figure; pdeplot(p_image,edge,tri,'xydata',q,'colorbar','on','xystyle','flat');
           new_mesh.node = jigglemesh(temp,new_mesh.edge,new_mesh.tri,...
               'opt','mean','iter',inf);
           new_mesh.dl = mesh.dl;
           clear temp; clear mesh; mesh = new_mesh; clear new_mesh;
        %    q=pdetriq(p_image,tri); 
        %    figure; pdeplot(p_image,edge,tri,'xydata',q,'colorbar','on','xystyle','flat');   
        end
    case 2
        % method 2 -- refine only, do not jigglemesh.
        for i = 1:num_refines
            [new_mesh.node,new_mesh.edge,new_mesh.tri] =...
                refinemesh(mesh.dl, mesh.node, mesh.edge, mesh.tri);
            new_mesh.dl = mesh.dl; clear mesh;
            mesh = new_mesh; clear new_mesh;
        end
end
new_mesh = mesh;
return;