% function mesh =create_mesh(boundary, varargin);
% para_name = {'boundary2', 'method'};
% default_v = {[], 1};
% method 1 --- create mesh using 1 boundary
% method 2 --- create mesh using multi boundaries
% mesh.node: node points with x and y values
% mesh.edge and mesh.tri defines the mesh.
% mesh.dl is probably the boundary condition

% Copyright: Shaoying Lu and Yingxiao Wang 2012
function mesh = create_mesh(boundary, varargin)
    para_name = {'method'};
    default_v = {1};
    method = parse_parameter(para_name, default_v, varargin);

    % create the geometric description and the initial mesh
    switch method
        case 1
            % method 1 --- create mesh using 1 boundary
            polygon =2; 
            num_nodes = size(boundary,1);
            if boundary(1,1)==boundary(num_nodes,1) &&...
                boundary(1,2)==boundary(num_nodes,2),
                new_boundary = boundary(1:num_nodes-1,:);
                num_nodes = num_nodes-1;
            else
                new_boundary = boundary;
            end
            gd = [polygon; num_nodes; new_boundary(:,1); new_boundary(:,2)];
            mesh.dl = decsg(gd);
            [temp,mesh.edge,mesh.tri] = initmesh(mesh.dl);
        case 2
            % method 2 --- create the mesh using multi boundaries
            % mesh = create_mesh_2(boundary);
            % Require the boundary to be clockwise and closed.
            mesh = create_mesh_2(boundary);
            temp = mesh.node;
    end

    % draw image and plot initial mesh
    % handle=figure; pdemesh(p_image,edge,tri);
    % improve mesh
    % q=pdetriq(p_image,tri); 
    % figure; pdeplot(p_image,edge,tri,'xydata',q,'colorbar','on','xystyle','flat')
    mesh.node=jigglemesh(temp,mesh.edge,mesh.tri,'opt','mean','iter',inf);
    clear temp;
    % q=pdetriq(p_image,tri); 
    % figure; pdeplot(p_image,edge,tri,'xydata',q,'colorbar','on','xystyle','flat');

    % % draw image and plot mesh
    % handle=figure; pdemesh(p_image,edge,tri);
return;

% mesh = create_mesh_2(boundary);
% Require the boundaries to be clockwise and closed.
function mesh = create_mesh_2(boundary)
    
    num_boundary = length(boundary);
    % num_boundary = 2;
    n_boundary = zeros(num_boundary,1);
    for i = 1:num_boundary
        n_boundary(i) = length(boundary{i});
    end
    max_n = max(n_boundary);
    gd = zeros(2*max_n, num_boundary);
    polygon = 2;
    for i = 1:num_boundary
        nn = n_boundary(i) - 1;
        x = floor(boundary{i}(1:nn,1)+0.5); y = floor(boundary{i}(1:nn,2)+0.5);
        % % Make the boundaries convex if needed.
        % dt = delaunayTriangulation(x,y);
        % index = convexHull(dt); clear dt;
        % nn = length(index)-1;
        % gd(1:2*(nn+1), i) = [polygon; nn; x(index(nn:-1:1)); y(index(nn:-1:1))]; 
        % clear index;
      gd(1:2*(nn+1), i) = [polygon; nn; x; y];
    end
    mesh.dl = decsg(gd);
    
    % gd - geometric description
    % sf - set formular
    % ns - name space matrix by ascii code
    % 80 - P, 49 - 1, 50 - 2, 51 - 3
    % sf = 'P1'; 
    % ns = [80 80 80 80  
    %       49 50 51 52]; 
    % mesh.dl = decsg(gd, sf, ns);

    [mesh.node,mesh.edge,mesh.tri] = initmesh(mesh.dl);
    
    % When there are multiple interfaces, initmesh can generate edges with
    % incorrect orientation. The code below correct the orientation problem
    % If the boundaries are clockwise, the subregion inside the edge should be
    % larger that that outside. So check this inequality : edge(6,:)<edge(7,:)
    % If edge(6,:)>edge(7,:), the edge is the wrong orientation
    % Switch the orientation of the edge
    temp = mesh.edge;
    ii = (temp(6,:) > temp(7,:)); % returns 1 if the orientation is incorrect.
    mesh.edge(1,ii) = temp(2,ii);
    mesh.edge(2,ii) = temp(1,ii);
    mesh.edge(3,ii) = temp(4,ii);
    mesh.edge(4,ii) = temp(3,ii);
    mesh.edge(6,ii) = temp(7,ii);
    mesh.edge(7,ii) = temp(6,ii);

return;




