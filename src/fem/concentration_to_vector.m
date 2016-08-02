% function u = concentration_to_vector(con,p_image)
% parameter_name = {'method'};
% default_value = {1};
%
% covert concentration map to solution vector
% using piecewise linear interpolation
% p_image representing the location of the vectors.
% It can be num_node x 2 or 2 x num_nodes.
%
% 'method' --- 1
% p_image represents the [row col] indices into the matrix con,
% which is the swapped version of the (x,y) location of 
% the nodes on the image [col row].
% 'method' --- 2
% p_image represents the [col row] indices, or the (x,y) location
% of the image con.
% 'interp' --- 'linear' or 'nearest'

% Copyright: Shaoying Lu and Yingxiao Wang 2012
function u = concentration_to_vector(con,p_image,varargin)
% If p_image is num_nodes x 2 
% transpose it.
if size(p_image,1)> size(p_image,2),
    temp = p_image'; clear p_image;
    p_image = temp; clear temp;
end;
% 
parameter_name = {'method','interp'};
default_value = {1, 'linear'};
[method, interp] = parse_parameter(parameter_name, default_value,varargin);
if method == 2,
    temp = [p_image(2,:); p_image(1,:)]; clear p_image;
    p_image = temp; clear temp;
end;
num_nodes = size(p_image, 2);
u = zeros(num_nodes,1);
for k = 1:num_nodes,
    x = p_image(1,k); y = p_image(2,k);
    switch interp,
        case 'linear',
            % linear interpolation
            xi = floor(x); yi = floor(y);
            xj = xi+1; yj = yi+1;
            fii = con(xi,yi); %con: concentration map
            if xj>size(con,1)|| yj>size(con,2), % at the edge of image
                u(k) = fii;
                continue;
            end;
            fij = con(xi,yj); 
            fji = con(xj, yi); fjj= con(xj,yj);
            if (x-xi) + (y-yi) <=1,
                u(k) = fii+(fji-fii)*(x-xi)+(fij-fii)*(y-yi);
            elseif (x-xi) +(y-yi)>1,
                u(k) = fjj+(fjj-fij)*(x-xj)+(fjj-fji)*(y-yj);
            end;
        case 'nearest', 
            % nearest neighbor interpolation
            xi = floor(x+0.5);
            yi = floor(y+0.5);
            if xi>size(con,1),
                xi = xi-1;
            end;
            if yi>size(con,2),
                yi = yi-1;
            end;
            u(k) = con(xi,yi);
    end; % switch interp
            
end;

