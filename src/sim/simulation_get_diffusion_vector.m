% function data = simulation_get_diffusion_vector(cell_name, data, tri_4)
% Initialize the diffusion coefficient vector and related tag, data.diff_coef 
% and data.diff_tag for  computer-generated simulations. 

% Copyright: Shaoying Lu and Yingxiao Wang 2016
function data = simulation_get_diffusion_vector(cell_name, data, tri4)
    num_tris = size(tri4, 2);
    if isfield(data, 'diff_coef')
        diff_coef = data.diff_coef;
        diff_tag = [];
    elseif isfield(data,'diff_vector')
        diff_coef = data.diff_vector;
        diff_tag = [];
    elseif isfield(data,'diff_map')
        
        figure; 
        imshow(data.diff_map); 
        title('Diffusion Map'); 
        colorbar; hold on;

        tc = data.tri_centroid;

        % Isotropic diffusion
        if strcmp(cell_name, 'layered_diffusion')
            diff_coef = zeros(num_tris,1);
            for i = 1:length(data.diff_const)
                ii = ( tri4 == i );
                diff_coef(ii) = data.diff_const(i);
            end
            diff_tag = [];

        elseif strcmp(cell_name, 'spot_diffusion')
            diff_coef = zeros(num_tris,1);
            diff_coef(tri4 == 1) = data.diff_const(1);
            diff_coef(tri4 >= 2) = data.diff_const(2);
            diff_tag = [];
            
        elseif strcmp(cell_name,'tensor_diffusion')
            diff_coef = zeros(num_tris,3);
            D = {(reshape(data.D{1},4,1));...
                 (reshape(data.D{2},4,1))};
            cartesian_diff_const = {[data.diff_const(1), 0, data.diff_const(1)], ...
                                    [D{1}(1), D{1}(2), D{1}(4)],...
                                    [D{2}(1), D{2}(2), D{2}(4)]};

            % set the diffusion coefficient                        
            for i = 1:3
                diff_coef(tri4==1,i)=cartesian_diff_const{1}(i); % the rest of the cell
                diff_coef(tri4==2,i)=cartesian_diff_const{2}(i); % upper filament
                diff_coef(tri4==3,i)=cartesian_diff_const{3}(i); % lower filament
            end
            diff_tag = [];

        elseif strcmp(cell_name, 'tensor_cross_2')
            num_tris = size(tri4, 2);
            diff_coef = zeros(num_tris, 3);
            D = {reshape(data.D{1}, 4, 1); ...
                 reshape(data.D{2}, 4, 1)};
            cartesian_diff_const = {[data.diff_const(1), 0, data.diff_const(1)];...
                                    [D{1}(1), D{1}(2), D{1}(4)];...
                                    [D{2}(1), D{2}(2), D{2}(4)];...
                                    [data.diff_const(2), 0, data.diff_const(2)]};

            % set the diffusion coefficient
            for i = 1 : 3
                diff_coef(tri4 == 1, i) = cartesian_diff_const{1}(i);               % background
                diff_coef((tri4 == 2 | tri4 == 3), i) = cartesian_diff_const{2}(i); % upper filament 
                diff_coef((tri4 == 4 | tri4 == 5), i) = cartesian_diff_const{3}(i); % lower filament
                diff_coef(tri4 == 6, i) = cartesian_diff_const{4}(i);               % cross area
            end
            diff_tag = [];
            
        elseif strcmp(cell_name,'tensor_cross')
            diff_map = double(data.diff_map)-128;
            num_tris = length(tc);
            diff_coef = zeros(num_tris, 3);
            for i = 1:3
               diff_coef(:,i) = concentration_to_vector(diff_map(:,:,i), [tc(2,:); tc(1,:)],'interp','nearest');
            end
            diff_tag = double(diff_coef(:,1)<=200-128)...
                     + 2*double(diff_coef(:,1)>200-128 & diff_coef(:,1)<=215-128)...
                     + 3*double(diff_coef(:,1)>215-128);

               
        elseif strcmp(cell_name, 'general_diffusion') % need to test 'layered_diffusion'
            % diff_coef in the triangular format 
            diff_coef = concentration_to_vector(data.diff_map,[tc(2,:); tc(1,:)],...
                                                'method',1,'interp','nearest');
            num_diff_const = length(diff_coef);
            diff_tag = (1:num_diff_const)';
            p = data.p_image'; 
            data.diff_coef_node = concentration_to_vector(data.diff_map, [p(:, 2), p(:, 1)]);
        else
            fprintf('\nFunction simulation_get_diffusion_vector(): incorrect cell_name! --- \n\n');
        end % if strcmp(cell_name, 'layered_diffusion'),
        clear tc;

    end % elseif isfield(data,'diff_map');

    data.diff_coef = diff_coef;
    data.diff_tag = diff_tag;
    
end

