% function data = init_diffuson_map(cell_name, data)
% Initialize the Diffusion Coefficient or Diffusion Maps for
% Computer-generated Simulations
%  parameter_name = {'load_file', 'save_file'};
%  default_value = { 0, 1};

% Copyright: Shaoying Lu and Yingxiao Wang 2014-2016
function data = init_diffusion_map(cell_name, data, varargin)
    parameter_name = {'load_file', 'save_file'};
    default_value = { 0, 1};
    [load_file, save_file] = parse_parameter(parameter_name, default_value, varargin);

    diff_const = data.diff_const;
    data.orientation = [];

    if exist(data.diffusion_map_file, 'file') && load_file
        data.diff_map = imread(data.diffusion_map_file);
        orientation_file = strcat(data.path, 'output/orientation.mat');
        if exist(orientation_file,'file')
        ori = load(orientation_file);
        x = ori.x; y = ori.y;
        data.orientation = [x y];
        clear ori;
        end
        return;
    end

    switch cell_name
        case 'test'
            [num_row, num_col, ~] = size(data.image_0);
            % Note 10/13/2014
            % The centroid of 1 triangle falls out of cell_bw
            % The resulting diffusion tag has a zero value in that triangle
            % The can be fixed by dilating cell_bw for a few pixels if
            % needed.

            
            diffusion_map    = F(num_col, num_row);

            % boundary = data.boundary{1}; % outermost boundary of the cell from simulation_get_boundary.m
            % cell_bw = poly2mask(boundary(:,1), boundary(:,2), num_row, num_col);
            % divide into 4 layers 1-4 from outside to inside
            % num_layers = 4;
            % [~, label_layer] = divide_layer(cell_bw, num_layers);
            % for i = 1:num_layers
            %      diffusion_map = diffusion_map + (label_layer==i) * diff_const(i);
            % end;

            data.diff_map = uint8(diffusion_map);

            if save_file
                imwrite(data.diff_map, data.diffusion_map_file, 'tiff');
            end
        case {'photobleach_cell', 'square_5_circle'}
            data.diff_coef = data.diff_const;

        case {'photobleach_cell_2', 'photobleach_cell_2b', ...
              'photobleach_cell_3', 'photobleach_cell_3b'}
            % tri_centroid = data.tri_centroid;
            is_tri = eval(data.diff_vector_statement);
            data.diff_vector = diff_const(1) * double(is_tri)...
                             + diff_const(2) * double(~is_tri);

        case 'layered_diffusion'
            boundary = data.boundary{1}; % outermost boundary of the cell from simulation_get_boundary.m
            [num_row, num_col, ~] = size(data.image_0);
            cell_bw = poly2mask(boundary(:,1), boundary(:,2), num_row, num_col);
            % Note 10/13/2014
            % The centroid of 1 triangle falls out of cell_bw
            % The resulting diffusion tag has a zero value in that triangle
            % The can be fixed by dilating cell_bw for a few pixels if
            % needed.

            % divide into 4 layers 1-4 from outside to inside
            num_layers = 4;

            [~, label_layer] = divide_layer(cell_bw, num_layers);
            diffusion_map    = zeros(num_row, num_col);

            for i = 1:num_layers
                 diffusion_map = diffusion_map + (label_layer==i) * diff_const(i);
            end

            data.diff_map = uint8(diffusion_map);

            if save_file
                imwrite(data.diff_map, data.diffusion_map_file, 'tiff');
            end
            
        case 'general_diffusion'
            % outermost boundary of the cell from simulation_get_boundary.m
            boundary = data.boundary{1}; 
            [num_row, num_col, ~] = size(data.image_0);
            cell_bw = poly2mask(boundary(:,1), boundary(:,2), num_row, num_col);
            % The diffusion map is a function of the spatial variables i
            % and j. 
            diffusion_map = zeros(num_row, num_col);
            for i = 1 : num_row
              for j = 1 : num_col
%                diffusion_map(i, j) = cell_bw(i,j) * sqrt((i - num_row/2)^2 + (j - num_col/2)^2);
               diffusion_map(i, j) = cell_bw(i,j) * sqrt((i - 175)^2 + (j - 375)^2);
%                diffusion_map(i, j) = cell_bw(i,j) * sqrt((i - 250)^2 + (j - 400)^2);
%                % Use sigmoid function: Changing too fast, not used. 
%                xx = ((i-175)^2 + (j-375)^2)/100; 
%                sigmoid = 100 * (1-1/(1 + exp(100 - 4 * xx))); 
%                diffusion_map(i, j) = cell_bw(i, j) * sigmoid; 
              end
            end

            data.diff_map = uint8(diffusion_map);

            if save_file
                imwrite(data.diff_map, data.diffusion_map_file, 'tiff');
            end          


        case 'spot_diffusion'
            temp = imread(strcat(data.path, 'diffusion_spot.png'));
            im_spot = imcrop(temp, data.rectangle); clear temp;
            [num_row, num_col, ~] = size(data.image_0);
            diffusion_map = zeros(num_row, num_col);
            se = strel('square',3);
            mask = imdilate((max(abs(double(im_spot(:,:,:)) -255), [],3)<=30),se);

            diffusion_map = diffusion_map + mask * diff_const(1) + ...
                            ((~mask)&(max(abs(double(im_spot(:,:,:))-125), [], 3)<=30)) * diff_const(2);

            data.diff_map = uint8(diffusion_map);

            if save_file
                 imwrite(data.diff_map, data.diffusion_map_file, 'tiff');
            end

        case 'tensor_diffusion'
            temp = imread(strcat(data.path, 'diffusion_spot.png'));
            im_spot = imcrop(temp, data.rectangle); clear temp;
            [num_row, num_col, ~] = size(data.image_0);

            % Define diffusion direction on the filaments.
            orientation_file = strcat(data.path, 'output/orientation.mat');

            if ~exist(orientation_file, 'file')
                figure; imagesc(im_spot); title('Diffusion Map'); color bar;
                hold on;
                waitfor(msgbox({'Define diffusion direction*:',' ',...
                    '1.Start with the highest filament and go downward to other filaments',' ',...
                    '2.Select two points on the filament to define the diffusion dirention', ' ',...
                    '3.Repeat Step 2 on other filaments',' ','4.Press Return when finished',' ',...
                    '* The filament with an higher leftmost point is defined as an higher filament'}));

                [x,y]=getpts;
                data.orientation = [x,y];
                save(orientation_file,'x','y');
            else
                ori = load(orientation_file);
                x = ori.x; y = ori.y;
                data.orientation = [x y];
                clear ori;
            end

            num_line=numel(x)/2;

            % diffusion tensor matrix, D{i} defines the direction of different filaments
            D = cell(num_line,1);
            a = zeros(num_line, 1); 
            b = zeros(num_line, 1); 
            for i = 1:num_line
                a(i)=x(2*i)-x(2*i-1);
                b(i)=y(2*i)-y(2*i-1);
                % Normalizing
                temp = sqrt(a(i)*a(i)+b(i)*b(i));
                a(i) = a(i)/temp;
                b(i) = b(i)/temp; clear temp;
                % Calculate the diffusion tensor with the transforming matrix T
                T = [a(i) b(i); -b(i) a(i)]';
                D{i} = T'*[diff_const(2) 0; 0 diff_const(3)]*T; % D{i} is symmetric
            end
            %[a, b] is a 2*n matrix with each column corresponding to a filament

            % diffusion_spot.png
            % White part: diff_coef = [5; 0; 5].
            % Gray part : diff_coef = [D11; D12; D22].
            diffusion_map = zeros(num_row, num_col,3);
            se = strel('square',3);
            tag_matrix = imdilate((max(abs(double(im_spot(:,:,:)) -255), [],3)<=30),se); %dilate the white part

            % diffusion map for background
            diffusion_map(:,:,1)= tag_matrix * diff_const(1);
            diffusion_map(:,:,3)= tag_matrix * diff_const(1);
            temp = diffusion_map; clear diffusion_map;

            mask_nonwhite = ~tag_matrix;
            tag_matrix = (mask_nonwhite) & (max(abs(double(im_spot(:,:,:))-125), [], 3)<=30);
            % nonwhite mask AND gray
            % The 1st filament lies at the region y<225
            % The 2nd filament lies at the region y>225
            mask = [ones(225, num_col); zeros(num_row-225, num_col)];
            diffusion_map(:,:,1) = temp(:,:,1) ...
                                 + (tag_matrix.* mask) * D{1}(1,1) ... % upper filament
                                 + (tag_matrix.*~mask) * D{2}(1,1);    % lower filament
            diffusion_map(:,:,2) = temp(:,:,2) ...
                                 + (tag_matrix.* mask) * D{1}(1,2) ... % upper filament
                                 + (tag_matrix.*~mask) * D{2}(1,2);    % lower filament
            diffusion_map(:,:,3) = temp(:,:,3) ...
                                 + (tag_matrix.* mask) * D{1}(2,2) ... % upper filament
                                 + (tag_matrix.*~mask) * D{2}(2,2);    % lower filament
            clear temp tag_matrix;

            data.diff_map = uint8(diffusion_map+2^7);  %make all the numbers positive
            data.D = D;

            if save_file
                imwrite(data.diff_map, data.diffusion_map_file, 'tiff');
            end

        case 'tensor_cross_2'
            temp = imread(strcat(data.path, 'diffusion_spot.png'));
            t1 = imread(strcat(data.path, 't1.png'));
            t2 = imread(strcat(data.path, 't2.png'));
            im_spot = imcrop(temp, data.rectangle); clear temp;
            tt1 = imcrop(t1, data.rectangle); clear t1;
            tt2 = imcrop(t2, data.rectangle); clear t2;
            [num_row, num_col, ~] = size(data.image_0);

            % Define diffusion direction on the filaments.
            orientation_file = strcat(data.path, 'output/orientation.mat');

            if ~exist(orientation_file, 'file') %11/24
                figure; imagesc(im_spot); title('Diffusion Map'); color bar;
                hold on;
                waitfor(msgbox({'Define diffusion direction*:',' ',...
                    '1.Start with the highest filament and go downward to other filaments',' ',...
                    '2.Select two points on the filament to define the diffusion dirention', ' ',...
                    '3.Repeat Step 2 on other filaments',' ','4.Press Return when finished',' ',...
                    '* The filament with an higher leftmost point is defined as an higher filament'}));

                [x,y]=getpts;
                data.orientation = [x,y];
                save(orientation_file,'x','y');
            else
                ori = load(orientation_file);
                x = ori.x; y = ori.y;
                data.orientation = [x y];
                clear ori;
            end

            ORI = data.orientation;
            disp('     x         y');
            disp(ORI);

            num_line=numel(x)/2;

            % diffusion tensor matrix, D{i} defines the direction of different filaments
            D = cell(num_line,1);
            a = zeros(num_line, 1);
            b = zeros(num_line, 1);
            for i = 1:num_line
                a(i)=x(2*i)-x(2*i-1);
                b(i)=y(2*i)-y(2*i-1);
                % Normalizing
                temp = sqrt(a(i)*a(i)+b(i)*b(i));
                a(i) = a(i)/temp;
                b(i) = b(i)/temp; clear temp;
                % Calculate the diffusion tensor with the transforming matrix T
                T = [a(i) b(i); -b(i) a(i)]';
                D{i} = T'*[diff_const(2) 0; 0 diff_const(3)]*T; % D{i} is symmetric
            end
            %[a, b] is a 2*n matrix with each column corresponding to a filament

            % diffusion_spot.png
            % In the white part: diff_coef = [5; 0; 5].
            % In the gray part : diff_coef = [D11; D12; D22].
            diffusion_map = zeros(num_row, num_col,3);
            se = strel('square',3);
            tag_matrix = imdilate((max(abs(double(im_spot(:,:,:)) -255), [],3)<=30),se); %dilate the white part

            % set diffusion map for background area
            diffusion_map(:,:,1)= tag_matrix* diff_const(1);
            diffusion_map(:,:,3)= tag_matrix* diff_const(1);

            temp = diffusion_map; clear diffusion_map;
            mask_nonwhite = ~tag_matrix;

            % matrix for different tensors and the cross area
            tensor1 = (mask_nonwhite) & (max(abs(double(tt1(:,:,:))-125), [], 3)<=30);
            tensor2 = (mask_nonwhite) & (max(abs(double(tt2(:,:,:))-125), [], 3)<=30);
            cross = tensor1 & tensor2;

            % set diffusion map for cross area
            temp(:,:,1)= temp(:,:,1) + cross * diff_const(2);
            temp(:,:,3)= temp(:,:,3) + cross * diff_const(2);

            % set diffusion map for tensor area
            diffusion_map(:,:,1) = temp(:,:,1) ...
                                 + (tensor1-cross) * D{1}(1,1) ... % tensor 1
                                 + (tensor2-cross) * D{2}(1,1);    % tensor 2
            diffusion_map(:,:,2) = temp(:,:,2) ...
                                 + (tensor1-cross) * D{1}(1,2) ... % tensor 1
                                 + (tensor2-cross) * D{2}(1,2);    % tensor 2
            diffusion_map(:,:,3) = temp(:,:,3) ...
                                 + (tensor1-cross) * D{1}(2,2) ... % tensor 1
                                 + (tensor2-cross) * D{2}(2,2);    % tensor 2

            clear temp tag_matrix;

            data.diff_map = uint8(diffusion_map+2^7);  %make all the numbers positive
            data.D = D;

            % save tensor and cross area for later use
            data.tensor1 = tensor1;
            data.tensor2 = tensor2;
            data.cross = cross;

            if save_file
                imwrite(data.diff_map, data.diffusion_map_file, 'tiff');
            end

        case 'tensor_cross' % different attempt for tensor cross case
            temp = imread(strcat(data.path, 'diffusion_spot.png'));
            im_spot = imcrop(temp, data.rectangle); clear temp;
            [num_row, num_col, ~] = size(data.image_0);
            temp1 = imread(strcat(data.path, 'diffusion_spot1.png'));
            im_spot1 = imcrop(temp1, data.rectangle); clear temp;
            temp2 = imread(strcat(data.path, 'diffusion_spot2.png'));
            im_spot2 = imcrop(temp2, data.rectangle); clear temp;

            % Define diffusion direction on the filaments.
            orientation_file = strcat(data.path, 'output/orientation.mat');
            if ~exist(orientation_file, 'file')
                figure; imagesc(im_spot); title('Diffusion Map'); color bar;
                hold on;
                waitfor(msgbox({'Define diffusion direction*:',' ',...
                    '1.Start with the highest filament and go downward to other filaments',' ',...
                    '2.Select two points on the filament to define the diffusion dirention', ' ',...
                    '3.Repeat Step 2 on other filaments',' ','4.Press Return when finished',' ',...
                    '* The filament with an higher leftmost point is defined as an higher filament'}));

                [x,y]=getpts;
                % Save the coordinates of the filaments.
                data.orientation = [x,y];
                save(orientation_file,'x','y');
            else
                ori = load(orientation_file);
                x = ori.x; y = ori.y;
                data.orientation = [x y];
                clear ori;
            end
            ORI = data.orientation;
            disp('     x         y');
            disp(ORI);

            num_line=numel(x)/2;
            D = cell(num_line,1); % diffusion tensor matrix
            a = zeros(num_line, 1);
            b = zeros(num_line, 1);
            for i = 1:num_line
                a(i)=x(2*i)-x(2*i-1);
                b(i)=y(2*i)-y(2*i-1);
                % Normalizing
                temp = sqrt(a(i)*a(i)+b(i)*b(i));
                a(i) = a(i)/temp;
                b(i) = b(i)/temp; clear temp;
                % Calculate the diffusion tensor with the transforming matrix T
                T = [a(i) b(i); -b(i) a(i)]';
                D{i} = T'*[diff_const(2) 0; 0 diff_const(3)]*T; % D{i} is symmetric
            end
            %[a, b] is a 2*n matrix with each column corresponding to a filament

            % diffusion_spot.png
            % In the white part: diff_coef = [5; 0; 5].
            % In the gray part : diff_coef = [D11; D12; D22].
            diffusion_map = zeros(num_row, num_col,3);
            % temp_map = zeros(num_row, num_col,2);
            se = strel('square',3);
            tag_matrix1 = imdilate((max(abs(double(im_spot1(:,:,:)) -255), [],3)<=30),se);
            tag_matrix2 = imdilate((max(abs(double(im_spot2(:,:,:)) -255), [],3)<=30),se);%dilate the white part

            diffusion_map(:,:,1)= tag_matrix1* diff_const(1);
            diffusion_map(:,:,3)= tag_matrix1*diff_const(1);
            temp = diffusion_map; clear diffusion_map;

            mask_nonwhite1 = ~tag_matrix1;
            mask_nonwhite2 = ~tag_matrix2;

            tag_matrix1 = (mask_nonwhite1) & (max(abs(double(im_spot1(:,:,:))-125), [], 3)<=30);
            tag_matrix2 = (mask_nonwhite2) & (max(abs(double(im_spot2(:,:,:))-125), [], 3)<=30);
            tag_matrix3 = (tag_matrix1) & (tag_matrix2);

            mask = ones(num_row, num_col);
            diffusion_map(:,:,1) = temp(:,:,1) ...
                                 + (tag_matrix1.*mask)*D{1}(1,1) ...
                                 + (tag_matrix2.*mask)*D{2}(1,1);
            diffusion_map(:,:,2) = temp(:,:,2) ...
                                 + ((tag_matrix1&(~tag_matrix3)).*mask)*D{1}(1,2) ...
                                 + ((tag_matrix2&(~tag_matrix3)).*mask)*D{2}(1,2);
            diffusion_map(:,:,3) = temp(:,:,3) ...
                                 + (tag_matrix1.*mask)*D{1}(2,2) ...
                                 + (tag_matrix2.*mask)*D{2}(2,2)+(tag_matrix3.*mask)*100;
            data.temp_map(:,:,1) = temp(:,:,2) ...
                                 + (tag_matrix1.*mask)*D{1}(1,2);
            data.temp_map(:,:,2) = temp(:,:,2) ...
                                 + (tag_matrix2.*mask)*D{2}(1,2);

            clear temp tag_matrix;
            data.temp_map = uint8(data.temp_map+2^7);
            data.diff_map = uint8(diffusion_map+2^7);  %make all the numbers positive
            data.D = D;

            if save_file
                imwrite(data.diff_map, data.diffusion_map_file, 'tiff');
            end


    end
return;
