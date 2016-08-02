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
    if exist(data.diffusion_map_file, 'file') && load_file,
        data.diff_map = imread(data.diffusion_map_file);
        
        orientation_file = strcat(data.path, 'output/orientation.mat');
        if exist(orientation_file,'file')
        ori = load(orientation_file);
        x = ori.x; y = ori.y;
        data.orientation = [x y];
        clear ori;
        end
        
        return;
    end;
    switch cell_name,
        case {'photobleach_cell', 'square_5_circle'},
            % length(diff_const) = 1;
            data.diff_coef = data.diff_const;
        case {'photobleach_cell_2', 'photobleach_cell_2b', ...
                'photobleach_cell_3', 'photobleach_cell_3b'},
            % length(diff_const) = 2;
            tri_centroid = data.tri_centroid;
            is_tri = eval(data.diff_vector_statement);
            data.diff_vector = diff_const(1)*double(is_tri)+diff_const(2)*double(~is_tri);
        case 'layered_diffusion', 
            %length>2
            boundary = data.boundary{1}; % outermost boundary of the cell from simulation_get_boundary.m
            [num_rows, num_cols, ~] = size(data.image_0);
            cell_bw = poly2mask(boundary(:,1), boundary(:,2), num_rows, num_cols);
            % Note 10/13/2014
            % The centroid of 1 triangle falls out of cell_bw
            % The resulting diffusion tag has a zero value in that triangle
            % The can be fixed by dilating cell_bw for a few pixels if
            % needed. 
            % divide into 4 layers 1-4 from outside to inside
            num_layers = 4;
            %num_layers = 5; %03/19/15
            [~, label_layer] = divide_layer(cell_bw, num_layers);
            %label_layer(185,55) = 1; % 03/27/15 temporary fix for the weird triangle problem
            diffusion_map = zeros(num_rows, num_cols);
            for i = 1:num_layers,
                 diffusion_map = diffusion_map+(label_layer==i)*diff_const(i);
            end;  
            data.diff_map = uint8(diffusion_map);
            if save_file, 
                imwrite(data.diff_map, data.diffusion_map_file, 'tiff');
            end;
        case 'spot_diffusion',
            temp = imread(strcat(data.path, 'diffusion_spot.png'));
            im_spot = imcrop(temp, data.rectangle); clear temp;
            [num_rows, num_cols, ~] = size(data.image_0);
            % Assign diff_const(1) to the white region and  diff_const(2)
            % to the grey region
            diffusion_map = zeros(num_rows, num_cols);
            %9/26/2014 
            se = strel('square',3);
            mask = imdilate((max(abs(double(im_spot(:,:,:)) -255), [],3)<=30),se);
            diffusion_map = diffusion_map + mask * diff_const(1)+...
                ((~mask)&(max(abs(double(im_spot(:,:,:))-125), [], 3)<=30))*diff_const(2);
            
            data.diff_map = uint8(diffusion_map);
            if save_file, 
                 imwrite(data.diff_map, data.diffusion_map_file, 'tiff');
            end;
        case 'tensor_diffusion',
            temp = imread(strcat(data.path, 'diffusion_spot.png'));
            im_spot = imcrop(temp, data.rectangle); clear temp;
            [num_rows, num_cols, ~] = size(data.image_0);
            
            % Define diffusion direction on the filaments.
            orientation_file = strcat(data.path, 'output/orientation.mat');
            if ~exist(orientation_file, 'file'), %11/24
                h = figure; imagesc(im_spot); title('Diffusion Map'); color bar; 
                hold on;
                waitfor(msgbox({'Define diffusion direction*:',' ',...
                    '1.Start with the highest filament and go downward to other filaments',' ',...
                    '2.Select two points on the filament to define the diffusion dirention', ' ',...
                    '3.Repeat Step 2 on other filaments',' ','4.Press Return when finished',' ',...
                    '* The filament with an higher leftmost point is defined as an higher filament'}));

                [x,y]=getpts;
                data.orientation = [x,y]; %11/18/14 
                save(orientation_file,'x','y'); 
            else
                ori = load(orientation_file);
                x = ori.x; y = ori.y;
                data.orientation = [x y];
                clear ori;
            end;
%             ORI = data.orientation;
%             disp('     x         y');
%             disp(ORI);    
            % 11/18/14 save the coordinates of the filaments. 
            num_lines=numel(x)/2;
            D = cell(num_lines,1); % diffusion tensor matrix
            for i = 1:num_lines,
                a(i)=x(2*i)-x(2*i-1);
                b(i)=y(2*i)-y(2*i-1);
                % Normalizing
                temp = sqrt(a(i)*a(i)+b(i)*b(i));
                a(i) = a(i)/temp; 
                b(i) = b(i)/temp; clear temp;
                %plot(x(2*i-1:2*i),y(2*i-1:2*i),'r-');
                % Calculate the diffusion tensor with the transforming matrix T
                T = [a(i) b(i); -b(i) a(i)]';
                D{i} = T'*[diff_const(2) 0; 0 diff_const(3)]*T; % D{i} is symmetric
            end
            %[a, b] is a 2*n matrix with each column corresponding to a filament
            
            % In the white part of diffusion_spot.png, diff_coef = [5; 0;
            % 5]. In the gray part, diff_coef = [D11; D12; D22]. 
            %add imdilate 09/27/2014
            diffusion_map = zeros(num_rows, num_cols,3);
            se = strel('square',3);
            tag_matrix = imdilate((max(abs(double(im_spot(:,:,:)) -255), [],3)<=30),se); %dilate the white part
            diffusion_map(:,:,1)= tag_matrix* diff_const(1);
            %diffusion_map(:,:,2) = tag_matrix*0
            diffusion_map(:,:,3)= tag_matrix*diff_const(1);
            temp = diffusion_map; clear diffusion_map;
            mask_nonwhite = ~tag_matrix; %tag_matrix is the dilated white part, so the rest nonwhite part is ~tag_matrix.
            tag_matrix = (mask_nonwhite) & (max(abs(double(im_spot(:,:,:))-125), [], 3)<=30);
            % nonwhite mask AND gray 
            % The first filament lies at the region y<225
            % The 2nd filamen lies at the region y>225
            mask = [ones(225, num_cols); zeros(num_rows-225, num_cols)];
            diffusion_map(:,:,1) = temp(:,:,1) +(tag_matrix.*mask)*D{1}(1,1)+...
                (tag_matrix.*~mask)*D{2}(1,1);
            diffusion_map(:,:,2) = temp(:,:,2)+ (tag_matrix.*mask)*D{1}(1,2)+...
                (tag_matrix.*~mask)*D{2}(1,2);
            diffusion_map(:,:,3) = temp(:,:,3) + (tag_matrix.*mask)*D{1}(2,2)+...
                (tag_matrix.*~mask)*D{2}(2,2);
            % Calculate the boundaries of the filaments
%             bd = bwboundaries(tag_matrix);
%             figure(6); hold on;
%             plot(bd{1}(:,2), bd{1}(:,1), 'w');
%             plot(bd{2}(:,2), bd{2}(:,1), 'w');
%             %
            clear temp tag_matrix;
            data.diff_map = uint8(diffusion_map+2^7);  %make all the numbers positive
            data.D = D; % The field D is only for tensor diffusion
            if save_file,
                imwrite(data.diff_map, data.diffusion_map_file, 'tiff');
            end;

        case 'tensor_cross',
            temp = imread(strcat(data.path, 'diffusion_spot.png'));
            im_spot = imcrop(temp, data.rectangle); clear temp;
            [num_rows, num_cols, ~] = size(data.image_0);
            temp1 = imread(strcat(data.path, 'diffusion_spot1.png'));
            im_spot1 = imcrop(temp1, data.rectangle); clear temp;
            temp2 = imread(strcat(data.path, 'diffusion_spot2.png'));
            im_spot2 = imcrop(temp2, data.rectangle); clear temp;
%            [num_rows, num_cols, ~] = size(data.image_0);
            % Define diffusion direction on the filaments.
            orientation_file = strcat(data.path, 'output/orientation.mat');
            if ~exist(orientation_file, 'file'), %11/24
                h = figure; imagesc(im_spot); title('Diffusion Map'); color bar; 
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
            end;
            ORI = data.orientation;
            disp('     x         y');
            disp(ORI);    

            num_lines=numel(x)/2;
            D = cell(num_lines,1); % diffusion tensor matrix
            for i = 1:num_lines,
                a(i)=x(2*i)-x(2*i-1);
                b(i)=y(2*i)-y(2*i-1);
                % Normalizing
                temp = sqrt(a(i)*a(i)+b(i)*b(i));
                a(i) = a(i)/temp; 
                b(i) = b(i)/temp; clear temp;
                %plot(x(2*i-1:2*i),y(2*i-1:2*i),'r-');
                % Calculate the diffusion tensor with the transforming matrix T
                T = [a(i) b(i); -b(i) a(i)]';
                D{i} = T'*[diff_const(2) 0; 0 diff_const(3)]*T; % D{i} is symmetric
            end
            %[a, b] is a 2*n matrix with each column corresponding to a filament
            
            % In the white part of diffusion_spot.png, diff_coef = [5; 0;
            % 5]. In the gray part, diff_coef = [D11; D12; D22]. 
            %add imdilate 09/27/2014
            diffusion_map = zeros(num_rows, num_cols,3);
            temp_map = zeros(num_rows, num_cols,2);
            se = strel('square',3);
           % tag_matrix = imdilate((max(abs(double(im_spot(:,:,:)) -255), [],3)<=30),se);
            tag_matrix1 = imdilate((max(abs(double(im_spot1(:,:,:)) -255), [],3)<=30),se);
            tag_matrix2 = imdilate((max(abs(double(im_spot2(:,:,:)) -255), [],3)<=30),se);%dilate the white part
             %%%%%%%%%%%%Spots%%%%%%%%%%%%%%%%
%            maskS = imdilate((max(abs(double(im_spotS(:,:,:)) -255), [],3)<=30),se);
%            temp(:,:,2) = temp(:,:,2) + maskS*D{1}(1,2);
%                ((~maskS)&(max(abs(double(im_spotS(:,:,:))-125), [], 3)<=30))*D{2}(1,2);
             
            diffusion_map(:,:,1)= tag_matrix1* diff_const(1);
            %diffusion_map(:,:,2) = tag_matrix*0
            diffusion_map(:,:,3)= tag_matrix1*diff_const(1);
            temp = diffusion_map; clear diffusion_map;
            %mask_nonwhite = ~tag_matrix;
            mask_nonwhite1 = ~tag_matrix1;
            mask_nonwhite2 = ~tag_matrix2;
            %tag_matrix is the dilated white part, so the rest nonwhite part is ~tag_matrix.
            tag_matrix1 = (mask_nonwhite1) & (max(abs(double(im_spot1(:,:,:))-125), [], 3)<=30);
            tag_matrix2 = (mask_nonwhite2) & (max(abs(double(im_spot2(:,:,:))-125), [], 3)<=30);
            tag_matrix3 = (tag_matrix1) & (tag_matrix2);
            % nonwhite mask AND gray 
            % The first filament lies at the region y<225
            % The 2nd filamen lies at the region y>225
            mask = ones(num_rows, num_cols);
            diffusion_map(:,:,1) = temp(:,:,1) +(tag_matrix1.*mask)*D{1}(1,1)+...
                (tag_matrix2.*mask)*D{2}(1,1);
            diffusion_map(:,:,2) = temp(:,:,2)+ ((tag_matrix1&(~tag_matrix3)).*mask)*D{1}(1,2)+...
                ((tag_matrix2&(~tag_matrix3)).*mask)*D{2}(1,2);
            diffusion_map(:,:,3) = temp(:,:,3) + (tag_matrix1.*mask)*D{1}(2,2)+...
                (tag_matrix2.*mask)*D{2}(2,2)+(tag_matrix3.*mask)*100;
            data.temp_map(:,:,1) = temp(:,:,2)+ (tag_matrix1.*mask)*D{1}(1,2);
            data.temp_map(:,:,2) = temp(:,:,2)+ (tag_matrix2.*mask)*D{2}(1,2);
            % Calculate the boundaries of the filaments
%             bd = bwboundaries(tag_matrix);
%             figure(6); hold on;
%             plot(bd{1}(:,2), bd{1}(:,1), 'w');
%             plot(bd{2}(:,2), bd{2}(:,1), 'w');
%             %
            clear temp tag_matrix;
            data.temp_map = uint8(data.temp_map+2^7);
            data.diff_map = uint8(diffusion_map+2^7);  %make all the numbers positive
            data.D = D; % The field D is only for tensor diffusion
            if save_file,
                imwrite(data.diff_map, data.diffusion_map_file, 'tiff');
            end;
            
             
    end; %  switch cell_name,
return;