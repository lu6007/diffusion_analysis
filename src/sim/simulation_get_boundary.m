% Generate the boundary data for computer_simulation 
% function boundary = simulation_get_boundary(name, data, varargin)
% parameter_name = {'option','load_file'};
% default_value = {1, false};
%
% Exaple and usage: (see computer_simulation.m)
% data.boundary = simulation_get_boundary(cell_name, data, 'load_file', load_file);
%
% Sub-functions:
% function boundary = init_boundary(data)

% Copyright: Shaoying Lu and Yingxiao Wang 2016function boundary = simulation_get_boundary(name, data, varargin)
function boundary = simulation_get_boundary(name, data, varargin)
parameter_name = {'option','load_file'};
default_value = {1, false};
[option, load_file] = parse_parameter(parameter_name, default_value, varargin);
% 1 --- initialize the boundary 
% 2 --- add interface to existing boundary in boundary_file.
% The interface is calculated based on data.diff_map

boundary_file = strcat(data.path, 'boundary.mat');
if exist(boundary_file, 'file') && load_file,
    res = load(boundary_file);
    boundary = res.boundary;
    % rectangle = res.rectangle;
    % Add the subregion interface to the mesh
    % if length(boundary) == 3 && option == 2, % Add Interface
    if length(boundary) == 2 && option == 2, % Add Interface
        switch name,
            case 'layered_diffusion',
                num_subregions = length(data.diff_const);
                cc = 10;
                for i = 2:num_subregions,
                    bd = bwboundaries(data.diff_map == data.diff_const(i), 4, 'noholes');
                    % coarsen by a factor of 10
                    nn = length(bd{1});
                    % bd{1}(nn,:) is different from bd{1}(1, :) 
                    % Turn the interface boundary from counter clockwsie to
                    % clockwise, which is the same between boundary{1} and
                    % boundary{2}. Add the first node to close the
                    % boundary. 
                    boundary{i} = [bd{1}(nn:-cc:2,2), bd{1}(nn:-cc:2,1); bd{1}(1,2), bd{1}(1,1)];
                    clear bd;
               end;
               
            case 'spot_diffusion',
                bd = bwboundaries(data.diff_map == data.diff_const(2), 4, 'noholes');
                num_subregions = length(bd)+1;
                cc = 10;
                % turn the boundary around and close it up
                for i = 2:num_subregions,
                    nn = length(bd{i-1});
                    boundary{i} = [bd{i-1}(nn:-cc:2, 2),bd{i-1}(nn:-cc:2,1); bd{i-1}(1,2), bd{i-1}(1,1)]; 
                end;
                clear bd;
                
%basic idea: D'=T'DT; D' and D have same eigenvalues and in this case are the 
%given diff_const 5 and 100. use this to identify the boundaries of the
%filaments.
            case 'tensor_diffusion', % a specialized case; will have to 
                % change the program if change the diffusion constants or 
                % the orientation of the filaments
                % num_subregions = 3;
                cc = 10;
                a = double(data.diff_map(:,:,1))-2^7;
                b = double(data.diff_map(:,:,2))-2^7;
                c = double(data.diff_map(:,:,3))-2^7;
                mask = zeros(size(a));
                
                % upper filament
                for j = 1:length(data.diff_map(:,1,1))
                    for k = 1:length(data.diff_map(1,:,1))
                        eigv = eig([a(j,k), b(j,k); b(j,k), c(j,k)]);
                        mask(j,k) = (eigv(2) < data.diff_const(2)+2 & eigv(2) > data.diff_const(2) -2 & j < 225); % upper filament
                    end
                end
                bd = bwboundaries(mask,4,'noholes');
                nn = length(bd{1});
                boundary{2} = [bd{1}(nn:-cc:2,2), bd{1}(nn:-cc:2,1);bd{1}(1,2),bd{1}(1,1)];
                clear bd
                
                %lower filament
                for j = 1:length(data.diff_map(:,1,1))
                    for k = 1:length(data.diff_map(1,:,1))
                        eigv = eig([a(j,k), b(j,k); b(j,k), c(j,k)]);
                        mask(j,k) = (eigv(2) < data.diff_const(2)+2 & eigv(2) > data.diff_const(2) -2 & j > 225); % lower filament
                    end
                end
                bd = bwboundaries(mask,4,'noholes');
                nn = length(bd{1});
                boundary{3} = [bd{1}(nn:-cc:2,2), bd{1}(nn:-cc:2,1);bd{1}(1,2),bd{1}(1,1)];
                clear bd

            case 'tensor_cross',
                num_subregions = 3;
                cc = 10;
                a = double(data.diff_map(:,:,1))-2^7;
                b1 = double(data.temp_map(:,:,1))-2^7;
                b2 = double(data.temp_map(:,:,2))-2^7;
                c = double(data.diff_map(:,:,3))-2^7;
                
                mask = zeros(size(a));
                
                
                
                bd = bwboundaries(b1,4,'noholes');
                nn = length(bd{1});
                boundary{2} = [bd{1}(nn:-cc:2,2), bd{1}(nn:-cc:2,1);bd{1}(1,2),bd{1}(1,1)];
                clear bd
                
                bd = bwboundaries(b2,4,'noholes');
                nn = length(bd{1});
                boundary{3} = [bd{1}(nn:-cc:2,2), bd{1}(nn:-cc:2,1);bd{1}(1,2),bd{1}(1,1)];
               
        end; % switch name
        % Check the orientation of the boundaries and interfaces. 
        figure(20); hold on;
        for i = 1:length(boundary),
            if i <=2,
                jj = 2;
            else %if i>=3,
                %jj = 1+cc;
                jj = 2;
            end;
%             for j = 1:size(boundary{i},1),
%                 text(boundary{i}(j,1), boundary{i}(j,2), num2str(j));
%             end;
            plot(boundary{i}(:,1), boundary{i}(:,2), 'r-');
            plot(boundary{i}(1,1), boundary{i}(1,2), 'ko');
            plot(boundary{i}(jj,1), boundary{i}(jj,2), 'k+');
        end;
        %

        
        save(boundary_file, 'boundary');
    end; %if option == 2

else
    switch name,
        case 'square_5_circle',
            rr = data.rectangle;
            t = (1:50:501)'; n = length(t); rt = (501:-50:1)';
            boundary = [ones(n,1), t, ...
                        t(2:n-1), rr(4)*ones(n-2,1), ...
                        rr(3)*ones(n,1), rt, ...
                        rt(2:n-1), ones(n-2,1)];
        case {'photobleach_cell', 'layered_diffusion', 'spot_diffusion', 'tensor_diffusion'},
            % photobleach_cell , layered_diffusion
            boundary = init_boundary(data); 
    end; % switch name,
    save(boundary_file, 'boundary');
end; % if exist(boundary_file, 'file'),

%     figure(20); imagesc(data.image_0); hold on;
%     for i = 1:length(boundary),
%         plot(boundary{i}(:,1), boundary{i}(:,2),'r-');
%     end;

% data.boundary = boundary;
%data.rectangle = rectangle;

return;

% Initialize the boundary file
function boundary = init_boundary(data)
%     temp = fileparts(boundary_file);
%     pa = strcat(temp,'/'); clear temp;
%     image_before = imread(strcat(pa,'cell_before_photobleach.PNG'));
%     rectangle_file = strcat(pa, 'image_rectangle.data');
%     rectangle = get_rectangle(image_before, rectangle_file, 'format', '-ascii');
%     temp = imread(strcat(data.path,'cell_after_photobleach.PNG'));
%     image_0 = imcrop(temp, rectangle); clear temp;
    pa = data.path;
    x1_file = strcat(pa, 'x1.data');
    y1_file = strcat(pa, 'y1.data');
    x2_file = strcat(pa, 'x2.data');
    y2_file = strcat(pa, 'y2.data');
    if ~exist(x1_file,'file'),
        figure; imshow(image_0); caxis auto;
        title('Please mark cell boundary');
        [~,~, ~, x1,y1] = roipoly();
        save(x1_file, 'x1', '-ascii');
        save(y1_file, 'y1', '-ascii');
    else
        x1 = load(x1_file);
        y1 = load(y1_file);
    %                 x2 = load(x2_file);
    %                 y2 = load(y2_file);
    end;
    if ~exist(x2_file,'file'),
        figure; imshow(image_0); caxis auto;
        hold on; plot(x1,y1,'r--');
        title('Please mark ROI boundary')
        [~,~,~,x2,y2] = roipoly();
        save(x2_file, 'x2', '-ascii');
        save(y2_file, 'y2', '-ascii'); 
    else
        x2 = load(x2_file);
        y2 = load(y2_file);
    end;

    boundary{1} = [x1 y1];
    boundary{2} = [x2 y2];

return;