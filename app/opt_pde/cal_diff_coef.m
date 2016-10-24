% Calculate the diffusion coefficient
% with different cases
function [A_sub, M_sub, diff_coef] = cal_diff_coef(cell_name, num_nodes, num_para, d, TT, tri, num_lines, p)

	coef_1 = {'layered_diffusion','spot_diffusion'};
    coef_3 = {'tensor_diffusion','tensor_cross','tensor_cross_2'};


    % Stiffness and mass matrix for different subregion

    if sum(ismember(coef_1, cell_name)) > 0
        diff_coef = ones(1, num_nodes);
        for i = 1 : num_para
            idx_tmp = (tri(4, :) == i);
            tri_tmp = tri(:, idx_tmp);
            diff_coef(idx_tmp) = d(i);
            [A_tmp, M_tmp] = assemble_matrix(p, tri_tmp, 'diff_coef', diff_coef(idx_tmp));
            A_sub{i} = A_tmp;
            M_sub{i} = M_tmp;
        end
				%----------------------
				% i = 1;
				% diff_coef = d;
				% [A_sub{i}, M_sub{i}] = assemble_matrix(p, tri, 'diff_coef', diff_coef);
				%----------------------
    elseif sum(ismember(coef_3, cell_name)) > 0

    	switch cell_name
    		case 'tensor_diffusion'
		        diff_coef  = ones(3, num_nodes);
		        diff_const = d;

		        for i = 1:num_lines,
		            T = TT{i};
		            D{i} = T'*[diff_const(2) 0; 0 diff_const(3)]*T;
		        end

		        % diff_coef = zeros(num_tris,3);
		        D = {[reshape(D{1},4,1)];...
		             [reshape(D{2},4,1)]};
		        cartesian_diff_const = {[diff_const(1),0,diff_const(1)];...
		                                [D{1}(1),D{1}(2),D{1}(4)];...
		                                [D{2}(1),D{2}(2),D{2}(4)]};

		        for i = 1:3
		            diff_coef(i, tri(4,:)==1)=cartesian_diff_const{1}(i); % the rest of the cell
		            diff_coef(i, tri(4,:)==2)=cartesian_diff_const{2}(i); % upper filament
		            diff_coef(i, tri(4,:)==3)=cartesian_diff_const{3}(i); % lower filament
		        end

		        for i = 1 : num_para
		            idx_tmp = (tri(4, :) == i);
		            tri_tmp = tri(:, idx_tmp);
		            [A_tmp, M_tmp] = assemble_matrix(p, tri_tmp, 'diff_coef', diff_coef(:, idx_tmp));
		            A_sub{i} = A_tmp;
		            M_sub{i} = M_tmp;
		        end

    		case 'tensor_cross_2'

	            diff_coef = ones(3, num_nodes);
	            diff_const = d;

	            for i = 1:num_lines,
		            T = TT{i};
		            D{i} = T'*[diff_const(2) 0; 0 diff_const(3)]*T;
		        end

	            D = {[reshape(D{1}, 4, 1)]; ...
	                 [reshape(D{2}, 4, 1)]};
	            cartesian_diff_const = {[diff_const(1), 0, diff_const(1)];...
	                                    [D{1}(1), D{1}(2), D{1}(4)];...
	                                    [D{2}(1), D{2}(2), D{2}(4)];...
	                                    [diff_const(2), 0, diff_const(2)]};

	            % set the diffusion coefficient
	            for i = 1 : 3
	                diff_coef(i, tri(4,:) == 1) = cartesian_diff_const{1}(i); % background
	                diff_coef(i, tri(4,:) == 2) = cartesian_diff_const{2}(i); % upper filament
	                diff_coef(i, tri(4,:) == 3) = cartesian_diff_const{3}(i); % lower filament
	                diff_coef(i, tri(4,:) == 4) = cartesian_diff_const{4}(i); % cross area
	            end

	            for i = 1 : num_para
		            idx_tmp = (tri(4, :) == i);
		            tri_tmp = tri(:, idx_tmp);
		            [A_tmp, M_tmp] = assemble_matrix(p, tri_tmp, 'diff_coef', diff_coef(:, idx_tmp));
		            A_sub{i} = A_tmp;
		            M_sub{i} = M_tmp;
		        end

    	end

    else
        diff_coef = ones(1, num_nodes);
        for i = 1 : num_para
            idx_tmp = (tri(4, :) == i);
            tri_tmp = tri(:, idx_tmp);
            diff_coef(idx_tmp) = d(i);
            [A_tmp, M_tmp] = assemble_matrix(p, tri_tmp, 'diff_coef', diff_coef(idx_tmp));
            A_sub{i} = A_tmp;
            M_sub{i} = M_tmp;
        end
    end
