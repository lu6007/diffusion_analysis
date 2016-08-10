% solve a large sparse quadratic programming problem
% Contains a linear solver and some optimization loops 
% Created by Yiwen Shi, Jul, 2016, UCSD

function [u,d, s_hist, err_hist, linerr_hist] = opt3(u0, u1, u2, d0, dt, num_nodes, tri, p, cell_name)
    % set the parameters
	gamma  = 1e-8;
    tol    = 1e-6;
    jnwtt  = 3e1;
    mxdamp = 20;
    bnorm0 = 0.0;
    blast  = 0.0;
    eeps   = 0.0;
    rp_52  = 1.0; % damping step s for newton's method
    rp_54  = 0.0; % relerr = ||du|| / ||u||
    rp_56  = 0.0; % relres = ||Gk|| / ||G0||
    rp_57  = 1.0; % ratio  = ||Gk|| / ||G(k-1)||
    rp_58  = 0.0; % dnew   = - <Gu du, G>
    rp_59  = 100.0; % bnorm0 = ||G0||
    snew   = 0.0;
    dnew   = 0.0;
    fnew   = 0.0;
    sold   = 0.0;
    dold   = 0.0;
    fold   = 0.0;
    sleft  = 0.0;
    sright = 0.0;

    num_para = size(unique(tri(4, :)), 2);

	% initialization [u; v; d]
	u = u0;
	v = zeros(num_nodes, 1);
	d = d0;

	Jd = zeros(num_para, 1);
    
    % Diffusion coefficient for each node
    coef_1 = {'layered_diffusion','spot_diffusion'};
    coef_3 = {'tensor_diffusion','tensor_cross','tensor_cross_2'};

    % load orientation data
    TT = {};
    num_lines = 0;
    if sum(ismember(coef_3, cell_name)) > 0

        path = '/Users/Yiwen/Desktop/for_yiwen/diffusion_data/diffusion/simulation/';
        orientation_file = strcat(path, cell_name, '/output/orientation.mat');
        ori = load(orientation_file);
        x   = ori.x; 
        y   = ori.y;
        clear ori;

        num_lines = numel(x)/2;
        D = cell(num_lines,1); % diffusion tensor matrix
        for i = 1 : num_lines
            a(i) = x(2*i) - x(2*i-1);
            b(i) = y(2*i) - y(2*i-1);
            % Normalizing
            temp = sqrt(a(i)^2 + b(i)^2);
            a(i) = a(i) / temp; 
            b(i) = b(i) / temp; 
            clear temp;
            % Calculate the diffusion tensor with the transforming matrix T
            TT{i} = [ a(i) b(i);...
                     -b(i) a(i)]';
        end

    end

    % Stiffness and mass matrix for different subregion
    [A_sub, M_sub, diff_coef] = cal_diff_coef(cell_name, num_nodes, num_para, d, TT, tri, num_lines, p);

    % A -- stiffness matris; M -- mass matrix
    [A, M] = assemble_matrix(p, tri, 'diff_coef', diff_coef);
    A = A + M / dt;

    % Evaluate J at point (u0,v0,d0), J = [Ju; Jv; Jd]
    Ju = M * (u - u2) + A * v;
    Jv = A * u - M * u1 / dt;

    dd0 = d - d0;
    for i = 1 : num_para
        Jd(i, 1) = gamma * sum(sum(M_sub{i})) * dd0(i) + u' * A_sub{i} * v;
    end
    
    F = zeros(num_nodes, num_nodes);
    
    for itnum = 1 : jnwtt

        % ---------- blk4 ----------
        % Set up for subregions 
        diff_coef = ones(1, num_nodes);

        [A_sub, M_sub, diff_coef] = cal_diff_coef(cell_name, num_nodes, num_para, d, TT, tri, num_lines, p);


        % calculate stiffness matrix and mass matrix according current settings
        [A, M] = assemble_matrix(p, tri, 'diff_coef', diff_coef);
        A = A + M / dt;
        
        % evaluate J for current [u; v; d]
		Ju = M * (u - u2) + A * v;
		Jv = A * u - M * u1 / dt;

        dd0 = d - d0;
        for i = 1 : num_para
            Jd(i, 1) = gamma * sum(sum(M_sub{i})) * dd0(i) + u' * A_sub{i} * v;
        end

        % format linear system
        for i = 1 : num_para
            Cv(:, i) = A_sub{i} * u;
            Cu(:, i) = A_sub{i} * v;
            DG(i)    = sum(sum(M_sub{i}));
        end
        G = gamma * diag(DG);
        
        % solve linear system
		Jv_bar = A \ Jv;
        n1 = norm(A * Jv_bar - Jv);

        Cv_bar = A \ Cv;
        n2 = norm(A * Cv_bar - Cv);

		Ju_bar = Ju - M * Jv_bar;
		Cu_bar = Cu - M * Cv_bar;

		G_bar = G - Cu' * Cv_bar - Cv_bar' * Cu_bar;
        tmp   = -(Jd - Cu' * Jv_bar - Cv_bar' * Ju_bar);    
		dd    = G_bar \ tmp;
        n3    = norm(G_bar * dd - tmp);

		du  = -(Jv_bar + Cv_bar * dd);
        tmp = -(Ju_bar + Cu_bar * dd);
		dv  = A \ tmp;
        n4  = norm(A * dv - tmp);
        err_hist(:, itnum) = [n1, n2, n3, n4]';

        HH = [M,   A,   Cu; 
              A,   F,   Cv; 
              Cu', Cv', G  ]; 
        pp = [du; dv; dd];
        JJ = [Ju; Jv; Jd]; 
        linerr_hist(itnum) = norm(HH * pp + JJ);
        % ---------- blk4 ---------- 

        % line search loop
        isw  = 0;
        % ---------- tpick ----------
        % norm4
        [isw,bnorm0,blast,eeps,rp_54,rp_56,rp_57,rp_58,rp_59] =...
         norm4(itnum,isw,bnorm0,blast,eeps,rp_54,rp_56,rp_57,rp_58,rp_59,M,A,Cu,Cv,du,dv,dd,Ju,Jv,Jd,u,v,d,G);        

        % cstep
        [isw,rp_52,rp_54,rp_56,rp_57,rp_58,tol,eeps,snew,sold,sleft,sright,dnew,dold,fnew,fold] =...
         cstep(isw,rp_52,rp_54,rp_56,rp_57,rp_58,tol,eeps,snew,sold,sleft,sright,dnew,dold,fnew,fold);

        if(isw == -1)
            u = u + rp_52 * du;
            v = v + rp_52 * dv;
            d = d + rp_52 * dd; 
            continue; % converged, go to next newton step        
        end
        % ---------- tpick ----------

        for iter = 1 : mxdamp
            % ---------- linsys ----------
            % format matrix and right hand side
            d_tmp = d + rp_52 * dd;
            u_tmp = u + rp_52 * du;
            v_tmp = v + rp_52 * dv;

            [A_sub, M_sub, diff_coef] = cal_diff_coef(cell_name, num_nodes, num_para, d_tmp, TT, tri, num_lines, p);
            
            % if sum(ismember(coef_1, cell_name)) > 0
            %     for i = 1 : num_para
            %         idx_tmp = (tri(4, :) == i); 
            %         tri_tmp = tri(:, idx_tmp);
            %         diff_coef(idx_tmp) = d(i);
            %         [A_tmp, M_tmp] = assemble_matrix(p, tri_tmp, 'diff_coef', diff_coef(idx_tmp));
            %         A_sub{i} = A_tmp;
            %         M_sub{i} = M_tmp;
            %     end
            % elseif sum(ismember(coef_3, cell_name)) > 0


            %     diff_coef  = ones(3, num_nodes);
            %     diff_const = d;

            %     for i = 1:num_lines
            %         T = TT{i};
            %         D{i} = T'*[diff_const(2) 0; 0 diff_const(3)]*T; % D{i} is symmetric
            %     end

            %     % diff_coef = zeros(num_tris,3);
            %     D = {[reshape(D{1},4,1)];[reshape(D{2},4,1)]};
            %     cartesian_diff_const = {[diff_const(1),0,diff_const(1)];...
            %                             [D{1}(1),D{1}(2),D{1}(4)];...
            %                             [D{2}(1),D{2}(2),D{2}(4)]};

            %     for i = 1:3
            %         diff_coef(i, tri(4,:)==1)=cartesian_diff_const{1}(i); % the rest of the cell
            %         diff_coef(i, tri(4,:)==2)=cartesian_diff_const{2}(i); % upper filament
            %         diff_coef(i, tri(4,:)==3)=cartesian_diff_const{3}(i); % lower filament
            %     end

            %     for i = 1 : num_para
            %         idx_tmp = (tri(4, :) == i); 
            %         tri_tmp = tri(:, idx_tmp);
            %         [A_tmp, M_tmp] = assemble_matrix(p, tri_tmp, 'diff_coef', diff_coef(:, idx_tmp));
            %         A_sub{i} = A_tmp;
            %         M_sub{i} = M_tmp;
            %     end

            % else
            %     for i = 1 : num_para
            %         idx_tmp = (tri(4, :) == i); 
            %         tri_tmp = tri(:, idx_tmp);
            %         diff_coef(idx_tmp) = d(i);
            %         [A_tmp, M_tmp] = assemble_matrix(p, tri_tmp, 'diff_coef', diff_coef(idx_tmp));
            %         A_sub{i} = A_tmp;
            %         M_sub{i} = M_tmp;
            %     end
            % end


            [A, M] = assemble_matrix(p, tri, 'diff_coef', diff_coef);
            A = A + M / dt;

            for i = 1 : num_para
                Cv(:, i) = A_sub{i} * u;
                Cu(:, i) = A_sub{i} * v;
                DG(i)    = sum(sum(M_sub{i}));
            end
            G = gamma * diag(DG);
            
            
            Ju = M * (u_tmp - u2) + A * v_tmp;
            Jv = A * u_tmp - M * u1 / dt;
            
            dd0 = d_tmp - d0;
            for i = 1 : num_para
                Jd(i, 1) = gamma * sum(sum(M_sub{i})) * dd0(i) + u' * A_sub{i} * v;
            end
            % ---------- linsys ----------

            % ---------- tpick ----------
            % norm4 
            [isw,bnorm0,blast,eeps,rp_54,rp_56,rp_57,rp_58,rp_59] =...
             norm4(itnum,isw,bnorm0,blast,eeps,rp_54,rp_56,rp_57,rp_58,rp_59,M,A,Cu,Cv,du,dv,dd,Ju,Jv,Jd,u_tmp,v_tmp,d_tmp,G);        

            % cstep
            [isw,rp_52,rp_54,rp_56,rp_57,rp_58,tol,eeps,snew,sold,sleft,sright,dnew,dold,fnew,fold] =...
             cstep(isw,rp_52,rp_54,rp_56,rp_57,rp_58,tol,eeps,snew,sold,sleft,sright,dnew,dold,fnew,fold);

            if(isw == -1)
                u = u + rp_52 * du;
                v = v + rp_52 * dv;
                d = d + rp_52 * dd; 
                continue; % converged, go to next newton step        
            end
            % ---------- tpick ----------

            if(isw >= 0) % need update
                if(iter < mxdamp)
                    break;
                else
                    fprintf('Newton damp reach max iter but still need update!!! \n');
                    break;
                end
            else 
                break;
            end
        end
        
        s_hist(itnum) = rp_52;
        
    end





