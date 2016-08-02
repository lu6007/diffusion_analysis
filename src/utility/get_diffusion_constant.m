% function [diff_const, est_u2, R] = ...
%    get_diffusion_const(u1,u2,M,K,dt)
% Estimate the diffusion constant

% Copyright: Shaoying Lu and Yingxiao Wang 2012
function [diff_const, est_u2, R, r] = ...
    get_diffusion_constant(u1,u2,M,K,dt, is_boundary)

need_gate = 0;
if nargin == 5,
    % set r to []
    r = [];
    % draw_figure = 0;
    Mdu = M*(u2-u1);
    laplace = -dt*K*u2;
    % n = size(Mdu,1);
    % diff_const_1 = (Mdu\laplace);
    % for i=1:n,
    %     if abs(laplace(i))>1.*abs(Mdu(i)),
    %         laplace(i) = diff_const_1*Mdu(i);
    %     end;
    % end;
    
    % apply a gate for including data points
    % The gate is a polygon with 4 nodes.
    % Only the data inside the gate is used for analysis.
    if need_gate,
        figure; plot(Mdu, laplace ,'b+'); hold on;
        title('Please choose 4 points (counterclockwise) to define a selection gate');
        nn = 4;
        xx = zeros(nn+1,1); yy = zeros(nn+1,1);
        for i = 1:nn,
            [xx(i), yy(i)] = ginput(1);
            if i == 1,
                plot(xx(i), yy(i), 'ro');
            else % i>1
                plot(xx(i-1:i), yy(i-1:i), 'r-', 'LineWidth', 0.5);
            end;
        end
        xx(nn+1) = xx(1); yy(nn+1) = yy(1);
        plot(xx(nn:nn+1), yy(nn:nn+1), 'r-','LineWidth', 0.5);
        normal = zeros(length(Mdu), nn);
        index = true(length(Mdu), 1);
        for i = 1:nn,
            % Mdu - x; laplace - y
            % The cross product of (x-xi, y-yi) and (xi+1-xi, yi+1-yi)  has
            % only a z component normal, which is negative for all lines 
            % if (x,y) is inside of the polygon. 
            normal(:,i) = (Mdu-xx(i))*(yy(i+1)-yy(i))-(laplace-yy(i))*(xx(i+1)-xx(i));
            index = index & (normal(:,i)<0);
        end;
        % figure; plot(Mdu(index), laplace(index), 'b+');
        temp = Mdu(index); clear Mdu; Mdu = temp; 
        temp = laplace(index); clear laplace; laplace = temp;
    end;
    
    %
    one_over_D = Mdu\laplace;
    diff_const = 1./one_over_D;
    % diff_const = laplace\Mdu;
    % note laplace\Mdu can be different from 1./(Mdu\laplace)
    G = (M+dt*diff_const*K);
    est_u2 = G\(M*u1);
    [x,y] = average_data(Mdu, laplace);
    [R_matrix,~] = corrcoef(x,y);
    R = R_matrix(1,2);
elseif nargin == 6,
    % use the 2nd method which utilize the bounary information
    [diff_const, est_u2, R, r] = ...
        get_diffusion_constant_2(u1, u2, M, K, dt, is_boundary);
end;
return;

% estimate_diffusion_coefficient
function [diff_const, est_u2, R, r] = ...
    get_diffusion_constant_2(u1,u2,M,K,dt, is_boundary)
    alpha = 0.5; 
    beta = 0.5;
    Mdu = M*(u2-u1);
    dtKu = -dt*K*(alpha*u1+beta*u2);
    one_over_D = Mdu\dtKu;
    diff_const = 1.0/one_over_D;
    % [b, stat] = robustfit(Mdu,dtKu);
    % diff_const = 1.0/b(2);
    r0 = Mdu-diff_const*dtKu;
    r1 = r0.*is_boundary;
    est_u2 = (M+beta*dt*diff_const*K)\...
        ((M-alpha*dt*diff_const*K)*u1+r1);
    [x,y] = average_data(Mdu, dtKu);
    % save average_x.data x -ascii;
    % save average_y.data y -ascii;
    [R,~]=corrcoef(x,y);
    R = R(1,2);
    r = r0-r1;
return;

% function [new_x, new_y] = average_data(x,y)
% average data according to the values of x

% Copyright: Shaoying Lu and Yingxiao Wang 2012
function [new_x, new_y] = average_data(x,y)
min_x = min(x);
max_x = max(x);
num_intervals  = 10;
dx = (max_x-min_x)/num_intervals;
new_x = zeros(num_intervals, 1);
new_y = zeros(num_intervals, 1);
has_data = ones(num_intervals, 1);
for i = 1:num_intervals,
    is_interval = ((x>=min_x+(i-1)*dx) & (x<min_x+i*dx));
    if i == num_intervals,
        is_interval = is_interval+ (x==max_x);
    end;
    n_int = sum(is_interval);
    if n_int ==0,
        has_data(i) = 0;
        continue;
    end;
    x_int = x.* is_interval;
    new_x(i) = sum(x_int)/n_int;
    y_int = y.*is_interval;
    new_y(i) = sum(y_int)/n_int;
end;

old_x = new_x; 
old_y = new_y;
num_data_points = num_intervals - sum(~has_data);
new_x = zeros(num_data_points,1);
new_y = zeros(num_data_points,1);
j = 1;
for i = 1:num_intervals;
    if(has_data(i)),
        new_x(j) = old_x(i);
        new_y(j) = old_y(i);
        j = j+1;
    end;
end;
