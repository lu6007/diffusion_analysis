% function est_u = simulate_diffusion(u, dt, K, M)

% Copyright: Shaoying Lu and Yingxiao Wang 2012
function est_u = simulate_diffusion(u, dt, K, M, varargin)
    parameter_name = {'method'};
    default_value = {1};
    method = parse_parameter(parameter_name, default_value, varargin);
    switch method
        case 1
            % use the Crank-Nicolson method
            ddt =1.0;
            num_steps = dt/ddt;
            alpha = 0.5;
            beta = 0.5;
            est_u = (M+beta*dt*K)\...
                ((M-alpha*ddt*K)*u(:,1));
            for i = 1:num_steps-1
                est_u = (M+beta*dt*K)\...
                    ((M-alpha*ddt*K)*est_u);
            end
        case 2
            % Backward Eular
            G = M+dt*K;
            est_u = G\(M*u);            
    end
return
