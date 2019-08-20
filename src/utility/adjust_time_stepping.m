% function [G, time_stepping] = ...
%    adjust_time_stepping(diff_const,dt,u_0,M,K)
% adjust time stepping 
% until the solution converges
% using the backward Eular method

% Copyright: Shaoying Lu and Yingxiao Wang 2013
function [G, time_stepping] =  adjust_time_stepping(dt,u_0,M,K)
    solution_error = 5.0e-2*dt;
    G = (M+dt*K);
    while 1
        u_1 = G\(M*u_0);
        dtt = dt/2.0;
        H = (M+dtt*K);
        v_1 = H\(M*(H\(M*u_0)));
        err = norm((v_1-u_1)/size(v_1,1));
        if err>solution_error
            dt = dtt;
            clear G;
            G = H; clear H;
            continue;
        else
            time_stepping = dt;
            break;
        end
    end
return
