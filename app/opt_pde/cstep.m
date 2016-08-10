% Calculate the stepsize for newton damping procedure
% Reference to Prof Bank's pltmg package
% By Yiwen Shi, 07/2016, UCSD

function [isw,rp_52,rp_54,rp_56,rp_57,rp_58,tol,eeps,snew,sold,sleft,sright,dnew,dold,fnew,fold] = cstep(isw,rp_52,rp_54,rp_56,rp_57,rp_58,tol,eeps,snew,sold,sleft,sright,dnew,dold,fnew,fold)

    if(isw <= 0)
        eps = 1.0e2;
        tol = 1.0e-2;
        snew = 0.0;
        sleft = 0.0;
        sright = 0.0;
        dnew = rp_58;
        fnew = rp_56^2 / 2.0;
        step = rp_52;
        ratio = rp_57;
        step = step / (step + (1.0 - step) * ratio / 100.0);

        isw = 1;
        ksw = 0;
        rp_52 = step;

        return
    end

    isw = isw + 1;
    sold = snew;
    snew = rp_52;
    dold = dnew;
    dnew = rp_58;
    fold = fnew;
    fnew = rp_56^2 / 2.0;
    relres = rp_56;
    ratio = rp_57;
    relerr = rp_54;

    if(sright <= 0 || dnew > 0.0 || ksw == 1)
        sright = snew;

        if(dnew <= 0.0)
            ksw = 1;
        else
            ksw = 0;
        end
    else
        sleft = snew;
    end

    % sufficient decrease
    ds = sright - sleft;
    if(ds <= tol && dnew <= 0.0)
        isw = -1;
    end
    if(ratio <= 1.0 - eeps * snew && dnew <= 0.0)
        isw = -1;
    end
    if(min(relerr, relres) <= eeps)
        isw = -1;
    end

    if(isw == -1)
        return
    end

    % bisection step
    rp_52 = (sleft + sright) / 2.0;

    if(ksw == 0)
        % secant step
        if(dold == dnew)
            return
        end
        step = snew - dnew * (snew - sold) / (dnew - dold);
    else
        % cubic interpolation step
        ff = -(fold - fnew) * 6.0 / (sold - snew);
        gg = (dold + dnew);
        a  = ff + gg * 3.0;
        b  = -(ff + 2.0 * (gg + dnew));
        c  = dnew;

        if(snew > sold)
            a = -a;
            b = -b;
            c = -c;
        end

        rr = max(abs(a), abs(b), abs(c)) * eeps;

        % quadratix case
        if(abs(a) < rr)
            % b > 0 for min
            if(b <= rr)
                return
            end
            step = snew - (c / b) * (sold - snew);
        else
            % cubic case
            b = b / (2.0 * a);
            c = c / a;
            discr = b^2 - c;

            if(discr <= 0.0)
                return
            end

            d = sqrt(discr);
            if(b < 0.0)
                % the min occrs for 2*a r+b > 0 (not b/2a above)
                if(a > 0.0)
                    r = -(b - d);
                else
                    r = -c / (b-d);
                end
            else
                if(a < 0.0)
                    r = -(b + d);
                else
                    r = -c / (b+d);
                end
            end

            step = snew + r * (sold - snew);
        end
    end

    % choose alternative
    dl = abs(step - sleft);
    dr = abs(step - sright);

    if(max(dl, dr) <= ds * (1.0 - tol))
        rp_52 = step;
    elseif (dl <= ds * tol)
        rp_52 = sleft + ds * tol;
    elseif (dr <= ds * tol)
        rp_52 = sright - ds * tol;
    end

    return
end

% cstep ----------------- <<<<