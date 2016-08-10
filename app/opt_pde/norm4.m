% Calculate some related norms which will be used later for optimization
% Refer to Prof Bank's pltmg package
% By Yiwen Shi, 07/2016, UCSD

function [isw,bnorm0,blast,eeps,rp_54,rp_56,rp_57,rp_58,rp_59] = norm4(itnum,isw,bnorm0,blast,eeps,rp_54,rp_56,rp_57,rp_58,rp_59,M,A,Cu,Cv,du,dv,dd,Ju,Jv,Jd,u,v,d,G)

    hdu = M * du + A * dv + Cu * dd;
    adu = A * du + Cv * dd;

    bnorm = norm([Ju; Jv; Jd]);
    if (isw <= 0)
        eeps = 1.0e2;
        uunorm = norm(u);
        umnorm = norm(v);
        eunorm = norm(du);
        emnorm = norm(dv);
        rlnorm = norm(d);
        ecnorm = norm(dd);

        q = (uunorm + umnorm + rlnorm) * 1.0e-3;

        rulerr = 1.0;
        if(uunorm + q > eunorm)
            rulerr = eunorm / (uunorm + q);
        end
        if(uunorm + eunorm <= 0)
            rulerr = 0.0;
        end

        rmlerr = 1.0;
        if(umnorm + q > emnorm)
            rmlerr = emnorm / (umnorm + q);
        end
        if(umnorm + emnorm <= 0.0)
            rmlerr = 0.0;
        end

        rclerr = 1.0;
        if(rlnorm + q > ecnorm)
            rclerr = ecnorm / (rlnorm + q);
        end
        if(rlnorm + ecnorm <= 0.0)
            rclerr = 0.0;
        end

        rp_54 = rulerr + rmlerr + rclerr;
        if(bnorm <= 0.0)
            bnorm = eeps;
        end

        if(itnum == 1)
            bnorm0 = max(bnorm, rp_59);
            rp_59 = bnorm0;
        end

    else
        rp_56 = bnorm / bnorm0;
        rp_57 = bnorm / blast;
    end

    gdrl = G * dd;
    uip  = norm(Jv' * adu);
    umip = norm(Ju' * hdu);
    rlip = norm(Jd' * gdrl);
    
    gdrl = gdrl + Cu' * du + Cv' * dv;
    ddnew = -(uip + umip + rlip) / bnorm0^2;
    rp_58 = ddnew;
    blast = bnorm;

