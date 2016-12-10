function D = F(lx, ly)

    [gx, gy] = meshgrid([1 : lx], [1 : ly]);

    D = exp(sin(gx/10) + cos(gy/10)) * 10;
    % D = 74 - exp(sin(gx/10) + cos(gy/10)) * 10;
