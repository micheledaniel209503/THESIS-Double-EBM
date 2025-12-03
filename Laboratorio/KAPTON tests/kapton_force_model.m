function F = kapton_force_model(h, Y, par)
% KAPTON_FORCE_MODEL: force as a function of h (full height of ebm), Y and
% parmeters. Force is computed analytically in Mathematica from the
% equations.
ri = par.ri;
ro = par.ro;
t0 = par.t0;
nu = par.nu;

%F = (h.*pi.*(2.*ri + sqrt(h.^2 + 4.*(ri - ro).^2) - 2.*ro).*(ri + ro).*t0.*Y)./(4.*(-1 + nu.^2).*sqrt(h.^2 + 4.*(ri - ro).^2).*(ri - ro));
F = (h.*pi.*(2.*ri + sqrt(h.^2 + 4.*(ri - ro).^2) - 2.*ro).*(ri + ro).*t0.*Y)./(2.*(-1 + nu.^2).*sqrt(h.^2 + 4.*(ri - ro).^2).*(ri - ro));


end