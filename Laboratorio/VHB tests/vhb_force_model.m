function F = vhb_force_model(x, mu, par)
% VHB_FORCE_MODEL  Force as function of x and mu, using analytical formulas
% for the derivatives, without using diff() to overcome noise
%
%   x   : [m] displacement vector
%   mu  : [Pa] Neo-Hooke parameter
%   par : struct with geometric parameters

lambda_p = par.lambda_p;
r1r      = par.r1r;
r2r      = par.r2r;
t0       = par.t0;

% unstretched reference geometry
dR0 = (r2r - r1r)/lambda_p;   % side length at rest
A0  = pi*((r2r/lambda_p)^2 - (r1r/lambda_p)^2);  % active area at rest

% side length at x!=0
l = sqrt(x.^2 + (r2r - r1r).^2);

% stretches
lambda1 = l ./ dR0;
lambda2 = lambda_p;

% dPsi/dlambda1
% Psi = mu/2 * (lambda1.^2 + lambda2.^2 + (lambda1.*lambda2).^(-2) - 3);
% dPsi/dlambda1 = mu*(lambda1 - lambda1.^(-3).*lambda2.^(-2));
dPsidlambda1 = mu * ( lambda1 - (lambda1.^-3) .* (lambda2^-2) );

% d lambda1 / dx = 1/dR0 * dlambda1/dx
dl_dx        = x ./ l;              % dl/dx
dlambda1_dx  = dl_dx / dR0;         % d(lambda1)/dx

% F = dU/dx = A0 * t0 * dPsi/dlambda1 * d lambda1/dx
F = A0 * t0 .* dPsidlambda1 .* dlambda1_dx;

end
