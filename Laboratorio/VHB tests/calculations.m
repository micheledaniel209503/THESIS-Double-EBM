%% VHB TESTS FOR NEO-HOOK 

R0 = 0.5*55e-3; % radius of VHB disc at unstretched config
R1 = 0.5*75e-3; % radius of VHB disc at stretched config
lambda_p = R1/R0; % pre-stretch (uniform) factor
r1r = 0.5*15e-3; % radius of indentor
r2r = R0; % external radius of !!ACTIVE!! VHB disc
dR0 = (r2r - r1r)/lambda_p; % ACTIVE VHB disc at unstretched config
t0 = 0.5e-3; % thickness of unstretched VHB

% x stroke vector
hr = 10e-3; % typical reservoir height
x_max = hr + 5e-3;
x = linspace(0, x_max, 500);

l = sqrt(x.^2 + (r2r-r1r)^2);
lambda_1 = l./dR0;
lambda_2 = lambda_p;

mu = 0.4e5; % [Pa]
Psi_vhb = mu/2*(lambda_1.^2 + lambda_2.^2 + (lambda_1.*lambda_2).^(-2) - 3);
Um_vhb = pi*((r2r/lambda_p)^2 - (r1r/lambda_p)^2)*t0.*Psi_vhb;
F_x = diff(Um_vhb)./diff(x);

figure
plot(x(1:end-1).*1e3, F_x)
xlabel('x [mm]')
ylabel('Force [N]')
