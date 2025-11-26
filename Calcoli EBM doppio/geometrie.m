clear all;
%close all;
clc;

ri = 8e-3;
ro = 15e-3;
h = 1.5e-3;
tf = 0;
h0 = 2e-3; % [mm] preload

b = sqrt((ro^2-ri^2) + ((h0+h)/2)^2); % base
a = linspace(0, b/2, 100); % altezza arco (freccia)
rr = (ro-ri)/2 + ri;
area = zeros(1, 100);
area_semp = zeros(1, 100);
volume = zeros(1, 100);
volume_semp = zeros(1, 100);
arco = zeros(1,100);
arco_semp = zeros(1, 100);

for i=1:length(area)
    area(i) = area_segmento(b, a(i)); % area
    volume(i) = area(i)*2*pi*rr; % volume
    arco(i) = arco_circolare(b, a(i)); % lunghezza

    % versioni espanse e semplificate
    area_semp(i) = 2/3*b*a(i) + 8/15*a(i)^3/b;
    volume_semp(i) = area_semp(i)*2*pi*rr;
    arco_semp(i) = b + 8*a(i)^2 / (3*b) - 32*a(i)^4 / (15*b^3);
 end

area_max = (pi*(b/2)^2)/2;
volume_max = area_max*2*pi*rr;
arco_max = pi/2 *b;

figure;
plot(a, area, 'DisplayName', 'true');
hold on
plot(a, area_semp, 'DisplayName', 'simplified')
plot(a, area_max*ones(size(a)), 'k--', 'DisplayName', 'Max Area');
legend show;
xlabel('altezza (m)');
ylabel('Area (m^2)');
title('A(a) vs. area max');
grid on;

figure;
plot(a, volume, 'DisplayName', 'true')
hold on
plot(a, volume_semp, 'DisplayName', 'simplified')
plot(a, volume_max*ones(size(a)), 'k--', 'DisplayName', 'Max Volume')
legend show;
xlabel('altezza (m)');
ylabel('Volume (m^3)');
title('V(a) vs. volume max');
grid on;

figure;
plot(a, arco, 'DisplayName', 'true')
hold on
plot(a, arco_semp, 'DisplayName', 'simplified')
plot(a, arco_max*ones(size(a)), 'k--', 'DisplayName', 'Max Arco')
legend show;
xlabel('altezza (m)');
ylabel('Lunghezza arco (m)');
title('l(a) vs. arco max');
grid on;




function A = area_segmento(B,h)
% Area del segmento circolare dato corda B e freccia h (radianti interni)
% Restituisce l'area del segmento "minore" se h<=B/2, altrimenti l'area maggiore.
a = B/2;
R = (a^2 + h^2)/(2*h);
theta = 2*atan2(2*a*h, a^2 - h^2);   % in radianti
Amin = 0.5*R^2*(theta - sin(theta));
if h <= a
    A = Amin;
else
    A = pi*R^2 - Amin;
end
end


function L = arco_circolare(B,h)
% Lunghezza arco di cerchio dato corda B e freccia h
a = B/2;
R = (a^2 + h^2)/(2*h);
theta = 2*atan2(2*a*h, a^2 - h^2);  % radianti, robusto con atan2
L = R*theta;
end


%% Valutazione alpha(h,rc)

h0 = 1.5e-3; % [mm] preload

params.ro = ro;
params.ri = ri;
params.h0 = h0;

% no displacement of EE, bottom actuator zips
alpha_sol(0, 9e-3, params)
% displacement of EE, bottom actuator doesn't zip
alpha_sol(1e-3, ro, params)
% both displacement of EE and zipping
alpha_sol(1e-3, ro-5e-3, params)

% analisi della soluzione: supponiamo che EE si sposti verso il basso di
% massimo 2 mm (fino a "fondo corsa")
% supponiamo che rc vari da ro = 15 mm fino a 9 mm
h_vec = linspace(0, h0, 100)';
rc_vec = linspace(ro, ro - 6e-3, 100)';

% RISULTATI
% scorri colonna per scorrere h
% scorri riga per scorrere rc
alpha_mat = zeros(100, 100);

for i = 1:length(rc_vec)
    rc_val = rc_vec(i);
    alpha_mat(:, i) = alpha_sol(h_vec, rc_val, params); % compute the column
end

% plot alpha(h) per rc fissato
figure;
rc_samples = [15e-3, 13e-3, 11e-3, 9e-3]; % valori campione di rc
hold on; grid on; box on;
for rc_val = rc_samples
    [~, idx] = min(abs(rc_vec - rc_val));
    plot(h_vec*1e3, alpha_mat(:, idx), 'LineWidth', 1.5, ...
        'DisplayName', sprintf('r_c = %.1f mm', rc_val*1e3));
end
xlabel('h [mm]');
ylabel('alpha(h)');
title('alpha(h) per diversi rc fissati');
legend('Location','best');
hold off;

% plot alpha(rc) per h fissato
figure;
h_samples = [0, 0.5e-3, 1e-3, h0]; % valori campione di h
hold on; grid on; box on;
for h_val = h_samples
    [~, idx] = min(abs(h_vec - h_val));
    plot(rc_vec*1e3, alpha_mat(idx, :), 'LineWidth', 1.5, ...
        'DisplayName', sprintf('h = %.1f mm', h_val*1e3));
end
xlabel('rc [mm]');
ylabel('alpha(rc)');
title('alpha(rc) per diversi h fissati');
legend('Location','best');
hold off;



function alpha = alpha_sol(h, rc, params)
% soluzione constraint su alpha da Wolfram Mathematica
ro = params.ro;
ri = params.ri;
h0 = params.h0;

alpha = (sqrt(h.^2 + 2.*h.*h0 + h0.^2 + 4.*(ri - ro).^2).*(-15.*(ri + ro).^2 + ...
    15.^(2./3).*((1./(h.^2 + 2.*h.*h0 + h0.^2 + 4.*(ri - ro).^2)).*(18.*h.*rc.^2.*ri.^2 - 18.*h0.*rc.^2.*ri.^2 + 18.*h.*rc.*ri.^3 - 18.*h0.*rc.*ri.^3 + 36.*h.*rc.^2.*ri.*ro - ...
        36.*h0.*rc.^2.*ri.*ro + 36.*h.*rc.*ri.^2.*ro - 36.*h0.*rc.*ri.^2.*ro - 18.*h.*ri.^3.*ro + 18.*h0.*ri.^3.*ro + 18.*h.*rc.^2.*ro.^2 - 18.*h0.*rc.^2.*ro.^2 + ...
        18.*h.*rc.*ri.*ro.^2 - 18.*h0.*rc.*ri.*ro.^2 - 54.*h.*ri.^2.*ro.^2 + 54.*h0.*ri.^2.*ro.^2 - 54.*h.*ri.*ro.^3 + 54.*h0.*ri.*ro.^3 - 18.*h.*ro.^4 + 18.*h0.*ro.^4 + ...
        sqrt(3).*(h.^2 + 2.*h.*h0 + h0.^2 + 4.*(ri - ro).^2).*sqrt((1./(h.^2 + 2.*h.*h0 + h0.^2 + 4.*(ri - ro).^2).^2).* ...
           ((ri + ro).^4.*(5.*h.^4.*(ri + ro).^2 + 20.*h.^3.*h0.*(ri + ro).^2 + 5.*h0.^4.*(ri + ro).^2 + 80.*(ri - ro).^4.*(ri + ro).^2 + ...
             4.*h.*h0.*(-54.*rc.^4 - 108.*rc.^3.*ri + 108.*rc.*ri.*ro.*(ri + ro) + (ri + ro).^2.*(5.*h0.^2 + 20.*ri.^2 - 40.*ri.*ro - 34.*ro.^2) - 54.*rc.^2.* ...
                (ri.^2 - 2.*ri.*ro - 2.*ro.^2)) + 4.*h0.^2.*(27.*rc.^4 + 54.*rc.^3.*ri - 54.*rc.*ri.*ro.*(ri + ro) + 27.*rc.^2.*(ri.^2 - 2.*ri.*ro - 2.*ro.^2) + ...
               (ri + ro).^2.*(10.*ri.^2 - 20.*ri.*ro + 37.*ro.^2)) + 2.*h.^2.*(54.*rc.^4 + 108.*rc.^3.*ri - 108.*rc.*ri.*ro.*(ri + ro) + 54.*rc.^2.* ...
                (ri.^2 - 2.*ri.*ro - 2.*ro.^2) + (ri + ro).^2.*(15.*h0.^2 + 20.*ri.^2 - 40.*ri.*ro + 74.*ro.^2))))))).^(2./3)))./ ...
  (12.*15.^(1./3).*(ri + ro).*((1./(h.^2 + 2.*h.*h0 + h0.^2 + 4.*(ri - ro).^2)).*(18.*h.*rc.^2.*ri.^2 - 18.*h0.*rc.^2.*ri.^2 + 18.*h.*rc.*ri.^3 - 18.*h0.*rc.*ri.^3 + ...
      36.*h.*rc.^2.*ri.*ro - 36.*h0.*rc.^2.*ri.*ro + 36.*h.*rc.*ri.^2.*ro - 36.*h0.*rc.*ri.^2.*ro - 18.*h.*ri.^3.*ro + 18.*h0.*ri.^3.*ro + 18.*h.*rc.^2.*ro.^2 - ...
      18.*h0.*rc.^2.*ro.^2 + 18.*h.*rc.*ri.*ro.^2 - 18.*h0.*rc.*ri.*ro.^2 - 54.*h.*ri.^2.*ro.^2 + 54.*h0.*ri.^2.*ro.^2 - 54.*h.*ri.*ro.^3 + 54.*h0.*ri.*ro.^3 - 18.*h.*ro.^4 + ...
      18.*h0.*ro.^4 + sqrt(3).*(h.^2 + 2.*h.*h0 + h0.^2 + 4.*(ri - ro).^2).*sqrt((1./(h.^2 + 2.*h.*h0 + h0.^2 + 4.*(ri - ro).^2).^2).* ...
         ((ri + ro).^4.*(5.*h.^4.*(ri + ro).^2 + 20.*h.^3.*h0.*(ri + ro).^2 + 5.*h0.^4.*(ri + ro).^2 + 80.*(ri - ro).^4.*(ri + ro).^2 + ...
           4.*h.*h0.*(-54.*rc.^4 - 108.*rc.^3.*ri + 108.*rc.*ri.*ro.*(ri + ro) + (ri + ro).^2.*(5.*h0.^2 + 20.*ri.^2 - 40.*ri.*ro - 34.*ro.^2) - ...
             54.*rc.^2.*(ri.^2 - 2.*ri.*ro - 2.*ro.^2)) + 4.*h0.^2.*(27.*rc.^4 + 54.*rc.^3.*ri - 54.*rc.*ri.*ro.*(ri + ro) + 27.*rc.^2.*(ri.^2 - 2.*ri.*ro - 2.*ro.^2) + ...
             (ri + ro).^2.*(10.*ri.^2 - 20.*ri.*ro + 37.*ro.^2)) + 2.*h.^2.*(54.*rc.^4 + 108.*rc.^3.*ri - 108.*rc.*ri.*ro.*(ri + ro) + ...
             54.*rc.^2.*(ri.^2 - 2.*ri.*ro - 2.*ro.^2) + (ri + ro).^2.*(15.*h0.^2 + 20.*ri.^2 - 40.*ri.*ro + 74.*ro.^2))))))).^(1./3));
end
