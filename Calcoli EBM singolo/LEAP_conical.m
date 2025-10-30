%% THIS IS NOT USED IN THE SOLUTION_SINGLE_LEAP FILE, i use a local function in the main instead
function [Uel,C,Omega]=LEAP_conical(h,rc,ri,ro,t0,E,ni,epsilon_p,epsilon_o,to_left)
% this function is all MECHANICS --> computes Uel, C, Omega values given a
% value of rc and a vector of h values
% given rc, compute [Uel, C, Omega] for many possible h
% the function only models 1/2 of the LEAP (i.e., polymeric membrane against flat rigid electrode)

uc=-h.^2/(2*(rc-ri)^2)/(1/(rc-ri)-(rc^2+ro^2)/(rc^2-ro^2)/rc);

r1=linspace(ri,rc,20)';
r2=linspace(rc,ro,20)';

Uel=0*h;
C=0*h;
Omega=0*h;
for i = 1:length(h) % for the values of h, given the rc  --> compute the associated [Uel, C, Omega]
    eps1_1=uc(i)/(rc-ri)+(h(i).^2+uc(i)^2)/(2*(rc-ri)^2);
    eps2_1=uc(i)/(rc-ri)*(1-ri./r1);
    eps1_2=-uc(i)*rc/(ro^2-rc^2)*(1+ro^2./r2.^2);
    eps2_2=uc(i)*rc/(ro^2-rc^2)*(-1+ro^2./r2.^2);   
    Uel(i)=pi*t0*E/(1-ni^2)*(trapz(r1,r1.*(eps1_1.^2+eps2_1.^2+2*ni*eps1_1.*eps2_1))+trapz(r2,r2.*(eps1_2.^2+eps2_2.^2+2*ni*eps1_2.*eps2_2)));
    
    eps3_1=-ni/(1-ni)*(eps1_1+eps2_1); 
    eps3_2=-ni/(1-ni)*(eps1_2+eps2_2);    
    w=h(i)./(rc-ri).*(r1-ri);
    Cp_cone=2*pi*epsilon_p*epsilon_o*trapz(r1,(1+eps1_1).*(1+eps2_1)./(1+eps3_1).*r1./(epsilon_p*w+epsilon_o*t0));
    
    Cp=2*pi*epsilon_p/t0*trapz(r2,(1+eps1_2).*(1+eps2_2)./(1+eps3_2).*r2) + Cp_cone*0; % polymer Capacitance
    
    if to_left>0 % there's an oil film
        Co=2*pi*epsilon_o/to_left*trapz(r2,(1+eps1_2).*(1+eps2_2).*r2); % oil capacitance
        C(i)=Cp*Co/(Cp+Co);
    else % no oil film
        C(i)=Cp;
    end
    Omega=pi/3*h*(rc^2+ri^2+ri*rc)+pi*(ro^2-rc^2)*to_left;
end

