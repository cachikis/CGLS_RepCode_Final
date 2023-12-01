%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filename: mps_call.m
% Author: Craig A. Chikis
% Date: 11/17/2022
% Note(s): 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function errvec = mps_call(x, grad, rho_vec, chi_vec, errvec, tvec, rhotilde, chi, ...
                           method, tenyrtarget)


    % Input parameters
    zeta_p = x(1);
    zeta_chi = x(2);
    rho0 = rhotilde + 1.25;
    chi0 = chi + 0.25;

    % Initialize
    rho_vec(:) = 0;
    chi_vec(:) = 0; 

    % Initialize
    rho_vec(1) = rho0;
    chi_vec(1) = chi0; 

    % AR function
    rho_t = @(rho, zeta_p, rho_tm1) rho + zeta_p*(rho_tm1 - rho); 

    for tt = 2:length(rho_vec)
        rho_vec(tt) = rho_t(rhotilde, zeta_p, rho_vec(tt-1)); 
        chi_vec(tt) = rho_t(chi, zeta_chi, chi_vec(tt-1)); 
    end

    errvec(1) = mean(rho_vec(1:10) - chi_vec(1:10)) - (rhotilde - chi) - tenyrtarget; 
    errvec(2) = 2*(chi_vec(5) - chi) - 0.2*0.50;

    if ~strcmp(method, "fsolve")
        errvec = sqrt(sum(errvec .^ 2));
    end 


end