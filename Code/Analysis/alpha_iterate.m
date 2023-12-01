%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filename: alpha_iterate.m
% Author: Craig A. Chikis
% Date: 09/11/2023
% Note(s): 
% MA-MFA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [] = alpha_iterate()
alpha_vec = [[0, 0.0025, 0.005, 0.0075, 0.01], 0.01 + [0.0025, 0.005, 0.0075, 0.01]]; 

s_max = 100; 

m = benchmark_new(s_max); 
m = workhorse_nested_robust(m);
chi = 2.143390;
rk = 5.959426;
rhotilde = rk - 100*((1+m.growth_rate)^12-1) - chi; 
chi_vec = sort([linspace(0, 0.0029, 50), (1+chi/100)^(1/12)-1, linspace(0.0029+1e-6, (1+(chi+2)/100)^(1/12)-1, 10) ]); 
rhotilde_vec = sort([linspace(8.46e-4, 0.004, 50+9), (1+rhotilde/100)^(1/12)-1, (1+(rhotilde+2)/100)^(1/12)-1]);

res_collect = cell(length(alpha_vec), length(chi_vec), 2);

iter1 = cell(1, length(chi_vec));
iter2 = cell(1, length(rhotilde_vec)); 
for ii = 1:length(alpha_vec)

    for kk = 1:length(chi_vec)
        mtmp = robust_eta(s_max);
        nu_phi = alpha_vec(ii);
        nu_zeta = 0;
        nu_phi_E = alpha_vec(ii);
        jump_max = mtmp.s_max;
        jump_max_p = mtmp.s_max;
        jump_max_E = mtmp.s_max;
        phi_bar = mtmp.phi_wt;
        zeta_bar = 0;
        phi_E_bar = mtmp.phi_wt_tilde;
        [prob, probE, prob_exog] = kernel_exp_s(phi_bar, zeta_bar, phi_E_bar,...
                                                    nu_phi, nu_zeta, nu_phi_E, ...
                                                    m.s_max, jump_max, jump_max_p, jump_max_E);
        mtmp.prob = prob;
        mtmp.probE = probE;
        mtmp.prob_exog = prob_exog; 
        
        mtmp.rho = chi_vec(kk) + (1+rhotilde/100)^(1/12)-1; 
        iter1{kk} = mtmp; 
    end

    parfor kk = 1:length(chi_vec)
        iter1{kk} = workhorse_nested_robust(iter1{kk});
        iter1{kk}.rf =(1+rhotilde/100)^(1/12)-1 - chi_vec(kk) + iter1{kk}.growth_rate; 
    end

    for kk = 1:length(rhotilde_vec)
        mtmp = robust_eta(s_max);
        nu_phi = alpha_vec(ii);
        nu_zeta = 0;
        nu_phi_E = alpha_vec(ii);
        jump_max = mtmp.s_max;
        jump_max_p = mtmp.s_max;
        jump_max_E = mtmp.s_max;
        phi_bar = mtmp.phi_wt;
        zeta_bar = 0;
        phi_E_bar = mtmp.phi_wt_tilde;
        [prob, probE, prob_exog] = kernel_exp_s(phi_bar, zeta_bar, phi_E_bar,...
                                                    nu_phi, nu_zeta, nu_phi_E, ...
                                                    m.s_max, jump_max, jump_max_p, jump_max_E);
        mtmp.prob = prob;
        mtmp.probE = probE;
        mtmp.prob_exog = prob_exog; 
        
        mtmp.rho = rhotilde_vec(kk) + (1+chi/100)^(1/12)-1; 
        iter2{kk} = mtmp; 
    end

    parfor kk = 1:length(rhotilde_vec)
        iter2{kk} = workhorse_nested_robust(iter2{kk});
        iter2{kk}.rf = rhotilde_vec(kk) - ((1+chi/100)^(1/12)-1) + iter2{kk}.growth_rate;
    end

    res_collect(ii, :, 1) = iter1;
    res_collect(ii, :, 2) = iter2; 

    % disp(ii)


end


gvec = zeros(size(res_collect));
markupvec = zeros(3, size(res_collect, 1), size(res_collect, 2), size(res_collect, 3));
rfvec = zeros(size(res_collect)); 
xvec = zeros(2*mtmp.s_max+1, size(res_collect, 1), size(res_collect, 2), size(res_collect, 3));
muvec = zeros(mtmp.s_max+1, size(res_collect, 1), size(res_collect, 2), size(res_collect, 3)); 
bindingvec = false(size(xvec)); 
for ii = 1:size(gvec, 1)
    for jj = 1:size(gvec, 2)
        for kk = 1:size(gvec, 3)
            gvec(ii,jj,kk) = res_collect{ii,jj,kk}.growth_rate;
            markupvec(1, ii,jj,kk) = res_collect{ii,jj,kk}.markup;
            rfvec(ii,jj,kk) = res_collect{ii,jj,kk}.rf;
            xvec(:, ii, jj, kk) = res_collect{ii,jj,kk}.investment;
            muvec(:, ii, jj, kk) = res_collect{ii,jj,kk}.firm_distributions;
            bindingvec(:, ii,jj,kk) = res_collect{ii,jj,kk}.fric_binding; 

            [p50_approx, ~, ~, ~, ~] =  markup_quant(res_collect{ii,jj,kk}.kappa, ...
                                                     res_collect{ii,jj,kk}.lambda, ...
                                                     res_collect{ii,jj,kk}.s_max, ...
                                                     res_collect{ii,jj,kk}.firm_distributions, ...
                                                     res_collect{ii,jj,kk}.nu_s, 0.5);
            [p90_approx, ~, ~, ~, ~] =  markup_quant(res_collect{ii,jj,kk}.kappa, ...
                                                     res_collect{ii,jj,kk}.lambda, ...
                                                     res_collect{ii,jj,kk}.s_max, ...
                                                     res_collect{ii,jj,kk}.firm_distributions, ...
                                                     res_collect{ii,jj,kk}.nu_s, 0.9);
                                                     
            markupvec(2,ii,jj,kk) = p50_approx-1; 
            markupvec(3,ii,jj,kk) = p90_approx-1; 
        end
    end
end


chi_idx = find(abs(chi_vec - ((1+chi/100)^(1/12)-1)) == min(abs(chi_vec - ((1+chi/100)^(1/12)-1))), 1); 
alpha1_idx = find(abs(alpha_vec - 0.003) == min(abs(alpha_vec - 0.003)), 1); 
alpha2_idx = find(abs(alpha_vec - 0.007) == min(abs(alpha_vec - 0.007)), 1);
alpha3_idx = find(abs(alpha_vec - 0.015) == min(abs(alpha_vec - 0.015)), 1); 
alpha4_idx = find(abs(alpha_vec - 0.02) == min(abs(alpha_vec - 0.02)), 1); 

chi1_idx = find(abs(chi_vec - ((1+1/100)^(1/12)-1)) == min(abs(chi_vec - ((1+1/100)^(1/12)-1))), 1); 
chi3_idx = find(abs(chi_vec - ((1+3/100)^(1/12)-1)) == min(abs(chi_vec - ((1+3/100)^(1/12)-1))), 1); 



close all 
figure; 

set(gcf, 'PaperUnits', 'inches');
x_width=7;
y_width=4.5;
set(gcf, 'PaperPosition', [0 0 x_width y_width]); %


subplot(2,2,1)
p1 = plot(0:mtmp.s_max, 12*squeeze(xvec((mtmp.s_max+1):end, alpha1_idx, chi_idx, 2)), ...
          '--r', 'LineWidth', 2);
hold on 
p3 = plot(0:mtmp.s_max, 12*squeeze(xvec((mtmp.s_max+1):end, alpha3_idx, chi_idx, 2)), ...
          ':b', 'LineWidth', 2);
xlim([0,50])
ylim([0.4, 1])
yticks([0,0.25,0.5,0.75,1])
xlabel('Technology gap, $s$', 'Interpreter', 'latex')
title('Leader innovation rates, $x_{|\sigma|}$', 'Interpreter', 'latex')


subplot(2,2,2)
p1 = plot(0:mtmp.s_max, 12*flip(squeeze(xvec(1:(mtmp.s_max+1), alpha1_idx, chi_idx, 2))), ...
          '--r', 'LineWidth', 2);
hold on 
p3 = plot(0:mtmp.s_max, 12*flip(squeeze(xvec(1:(mtmp.s_max+1), alpha3_idx, chi_idx, 2))), ...
          ':b', 'LineWidth', 2);
xlim([0,50])
ylim([0.1,1])
yticks([0,0.25,0.5,0.75,1])
xlabel('Technology gap, $s$', 'Interpreter', 'latex')
title('Laggard innovation rates, $x_{|\sigma|}$', 'Interpreter', 'latex')
l1 = legend("$\alpha = $" + num2str(alpha_vec(alpha1_idx), '%0.3f'), ... 
            "$\alpha = $" + num2str(alpha_vec(alpha3_idx), '%0.3f'), ...
            'Interpreter', 'latex', 'Location', 'northeast');
set(l1, 'box', 'off')

subplot(2,2,4)
p1 = plot(100*((1+chi_vec).^12-1), 100*((1+squeeze(gvec(alpha1_idx, :, 1))).^12-1), '--r', 'LineWidth', 2);
hold on 
p2 = plot(100*((1+chi_vec).^12-1), 100*((1+squeeze(gvec(alpha2_idx, :, 1))).^12-1), ':b', 'LineWidth', 2);
xlabel('Excess return, $\chi$', 'Interpreter', 'latex')
title('Productivity growth, $g$', 'Interpreter', 'latex')
ylabel('\%', 'Interpreter', 'latex')

subplot(2,2,3)
p1 = plot(0:mtmp.s_max, squeeze(muvec(:, alpha1_idx, chi_idx, 2)), ...
          '--r', 'LineWidth', 2);
hold on 
p3 = plot(0:mtmp.s_max, squeeze(muvec(:, alpha3_idx, chi_idx, 2)), ...
          ':b', 'LineWidth', 2);
xlim([0, 50])
xlabel('Technology gap, $s$', 'Interpreter', 'latex')
title('Distribution of technology gaps, $\mu_s$', 'Interpreter', 'latex')



saveas(gcf, "Output/Figures_Paper/vanishingQCU_main_2.png")
saveas(gcf, "Output/Figures_Paper/vanishingQCU_main_2.eps", 'epsc')



end