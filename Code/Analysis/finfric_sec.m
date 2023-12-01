%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filename: finfric_sec.m
% Author: Craig A. Chikis
% Date: 09/11/2022
% Note(s): 
% MA-MFA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [] = fincfric_sec()
beta_vec = sort([[1, 0.05, 0.01, 0.001] / 12, linspace(0.002, 0.1, 10) / 12]); 

m = benchmark_new(); 
m = workhorse_nested_robust(m);
chi = 2.143390;
rk = 5.959426;
rhotilde = rk - 100*((1+m.growth_rate)^12-1) - chi; 
chi_vec = sort([linspace(0, 0.0029, 50), (1+chi/100)^(1/12)-1, linspace(0.0029+1e-6, (1+(chi+2)/100)^(1/12)-1, 10) ]); 
rhotilde_vec = sort([linspace(8.46e-4, 0.004, 50+9), (1+rhotilde/100)^(1/12)-1, (1+(rhotilde+2)/100)^(1/12)-1]);

res_collect = cell(length(beta_vec), length(chi_vec), 4);

iter1 = cell(1, length(chi_vec));
iter2 = cell(1, length(rhotilde_vec)); 
iter3 = cell(size(iter1));
iter4 = cell(size(iter2)); 
for ii = 1:length(beta_vec)

    for kk = 1:length(chi_vec)
        mtmp = robust_finfric();
        mtmp.alpha_fric = beta_vec(ii); 
        mtmp.rho = chi_vec(kk) + (1+rhotilde/100)^(1/12)-1; 
        iter1{kk} = mtmp; 
    end

    parfor kk = 1:length(chi_vec)
        iter1{kk} = workhorse_nested_robust(iter1{kk});
        iter1{kk}.rf =(1+rhotilde/100)^(1/12)-1 - chi_vec(kk) + iter1{kk}.growth_rate; 
    end

    for kk = 1:length(rhotilde_vec)
        mtmp = robust_finfric();
        mtmp.alpha_fric = beta_vec(ii); 
        mtmp.rho = rhotilde_vec(kk) + (1+chi/100)^(1/12)-1; 
        iter2{kk} = mtmp; 
    end

    parfor kk = 1:length(rhotilde_vec)
        iter2{kk} = workhorse_nested_robust(iter2{kk});
        iter2{kk}.rf = rhotilde_vec(kk) - ((1+chi/100)^(1/12)-1) + iter2{kk}.growth_rate;
    end

    for kk = 1:length(chi_vec)
        mtmp = robust_finfric();
        mtmp.alpha_fric = beta_vec(ii); 
        mtmp.rho = chi_vec(kk) + (1.001)^(1/12)-1; 
        iter3{kk} = mtmp; 
    end

    parfor kk = 1:length(chi_vec)
        iter3{kk} = workhorse_nested_robust(iter3{kk});
        iter3{kk}.rf = 1.001^(1/12) -1 - chi_vec(kk) + iter3{kk}.growth_rate; 
    end

    for kk = 1:length(rhotilde_vec)
        mtmp = robust_finfric();
        mtmp.alpha_fric = beta_vec(ii); 
        mtmp.rho = rhotilde_vec(kk) + (1.001)^(1/12)-1; 
        iter4{kk} = mtmp; 
    end

    parfor kk = 1:length(rhotilde_vec)
        iter4{kk} = workhorse_nested_robust(iter4{kk});
        iter4{kk}.rf = rhotilde_vec(kk) - ((1.001)^(1/12)-1) + iter4{kk}.growth_rate;
    end

    res_collect(ii, :, 1) = iter1;
    res_collect(ii, :, 2) = iter2; 
    res_collect(ii, :, 3) = iter3; 
    res_collect(ii, :, 4) = iter4; 



end

chi1_idx = find(abs(chi_vec - (1.01^(1/12)-1)) == min(abs(chi_vec - (1.01^(1/12)-1))), 1); 
chi3_idx = find(abs(chi_vec - (1.03^(1/12)-1)) == min(abs(chi_vec - (1.03^(1/12)-1))), 1); 
chi_idx = find(abs(chi_vec - ((1+chi/100)^(1/12)-1)) == min(abs(chi_vec - ((1+chi/100)^(1/12)-1))), 1); 


beta1_idx = find(abs(beta_vec - 0.1/12) == min(abs(beta_vec - 0.1/12)), 1); 
beta2_idx = find(abs(beta_vec - 0.05/12) == min(abs(beta_vec - 0.05/12)), 1);
beta3_idx = find(abs(beta_vec - 0.01/12) == min(abs(beta_vec - 0.01/12)), 1); 



gvec = zeros(size(res_collect));
markupvec = zeros(3, size(res_collect, 1), size(res_collect, 2), size(res_collect, 3));
rfvec = zeros(size(res_collect)); 
xvec = zeros(2*mtmp.s_max+1, size(res_collect, 1), size(res_collect, 2), size(res_collect, 3));
muvec = zeros(mtmp.s_max+1, size(res_collect, 1), size(res_collect, 2), size(res_collect, 3)); 
bindingvec = zeros(size(xvec)); 
for ii = 1:size(gvec, 1)
    for jj = 1:size(gvec, 2)
        for kk = 1:size(gvec, 3)
            gvec(ii,jj,kk) = res_collect{ii,jj,kk}.growth_rate;
            markupvec(1, ii,jj,kk) = res_collect{ii,jj,kk}.markup;
            rfvec(ii,jj,kk) = res_collect{ii,jj,kk}.rf;
            xvec(:, ii, jj, kk) = res_collect{ii,jj,kk}.investment;
            muvec(:, ii, jj, kk) = res_collect{ii,jj,kk}.firm_distributions;
            bindingvec(:,ii,jj,kk) = res_collect{ii,jj,kk}.fric_binding; 
            try
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
            catch
                markupvec(2,ii,jj,kk) = nan;
                markupvec(3,ii,jj,kk) = nan; 
            end
        end
    end
end

close all
figure; 


set(gcf, 'PaperUnits', 'inches');
x_width=7.5;
y_width=2.75;
set(gcf, 'PaperPosition', [0 0 x_width y_width]); %


subplot(1,3,1)
p1 = plot(-mtmp.s_max:mtmp.s_max, squeeze(xvec(:, beta1_idx, chi_idx, 1))*12, '-b', 'LineWidth', 2);
hold on 
p2 = plot(-mtmp.s_max:mtmp.s_max, squeeze(xvec(:, beta2_idx, chi_idx, 1))*12, '--r', 'LineWidth', 2);
p2 = plot(-mtmp.s_max:mtmp.s_max, squeeze(xvec(:, beta3_idx, chi_idx, 1))*12, 'color', [0,0.5,0],  ...
          'LineStyle', ':', 'LineWidth', 2);
xlabel('Technology position, $\sigma$', 'Interpreter', 'latex')
xlim([-10,25])
title({'Firm innovation rates, $x_\sigma$'; ''}, 'Interpreter', 'latex')
l1 = legend("$\beta = $" + num2str(12*beta_vec(beta1_idx), '%0.2f'), ...
       "$\beta = $" + num2str(12*beta_vec(beta2_idx), '%0.2f'), ...
       "$\beta = $" + num2str(12*beta_vec(beta3_idx), '%0.2f'), ...
       'Interpreter', 'latex', 'Location', 'southeast');
set(l1, 'box', 'off')

subplot(1,3,2)
p1 = plot(1-(beta_vec(12*beta_vec >= 0.002))*12, squeeze(markupvec(1, 12*beta_vec >= 0.002, chi1_idx, 1))*100, ...
          'LineStyle', '-', 'color', [0.5,0.5,0.5], 'LineWidth', 2);
hold on 
p2 = plot(1-(beta_vec(12*beta_vec >= 0.002))*12, squeeze(markupvec(1, 12*beta_vec >= 0.002, chi3_idx, 1))*100, ...
          'LineStyle', '--', 'color', 'k', 'LineWidth', 2);
xlabel('Non-pledgeable share, $1-\beta$', 'Interpreter', 'latex')
title({'Mean markup'; ''}, 'Interpreter', 'latex')
xlim(1-flip([0, 0.1]))
l1 = legend("$\chi = $" + num2str(100*((1+chi_vec(chi1_idx))^12-1), '%0.0f') + "\%", ...
            "$\chi = $" + num2str(round(100*((1+chi_vec(chi3_idx))^12-1)), '%0.0f') + "\%", ...
            'Interpreter', 'latex', 'Location', 'northwest');
set(l1, 'box', 'off'); 

subplot(1,3,3)
p1 = plot(100*((1+chi_vec).^12-1), 100*((1+squeeze(gvec(beta1_idx, :, 3))).^12-1), ...
          '-b', 'LineWidth', 2);
hold on 
p2 = plot(100*((1+chi_vec).^12-1), 100*((1+squeeze(gvec(beta2_idx, :, 3))).^12-1), ...
          '--r', 'LineWidth', 2); 
p3 = plot(100*((1+chi_vec).^12-1), 100*((1+squeeze(gvec(beta3_idx, :, 3))).^12-1), ...
          'LineStyle', ':', 'color', [0,0.5,0], 'LineWidth', 2);
set(l1, 'box', 'off')
title({'Productivity growth, $g$'; 'for low $\rho$'}, 'Interpreter', 'latex')
xlabel('Excess return, $\chi$', 'Interpreter', 'latex')


saveas(gcf, "Output/Figures_Paper/finfric_writeup.png")
saveas(gcf, "Output/Figures_Paper/finfric_writeup.eps", 'epsc')
end