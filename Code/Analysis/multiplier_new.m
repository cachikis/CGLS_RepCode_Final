%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filename: multiplier_new.m
% Author: Craig A. Chikis
% Date: 09/14/2023
% Note(s): 
% MA-MFA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [] = multiplier_new()


s_max = 100;

m = benchmark_new(s_max); 
m = workhorse_nested_robust(m);

chi = 2.143390;
rk = 5.959426;
rhotilde = rk - 100*((1+m.growth_rate)^12-1) - chi; 
chi_vec = (1+linspace(0, 0.05, 200)).^(1/12) - 1; 

mcell = cell(size(chi_vec));
for ii = 1:length(chi_vec)
    mtmp = benchmark_new(s_max);
    mtmp.rho = (1+rhotilde/100)^(1/12)-1 + chi_vec(ii);
    mcell{ii} = mtmp; 
end

parfor ii = 1:length(chi_vec)
    mcell{ii} = workhorse_nested_robust(mcell{ii}); 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
gvec = zeros(size(mcell)); 
xvec = zeros(2*m.s_max+1, size(mcell, 2));
muvec = zeros(m.s_max+1, size(mcell, 2)); 
for ii = 1:length(mcell)
    gvec(ii) = mcell{ii}.growth_rate; 
    xvec(:, ii) = mcell{ii}.investment;
    muvec(:, ii) =  mcell{ii}.firm_distributions; 
end

idx_choose = round(linspace(1, length(gvec), 100));
gspline = spline(chi_vec(idx_choose), gvec(idx_choose));
dgspline = fnder(gspline, 1);
dgdchi = ppval(dgspline, (1+chi/100)^(1/12)-1);

dxdchi = zeros(1, 2*m.s_max+1); 
for ii = -m.s_max:m.s_max 
    xspline = spline(chi_vec(idx_choose), xvec(m.s_max+1+ii, idx_choose));
    xder = fnder(xspline, 1); 
    dxdchi(m.s_max+1+ii) = ppval(xder, (1+chi/100)^(1/12)-1); 
end

dmudchi = zeros(1, m.s_max+1);
for ii = 0:m.s_max 
    muspline = spline(chi_vec(idx_choose), muvec(1+ii, idx_choose));
    muder = fnder(muspline, 1); 
    dmudchi(1+ii) = ppval(muder, (1+chi/100)^(1/12)-1); 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Dvs = zeros(1,2*m.s_max+1);
for (ii = 1:m.s_max)
    Dvs(m.s_max+1+ii) = sum(m.value_functions.*m.prob(m.s_max+1+ii, :)) - m.value_functions(m.s_max+1+ii);
    Dvs(m.s_max+1-ii) = sum(m.value_functions.*m.prob(m.s_max+1-ii, :)) - m.value_functions(m.s_max+1-ii);
end
Dvs(m.s_max+1) = sum(m.value_functions.*m.prob(m.s_max+1,:)) - m.value_functions(m.s_max+1);
Dv_omega = Dvs/m.wage_share;



m.dxdrho = dxdchi;
m.Dv_omega = Dv_omega;

m = PE_sys_eq_w_tp(m);
m = PE_entrant(m);
m = incumbent_entrant_strategic(m);
m = strategic_general_tp(m);
m = Thm1_robust(m, "all");
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
direct_vec = zeros(2*m.s_max+1, m.s_max+1);
strategic_vec = zeros(size(direct_vec));
composition_vec = zeros(size(direct_vec)); 
for ii = 1:size(direct_vec, 2)
    direct_vec(:, ii) = log(m.lambda) * (1 + (ii==1)) * m.firm_distributions(ii) * ...
                                [zeros(1,m.s_max+ii-1), 1, zeros(1,2*m.s_max+1-1-m.s_max-ii+1)]; 
    strategic_vec(:, ii) = log(m.lambda) * (1 + (ii==1)) * m.firm_distributions(ii)*(m.M(m.s_max+2+ii, :) - ...
            [zeros(1,m.s_max+ii-1), 1, zeros(1,2*m.s_max+1-1-m.s_max-ii+1)] );
    composition_vec(:, ii) = log(m.lambda) * (1 + (ii==1)) * m.investment(m.s_max+ii) * m.M(2*m.s_max+2+ii, :); 
end
Mgeffect = (sum(direct_vec + strategic_vec + composition_vec, 2))' * [m.parx_parrho(~isnan(m.parx_parrho)), zeros(1, sum(isnan(m.parx_parrho)))]';

close all 
figure; 

set(gcf, 'PaperUnits', 'inches');
x_width=8.75;
y_width=2.9;
set(gcf, 'PaperPosition', [0 0 x_width y_width]); %


subplot(1,3,1)
p1 = plot(-m.s_max:m.s_max, m.parx_parrho/100, '-r', 'LineWidth', 2);
hold on 
xline(0, 'LineStyle', ':', 'color', [0.75,0.75,0.75], 'LineWidth', 0.5);
ax = ancestor(p1, 'axes');
ax.YAxis.Exponent = 0;
title({'$\mathbf{\partial x}$'; ''}, 'Interpreter', 'latex')
xlabel('Technology position, $\sigma$', 'Interpreter', 'latex')
xlim([-35,35])

subplot(1,3,2)
p2 = plot(-m.s_max:m.s_max, m.M(1, :), '-k', 'LineWidth', 2);
hold on 
yline(0, '-k', 'LineWidth', 1.25);
xline(0, 'LineStyle', ':', 'color', [0.75,0.75,0.75], 'LineWidth', 0.5);
ax = ancestor(p2, 'axes');
ax.YAxis.Exponent = 0;
title({'Growth multiplier, $M_g$'; ''}, 'Interpreter', 'latex')
xlabel('Technology position, $\sigma$', 'Interpreter', 'latex')
xlim([-35,35])

subplot(1,3,3)
p3 = plot(-m.s_max:m.s_max, sum(direct_vec, 2), 'LineStyle', '-', 'color', [0, 0.5, 0], 'LineWidth', 2);
hold on
p4 = plot(-m.s_max:m.s_max, sum(strategic_vec, 2), ':b', 'LineWidth', 2);
p5 = plot(-m.s_max:m.s_max, sum(composition_vec, 2), '--r', 'LineWidth', 2);
xlim([-35,35])
yline(0, '-k', 'LineWidth', 1.25);
xline(0, 'LineStyle', ':', 'color', [0.75,0.75,0.75], 'LineWidth', 0.5);
ax = ancestor(p3, 'axes');
ax.YAxis.Exponent = 0;
xlabel('Technology position, $\sigma$', 'Interpreter', 'latex')
title({'Growth multiplier'; 'components'}, 'Interpreter', 'latex')
l1 = legend([p3,p4,p5],'Direct', 'Strategic', 'Composition', 'Interpreter', 'latex', 'Location', 'northeast');
l1.ItemTokenSize = l1.ItemTokenSize*0.75;
l1.FontSize = l1.FontSize*0.75;
set(l1, 'box', 'off')

saveas(gcf, "Output/Figures_Paper/mult1.eps", 'epsc')
saveas(gcf, "Output/Figures_Paper/mult1.png")


close all 
figure; 

set(gcf, 'PaperUnits', 'inches');
x_width=4;
y_width=3.75;
set(gcf, 'PaperPosition', [0 0 x_width y_width]); %



subplot(2,1,1)
p1 = plot(-m.s_max:m.s_max, m.M(m.s_max+2+4, :), '-k', 'LineWidth', 2);
hold on 
xline(-4, 'LineStyle', ':', 'color', [0.75,0.75,0.75], 'LineWidth', 1.25);
xline(4, 'LineStyle', ':', 'color', [0.75,0.75,0.75], 'LineWidth', 1.25);
xlabel('Technology position, $\sigma$', 'Interpreter', 'latex')
ylabel('Multiplier', 'Interpreter', 'latex')
title({'$M_{x_4}$'; ''}, 'Interpreter', 'latex')
xlim([-11, 11])
xticks([-10, -5 ,0, 5, 10]);

subplot(2,1,2)
p1 = plot(-m.s_max:m.s_max, [zeros(1,m.s_max), 0, 0, 0, 0, 1, zeros(1,2*m.s_max+1-s_max-5)], ...
          'LineStyle', '-', 'color', [0,0.5,0], 'LineWidth', 2);
hold on 
p2 = plot(-m.s_max:m.s_max, m.parx_parxsig(:, m.s_max+1+4), ':b', 'LineWidth', 2);
xline(-4, 'LineStyle', ':', 'color', [0.75,0.75,0.75], 'LineWidth', 1.25);
xline(4, 'LineStyle', ':', 'color', [0.75,0.75,0.75], 'LineWidth', 1.25);
xlabel('Technology position, $\sigma$', 'Interpreter', 'latex')
ylabel('Multiplier', 'Interpreter', 'latex')
title({'$\frac{\partial x_{4}}{\partial x^c_{\sigma}}$'; ''}, 'Interpreter', 'latex')
xlim([-11, 11])
xticks([-10, -5 ,0, 5, 10]);

saveas(gcf, "Output/Figures_Paper/mult2.eps", 'epsc')
saveas(gcf, "Output/Figures_Paper/mult2.png")

% Text
intensive = log(m.lambda)*sum([2, ones(1, m.s_max)] .* m.firm_distributions .* dxdchi((m.s_max+1):end)); 
extensive = log(m.lambda)*sum([2, ones(1, m.s_max)] .* dmudchi .* m.investment((m.s_max+1):end));

s1 = "we obtain that the intensive margin effect is " + num2str(100*intensive, '%0.1f') + ...
        " basis points and the extensive margin effect is " + num2str(100*extensive, '%0.1f') + " basis points.";
writematrix(s1, "Output/LaTeX_Output/IA_intensive_extensive.txt", 'QuoteStrings', false);


s2 = "with only the valuation channel (i.e., absent strategic and composition effects), a 100 basis point rise in the risk premium " + ... 
        "decreases aggregate productivity by " + ...
        "$-[\ln \lambda \sum_{\sigma \in S^+} (1 + \mathbbm{1}_{\sigma=0}) \mu_\sigma e_{\sigma}] \bm{\partial x}$, " + ...
        "or approximately " + num2str(-100 * (sum(direct_vec, 2))' * [m.parx_parrho(~isnan(m.parx_parrho)), zeros(1, sum(isnan(m.parx_parrho)))]', '%0.0f') + ... 
        " basis points. However, in general equilibrium, the decrease in aggregate growth is only about " + ...
        num2str(-100*dgdchi, '%0.0f') + " basis points. Overall, the strategic and composition channels \emph{dampen}, but do not overturn, " + ...
        "the decline in growth coming from the valuation channel. On balance, the strategic channel increases growth by " + ...
        "$[\ln \lambda \sum_{\sigma \in S^+} (1 + \mathbbm{1}_{\sigma=0})  (\mu_\sigma (\mathbb{M}_{x_\sigma}-e_{\sigma}))] \bm{\partial x}$, " + ...
        "or " + num2str(100*  (sum(strategic_vec, 2))' * [m.parx_parrho(~isnan(m.parx_parrho)), zeros(1, sum(isnan(m.parx_parrho)))]', '%0.1f') + ...
        " basis points.  The composition channel decreases growth by  " + ...
        "$-[\ln \lambda \sum_{\sigma \in S^+} (1 + \mathbbm{1}_{\sigma=0})  (\mathbb{M}_{\mu_\sigma} x_\sigma ) ] \bm{\partial x},$ or " + ...
        num2str(-100*  (sum(composition_vec, 2))' * [m.parx_parrho(~isnan(m.parx_parrho)), zeros(1, sum(isnan(m.parx_parrho)))]', '%0.1f') + ...
        " basis points, because of the anti-competitive shift induced by a lower discount rate.";
writematrix(s2, "Output/LaTeX_Output/IA_M_accounting.txt", 'QuoteStrings', false);


end