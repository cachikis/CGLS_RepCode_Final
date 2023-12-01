%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filename: error_figure.m
% Author: Craig A. Chikis
% Date: 09/14/2023
% Note(s): 
% MA-MFA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [] = error_figure()


models = {@benchmark_new, @entry_new, @LMS_param_paramcorrect, @AA2012_scu, ...
          @AA2012_scu}; 

m = benchmark_new(); 
m = workhorse_nested_robust(m); 
chi = 2.148960;
rk = 5.959426;
rhotilde = rk - 100*((1+m.growth_rate)^12-1) - chi; 

chivec = sort([(1+linspace(0, 0.05, 25)).^(1/12) - 1, 1.001^(1/12)-1]); 
mcell = cell(1, length(chivec));
mcollect = cell(length(models), length(chivec)); 
phi_wt_exog_vec = [nan(1, 3), 1, 1/2]; 
eta_vec = [nan(1, 3), 0.024/12, 0.5*0.024/12]; 
rho_tilde_vec = (1+[rhotilde, rhotilde, rhotilde, 0.1, 0.1]/100).^(1/12) - 1; 
for ii = 1:size(mcollect, 1)

    for jj = 1:length(chivec)
        if ~ismember(ii, [4,5])
            mtmp = models{ii}();
        else
            mtmp = models{ii}(mcollect{ii-1, 1}.s_max, phi_wt_exog_vec(ii), [0, eta_vec(ii)*ones(1, mcollect{ii-1,1}.s_max)] ); 
        end
        mtmp.rho = rho_tilde_vec(ii) + chivec(jj); 
        mcell{jj} = mtmp; 
    end


    parfor jj = 1:length(chivec)
        mcell{jj} = workhorse_nested_robust(mcell{jj});
    end

    mcollect(ii, :) = mcell;

end

errvec = zeros(size(mcollect)); 
gvec = zeros(size(mcollect)); 
xvec =  zeros(2*mcollect{1,1}.s_max+1, size(mcollect, 1), size(mcollect, 2));
muvec = zeros(mcollect{1,1}.s_max+1, size(mcollect, 1), size(mcollect, 2)); 
for ii = 1:size(errvec, 1)
    for jj = 1:size(errvec, 2)
        errvec(ii, jj) = mcollect{ii,jj}.end_hjb_error_rho;
        gvec(ii, jj) = mcollect{ii,jj}.growth_rate; 
        xvec(:, ii, jj) = mcollect{ii,jj}.investment;
        muvec(:, ii, jj) = mcollect{ii,jj}.firm_distributions; 
    end
end

chi1 = find(abs(chivec - ((1.015)^(1/12)-1)) == min(abs(chivec - ((1.015)^(1/12)-1))), 1);
chi2 = find(abs(chivec - ((1.001)^(1/12)-1)) == min(abs(chivec - ((1.001)^(1/12)-1))), 1);

close all 
figure; 


set(gcf, 'PaperUnits', 'inches');
x_width=6.3;
y_width=2.5;
set(gcf, 'PaperPosition', [0 0 x_width y_width]); %


subplot(1,2,1)
p1 = plot(100*((1+chivec).^12-1), log(errvec(1, :)) ./ log(10), '-k', 'LineWidth', 2);
xlabel('Excess return, $\chi$', 'Interpreter', 'latex')
ylabel('$\log_{10} RMSE_{HJB}$', 'Interpreter', 'latex')
title({'Benchmark'; ''}, 'Interpreter', 'latex')
ylim([-20, -10])
xlim([0, 5])
xticks(0:5)


subplot(1,2,2)
p1 = plot(100*((1+chivec).^12-1), log(errvec(2, :)) ./ log(10), '-k', 'LineWidth', 2);
xlabel('Excess return, $\chi$', 'Interpreter', 'latex')
ylabel('$\log_{10} RMSE_{HJB}$', 'Interpreter', 'latex')
title({'Entry'; ''}, 'Interpreter', 'latex')
ylim([-20, -10])
xlim([0, 5])
xticks(0:5)


saveas(gcf, "Output/Figures_Paper/error_figure_chi_out.png");
saveas(gcf, "Output/Figures_Paper/error_figure_chi_out.eps", 'epsc');

close all 
figure;

set(gcf, 'PaperUnits', 'inches');
x_width=6.3;
y_width=2.5;
set(gcf, 'PaperPosition', [0 0 x_width y_width]); %


subplot(1,2,1)
p1 = plot(100*((1+chivec).^12-1), 100*((1+gvec(4, :)).^12-1), '-k', 'LineWidth', 2);
ylabel('Growth rate, $g$', 'Interpreter','latex')
xlabel('Excess return, $\chi$', 'Interpreter', 'latex')
title({'Acemoglu and Akcigit (2012) Slow Catch-up'; ''}, 'Interpreter', 'latex')
xlim([0,5])
xticks(0:5)

subplot(1,2,2)
p2 = plot(100*((1+chivec(~isnan(gvec(5, :)))).^12-1), 100*((1+gvec(5, ~isnan(gvec(5,:)))).^12-1), '--r', 'LineWidth', 2);
ylabel('Growth rate, $g$', 'Interpreter', 'latex')
xlabel('Excess return, $\chi$', 'Interpreter', 'latex')
title({'Version with severe constraints'; 'on creative destruction'}, 'Interpreter', 'latex')
xlim([0,5])
xticks(0:5)


saveas(gcf, "Output/Figures_Paper/chi_aascu1.png");
saveas(gcf, "Output/Figures_Paper/chi_aascu1.eps", 'epsc');

close all 
figure;


set(gcf, 'PaperUnits', 'inches');
x_width=6.3;
y_width=4.5;
set(gcf, 'PaperPosition', [0 0 x_width y_width]); %


subplot(2,2,1)
p1 = plot(0:mcollect{1,1}.s_max, 12*squeeze(xvec((mcollect{1,1}.s_max+1):end, 5, chi2)), ':b', 'LineWidth', 2);
hold on 
p2 = plot(0:mcollect{1,1}.s_max, 12*flip(squeeze(xvec(1:(mcollect{1,1}.s_max+1), 5, chi2))), 'LineStyle', '-.', ...
          'LineWidth', 2, 'color', [0, 0.5, 0]);
xlabel('Technology position, $\sigma$', 'Interpreter', 'latex');
ylabel('$x_\sigma$', 'Interpreter', 'latex')
title({"Innovation rates ($\chi = " + num2str(100*((1+chivec(chi2))^12-1), '%0.1f') + "\%$)"; ""}, 'Interpreter', 'latex')
xlim([0,25])
ylim([0,6])
l1 = legend('Leader', 'Laggard', 'Interpreter', 'latex', 'Location', 'northeast'); 
set(l1, 'box', 'off')

subplot(2,2,2)
p1 = plot(0:mcollect{1,1}.s_max, 12*squeeze(xvec((mcollect{1,1}.s_max+1):end, 5, chi1)), '--r', 'LineWidth', 2);
hold on 
p2 = plot(0:mcollect{1,1}.s_max, 12*flip(squeeze(xvec(1:(mcollect{1,1}.s_max+1), 5, chi1))), 'LineStyle', '-', ...
          'LineWidth', 2, 'color', [0.75, 0.75, 0.75]);
xlabel('Technology position, $\sigma$', 'Interpreter', 'latex');
ylabel('$x_\sigma$', 'Interpreter', 'latex')
title({"Innovation rates ($\chi = " + num2str(100*((1+chivec(chi1))^12-1), '%0.1f') + "\%$)"; ""}, 'Interpreter', 'latex')
xlim([0,25])
ylim([0,6])
l1 = legend('Leader', 'Laggard', 'Interpreter', 'latex', 'Location', 'northeast'); 
set(l1, 'box', 'off')

subplot(2,2,3)
p1 = plot(-mcollect{1,1}.s_max:mcollect{1,1}.s_max, 12*(squeeze(xvec(:, 5, chi1)) - ...
                                                        squeeze(xvec(:, 5, chi2))), ...
          '-k', 'LineWidth', 2); 
ylabel("$\Delta_{" + num2str(100*((1+chivec(chi2))^12-1), '%0.1f') + "\% \to " + ...
       num2str(100*((1+chivec(chi1))^12-1), '%0.1f') + "\%} x_\sigma$", 'Interpreter', 'latex');
xlabel('Technology position, $\sigma$', 'Interpreter', 'latex')
xlim([-25,25])
title({'Change in innovation rate'; ''}, 'Interpreter', 'latex')


subplot(2,2,4)
p1 = plot(0:mcollect{1,1}.s_max, 12*(squeeze(muvec(:, 5, chi1)) - ...
                                     squeeze(muvec(:, 5, chi2))), ...
          '-k', 'LineWidth', 2); 
ylabel("$\Delta_{" + num2str(100*((1+chivec(chi2))^12-1), '%0.1f') + "\% \to " + ...
       num2str(100*((1+chivec(chi1))^12-1), '%0.1f') + "\%} \mu_s$", 'Interpreter', 'latex');
xlabel('Technology gap, $s$', 'Interpreter', 'latex')
xlim([-25,25])
title({'Change in gap distribution'; ''}, 'Interpreter', 'latex')
xlim([0,10])




saveas(gcf, "Output/Figures_Paper/chi_aascu2.png");
saveas(gcf, "Output/Figures_Paper/chi_aascu2.eps", 'epsc');


end