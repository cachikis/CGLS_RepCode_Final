%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filename: mps.m
% Author: Craig A. Chikis
% Date: 01/06/2022
% Note(s): 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [] = mps() 
     mvec = repmat({@benchmark_new, @robust_finfric}, [1,1]); 
     model = repmat(["Benchmark", "Financial_Frictions"], [1,1]);
     low_df = repmat([false, true], [1,1]);
     constant_surprise_vec = true(1,2); 
     ar_coef = zeros(length(mvec), 2); 

     mcollect = cell(size(mvec));
     markup_vec_collect = zeros(length(mvec), 4, 7201); 
     quality_collect = zeros(length(mvec), 7201);
     prod_collect = zeros(length(mvec), 7201);
     rhovec_monthly_collect = zeros(length(mvec), 1200);
     chivec_monthly_collect = zeros(length(mvec), 1200);
     rho_path_collect = zeros(length(mvec), 7201);
     rvec_collect = zeros(length(mvec), 2, 7201);
     for oi = 1:length(mvec)

          m = mvec{oi}(); 
          if low_df(oi)
               m.rho = (1.001)^(1/12)-1;
          end

          m = workhorse_nested_robust(m);
          [p10_approx, ~, ~, ~, ~] =  markup_quant(m.kappa, m.lambda, m.s_max, ...
                                                  m.firm_distributions, m.nu_s, 0.1);
          [p50_approx, ~, ~, ~, ~] =  markup_quant(m.kappa, m.lambda, m.s_max, ...
                                                  m.firm_distributions, m.nu_s, 0.5);
          [p90_approx, ~, ~, ~, ~] =  markup_quant(m.kappa, m.lambda, m.s_max, ...
                                                  m.firm_distributions, m.nu_s, 0.9);
          m.p10_markup = p10_approx;
          m.p50_markup = p50_approx;
          m.p90_markup = p90_approx;
          if low_df(oi)
               chi = 0.05;
               rhotilde = 0.05;
               rk = nan;
          else
               chi = 2.143390;
               rk = 5.959426;
               rhotilde = rk - 100*((1+m.growth_rate)^12-1) - chi; 
          end


          x = [0.1, 0.1]; 
          grad = []; 
          tvec = 0:99; 
          rho_vec = zeros(size(tvec)); 
          chi_vec = zeros(size(tvec));
          errvec = zeros(size(x)); 
          method = "fsolve"; 

          function_optimize = @(x) mps_call(x, grad, rho_vec, chi_vec, errvec, tvec, rhotilde, chi, ... 
                                        method, 0.45); 

               
          optimal_ar_coefficients = fsolve(function_optimize, x, ...
                                        optimset('Display', 'off')); 

          ar_coef(oi, :) = optimal_ar_coefficients;

          % Input parameters
          zeta_p = optimal_ar_coefficients(1);
          zeta_chi = optimal_ar_coefficients(2);
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

          for tt = 2:length(tvec)
               rho_vec(tt) = rho_t(rhotilde, zeta_p, rho_vec(tt-1)); 
               chi_vec(tt) = rho_t(chi, zeta_chi, chi_vec(tt-1)); 
          end

          rhovec_monthly = zeros(1, size(tvec, 2)*12);
          chivec_monthly = zeros(size(rhovec_monthly)); 

          count = 1;
          for ii = 1:length(rho_vec)
               if ii < length(rho_vec)
                    rhovec_monthly(count:(count+12-1)) = linspace(rho_vec(ii), rho_vec(ii+1), 12);
                    chivec_monthly(count:(count+12-1)) = linspace(chi_vec(ii), chi_vec(ii+1), 12);
               else
                    rhovec_monthly(count:end) = rho_vec(end);
                    chivec_monthly(count:end) = chi_vec(end); 
               end
               count = count+12; 
          end

          tvec_monthly = zeros(size(rhovec_monthly));
          count = 1; 
          count2 = 0; 
          while count+12-1 <= length(tvec_monthly)
               tvec_monthly(count:(count+12-1)) = linspace(count2, count2+1, 12); 
               count = count+12; 
               count2 = count2+1; 
          end
          tmp = linspace(count2, count2+1, 12);
          tvec_monthly(count:end) = tmp(1:(length(tvec_monthly)-count+1));
          clear tmp; 

          
          
          close all 
          figure; 

          subplot(1,2,1)
          plot(tvec, rho_vec)
          hold on 
          plot(tvec_monthly, rhovec_monthly)

          subplot(1,2,2)
          plot(tvec, chi_vec)
          hold on 
          plot(tvec_monthly, chivec_monthly)


          rhostart = (rhotilde + chi)/100; 
          rhoend = (rhotilde + chi)/100;
          years =  500; 
          feedmu = nan; 
          rho_t = power(1 + (rhovec_monthly + chivec_monthly)/100, 1/12) - 1; 
          justprodlabor = true; 
          no_surprise = ~constant_surprise_vec(oi); 

          [prod, quality, dL_t_dt, ...
                    muder, xder, rdlabor, prodlabor, ...
                    prodlaborder, markup_vec, ...
                    rho_path, tplot, ximpact, ...
                    muimpact, ximpactback, muimpactback, tplotback, ...
                    quality2, ...
                    prod2, ...
                    prod3] = transition_wrapper2(m, rhostart, rhoend, years, ...
                                             justprodlabor, ...
                                             feedmu, no_surprise, ...
                                             rho_t, m.alpha_fric);


          markup_vec = nan(4, size(muimpact, 1)); 
          for ii = 1:(size(muimpact, 1)-1)

               [p10_approx, ~, ~, ~, ~] =  markup_quant(m.kappa, m.lambda, m.s_max, muimpact(ii, :), ...
                                                            m.nu_s, 0.1); 
               [p50_approx, ~, ~, ~, ~] =  markup_quant(m.kappa, m.lambda, m.s_max, muimpact(ii, :), ...
                                                            m.nu_s, 0.5);
               [p90_approx, ~, ~, ~, ~] =  markup_quant(m.kappa, m.lambda, m.s_max, muimpact(ii, :), ...
                                                            m.nu_s, 0.9);
               [~, markup_array, ~, ~, ~] =  markup_quant(m.kappa, m.lambda, m.s_max, muimpact(ii, :), ...
                                                            m.nu_s, 0.5); 


               markup_vec(1, ii) = sum(markup_array(2, :) .* markup_array(1, :)); 
               markup_vec(2, ii) = p10_approx; 
               markup_vec(3, ii) = p50_approx;
               markup_vec(4, ii) = p90_approx; 
          end
          
          rvec = zeros(2, size(muimpact, 1)); 
          rvec(1, :) = rho_path + 100*((1+prod).^12-1);
          rvec(2, :) = rho_path + 100*((1+prod).^12-1) - ...   
                              2*[chi, chivec_monthly, repelem(chi, 12*years)]; 



          mcollect{oi} = m; 
          markup_vec_collect(oi, :, :) = markup_vec; 
          quality_collect(oi, :) = quality; 
          prod_collect(oi, :) = prod; 
          rhovec_monthly_collect(oi, :) = rhovec_monthly;  
          chivec_monthly_collect(oi, :) = chivec_monthly; 
          rho_path_collect(oi, :) = rho_path; 

          rvec_collect(oi, :, :) = rvec; 

     end

     m = benchmark_new();
     m = workhorse_nested_robust(m);
     chi = 2.143390;
     rk = 5.959426;
     rhotilde = rk - 100*((1+m.growth_rate)^12-1) - chi; 
     initrf = zeros(1, 2); 
     initg = zeros(1, 2);
     rfarray = zeros(size(rho_path_collect));
     rhotildevec = zeros(size(rfarray));
     chivec = zeros(size(rfarray));

     initrf(1) = 100*((1+m.growth_rate)^12-1) + rhotilde -chi; 
     initg(1) = 100*((1+m.growth_rate)^12-1);

     rhotildevec(1, :) = [rhovec_monthly_collect(1, :), rhotilde*ones(1, 12*years+1)];
     chivec(1, :) = [chivec_monthly_collect(1, :), chi*ones(1, 12*years+1)]; 

     rfarray(1, :) = rhotildevec(1, :) - chivec(1, :) + 100*((1+prod_collect(1, :)).^12-1); 

     m = robust_finfric();
     m.rho = 1.001^(1/12) - 1; 
     m = workhorse_nested_robust(m);
     chi = 0.05;
     rhotilde = 0.05;

     initrf(2) = 100*((1+m.growth_rate)^12-1) + rhotilde - chi; 
     initg(2) = 100*((1+m.growth_rate)^12-1); 

     rhotildevec(2, :) = [rhovec_monthly_collect(2, :), rhotilde*ones(1,12*years+1)];
     chivec(2, :) = [chivec_monthly_collect(2, :), chi*ones(1,12*years+1)]; 

     rfarray(2, :) = rhotildevec(2, :) - chivec(2, :) + 100*((1+prod_collect(2, :)).^12-1); 

     xlims = [-2;
          10]; 

     for row = 1:size(xlims, 2)
          close all;
          figure; 
          set(gcf, 'PaperUnits', 'inches');
          x_width=8;
          y_width=0.7*3.5;
          set(gcf, 'PaperPosition', [0 0 x_width y_width]); %

          subplot(1,3,1)
          p1 = plot([linspace(-5*12, -0.00000001*12, 5),  0:(length(prod)-1)] ./ 12, ...
                    [0*ones(1,5), rfarray(1, :) - initrf(1)], ...
                    '-b', 'LineWidth', 2);
          hold on 
          p2 = plot([linspace(-5*12, -0.00000001*12, 5),  0:(length(prod)-1)] ./ 12, ...
                    [0*ones(1,5), rfarray(2, :) - initrf(2)], ...
                    '--r', 'LineWidth', 2);
          title({'Risk-free rate, $r^f$'; ''}, 'Interpreter', 'latex')
          xlabel('Years', 'Interpreter', 'latex')
          ylabel('Percentage points', 'Interpreter', 'latex')
          xlim(([xlims(1, row), xlims(2, row)]))
          l1 = legend('Benchmark', "Financial" + newline + "frictions", ... 
                    'Interpreter', 'latex', ...
                    'Location', 'northeast');
          set(l1, 'box', 'off')

          subplot(1,3,2)
          p1 = plot([linspace(-5*12, -0.00000001*12, 5),  0:(length(prod)-1)] ./ 12, ...
                    100*[0*ones(1,5), 100*((1+prod_collect(1, :)).^12-1) - initg(1)], ...
                    '-b', 'LineWidth', 2);
          hold on 
          p2 = plot([linspace(-5*12, -0.00000001*12, 5),  0:(length(prod)-1)] ./ 12, ...
                    100*[0*ones(1,5), 100*((1+prod_collect(2, :)).^12-1) - initg(2)], ...
                    '--r', 'LineWidth', 2);
          title({'Productivity growth, $g$'; ''}, 'Interpreter', 'latex')
          xlabel('Years', 'Interpreter', 'latex')
          ylabel('Basis points', 'Interpreter', 'latex')
          xlim(([xlims(1, row), xlims(2, row)]))


          subplot(1,3,3)
          Ylevel = ones(size(prod_collect, 1), size(prod_collect, 2)+5); 
          Ylevelcounterfactual = ones(size(Ylevel)); 
          for oi = 1:size(Ylevel, 1)
               for ii = 2:size(Ylevel, 2)
                    if ii <= 5
                         Ylevel(oi, ii) = Ylevel(oi, ii-1)*(1 + ((1+initg(oi)/100)^(1/12)-1));
                    else
                         Ylevel(oi, ii) = Ylevel(oi, ii-1)*(1 + prod_collect(oi, ii-5));
                    end
                    Ylevelcounterfactual(oi, ii) = Ylevelcounterfactual(oi, ii-1)*(1 + ((1+initg(oi)/100)^(1/12)-1)); 
               end
          end
          p1 = plot([linspace(-5*12, -0.00000001*12, 5),  0:(length(prod)-1)] ./ 12, ...
                    100*(Ylevel(1, :) ./ Ylevelcounterfactual(1, :) - 1), ...
                    '-b', 'LineWidth', 2);
          hold on 
          p2 = plot([linspace(-5*12, -0.00000001*12, 5),  0:(length(prod)-1)] ./ 12, ...
                    100*(Ylevel(2, :) ./ Ylevelcounterfactual(2, :) - 1), ...
                    '--r', 'LineWidth', 2);
          title({'Productivity level'; ''}, 'Interpreter', 'latex')
          xlabel('Years', 'Interpreter', 'latex')
          ylabel('Percent', 'Interpreter', 'latex')
          xlim([xlims(1, row), xlims(2, row)])


          saveas(gcf, "Output/Figures_Paper/mps_" + row + ".eps",'epsc')
          saveas(gcf, "Output/Figures_Paper/mps_" + row + ".png")
     end
end