%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filename: transition_solve.m
% Author: Craig A. Chikis
% Date: 09/12/2022
% Note(s):
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [] = transition_solve() 

    m = benchmark_new(); 
    m = workhorse_nested_robust(m); 
    
    chi = 2.143390;
    rk = 5.959426;
    rhotilde = rk - 100*((1+m.growth_rate)^12-1) - chi; 

    years = 500;
    yearsrhochange_vec = [1/12, 10];
    no_surprise_vec = [true, false]; 
    rhoend = 0.03;
    rhostart = 0.01;
    justprodlabor = true;
    res_cell = cell(1, length(yearsrhochange_vec) + 1);
    feedmu = nan;
    for (ii = 1:(length(res_cell)-1))
        yearsrhochange = yearsrhochange_vec(ii); 
        no_surprise = no_surprise_vec(ii); 

        [prod, quality, dL_t_dt, ...
        muder, xder, rdlabor, prodlabor, ...
        prodlaborder, markup_vec, ...
        rho_path, tplot, ximpact, ...
        muimpact, ximpactback, muimpactback, tplotback, ...
        quality2, ...
        prod2, prod3] = ...
                    transition_wrapper(m, rhostart, rhoend, years, yearsrhochange, justprodlabor, feedmu, no_surprise, ...
                                       (1 + rhotilde/100)^(1/12) - 1);

        res_cell{ii} = {prod, quality, dL_t_dt, ...
                            muder, xder, rdlabor, prodlabor, ...
                            prodlaborder, markup_vec, ...
                            rho_path, tplot, ximpact, ...
                            muimpact, ximpactback, muimpactback, tplotback, ...
                            quality2, ...
                            prod2, ...
                            prod3};
    end



    close all 
    figure;

        
    set(gcf, 'PaperUnits', 'inches');
    x_width=8;
    y_width=0.7*3.5;
    set(gcf, 'PaperPosition', [0 0 x_width y_width]); %

    subplot(1,3,1);
    p1 = plot([-5:-1, res_cell{1}{11}(1:(end-1))], ...
              [res_cell{1}{10}(1)*ones(1,5), res_cell{1}{10}(2:end)] - rhotilde,'LineStyle', '-', 'LineWidth', 2, ...
            'color', [0.4, 0.4, 0.4]);
    hold on ;
    p2 = plot([-5:-1, res_cell{2}{11}(1:(end-1))], ...
              [res_cell{2}{10}(1)*ones(1,5), res_cell{2}{10}(2:end)] - rhotilde, ...
              'LineStyle', '--', 'LineWidth', 2, 'color', [0, 0.5, 0]);
    ylabel('Percent', 'Interpreter', 'latex');
    xlabel('Years', 'Interpreter', 'latex');
    title({'Excess return, $\chi$'; ''}, 'Interpreter', 'latex');
    xlim([-5,40]);
    ylim([0.5,3.5]);

    subplot(1,3,2);
    p1 = plot([-5:-1, res_cell{1}{11}(1:(end-1))], ...
              100*(power(1+[res_cell{1}{1}(1)*ones(1,5), res_cell{1}{1}(2:end)],12)-1), 'LineStyle', ...
            '-', 'LineWidth', 2, 'color', [0.4, 0.4, 0.4]);
    hold on ;
    p2 = plot([-5:-1, res_cell{2}{11}(1:(end-1))], ...
             100*(power(1+[res_cell{2}{1}(1)*ones(1,5), res_cell{2}{1}(2:end)],12)-1), ...
             'LineStyle', '--', 'color', [0, 0.5, 0], 'LineWidth', 2);
    ylabel('Percent', 'Interpreter', 'latex');
    xlabel('Years', 'Interpreter', 'latex');
    title({'Productivity growth'; ''}, 'Interpreter', 'latex');
    xlim([-5,40]);

    subplot(1,3,3);
    p1 = plot([-5:-1, res_cell{1}{11}(1:(end-1))], ...
              (power(1+[res_cell{1}{9}(2,1)*ones(1,5), res_cell{1}{9}(2,2:end)],1)-1), 'LineStyle', ...
            '-', 'LineWidth', 2, 'color', [0.4, 0.4, 0.4]);
    hold on ;
    p2 = plot([-5:-1, res_cell{1}{11}(1:(end-1))], ...
              (power(1+[res_cell{1}{9}(3,1)*ones(1,5), res_cell{1}{9}(3,2:end)],1)-1), 'LineStyle', ...
            '--', 'LineWidth', 2, 'color', [0.4, 0.4, 0.4]);
    p3 = plot([-5:-1, res_cell{2}{11}(1:(end-1))], ...
              (power(1+[res_cell{2}{9}(2,1)*ones(1,5), res_cell{2}{9}(2,2:end)],1)-1), ...
              'LineStyle', '-', 'color', [0, 0.5, 0], 'LineWidth', 2);
    p4 = plot([-5:-1, res_cell{2}{11}(1:(end-1))], ...
              (power(1+[res_cell{2}{9}(3,1)*ones(1,5), res_cell{2}{9}(3,2:end)],1)-1), ...
              'LineStyle', '--', 'color', [0, 0.5, 0], 'LineWidth', 2);
    ylabel('Percent', 'Interpreter', 'latex');
    xlabel('Years', 'Interpreter', 'latex');
    title({'Net markup'; ''}, 'Interpreter', 'latex');
    xlim([-5,40]);
    text(10,38,'90th percentile', 'Interpreter', 'latex');
    text(10,14,'Median', 'Interpreter', 'latex');


    saveas(gcf, "Output/Figures_Paper/transition_1x3_v2.png");
    saveas(gcf, "Output/Figures_Paper/transition_1x3_v2.eps", 'epsc');




end