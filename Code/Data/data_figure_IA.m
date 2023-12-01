%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filename: data_figure_IA.m
% Author: Craig A. Chikis
% Date: 10/09/2023
% Note(s):
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function data_figure_IA() 
    
    plotdata = readtable("Output/Store_Data/yields.csv");
    plotdata = plotdata(year(plotdata.month) >= 1980, :);

    close all 
    figure; 


    set(gcf, 'PaperUnits', 'inches');
    x_width=8;
    y_width=3;
    set(gcf, 'PaperPosition', [0 0 x_width y_width]);


    subplot(1,3,1)
    p1 = plot(plotdata.month(plotdata.series == "One year"), plotdata.real_yield(plotdata.series == "One year"), ...
            'LineStyle', ':', 'LineWidth', 2, 'color', [0.75, 0.75, 0.75]);
    hold on 
    p2 = plot(plotdata.month(plotdata.series == "Five year"), plotdata.real_yield(plotdata.series == "Five year"), ...
            'LineStyle', '--', 'LineWidth', 2, 'color', [1, 0, 0]);
    p3 = plot(plotdata.month(plotdata.series == "Ten year"), plotdata.real_yield(plotdata.series == "Ten year"), ...
            'LineStyle', '-.', 'LineWidth', 2, 'color', [0, 0, 1]);
    ylabel('\%', 'Interpreter', 'latex')
    xlim([datetime(1980,1,1), datetime(2020,1,1)])
    xticks([datetime(1980,1,1), datetime(1990,1,1), datetime(2000,1,1), datetime(2010,1,1), datetime(2020,1,1)])
    title({'Risk-free rate'; ''}, 'Interpreter', 'latex')
    l1 = legend([p1, p2, p3], {'One year', 'Five year', 'Ten year'}, 'Interpreter', 'latex');
    set(l1, 'box', 'off')

    subplot(1,3,2) 
    p1 = plot(plotdata.month(plotdata.series2 == "Return on capital"), plotdata.real_yield(plotdata.series2 == "Return on capital"), ...
            '-k', 'LineWidth', 2); 
    ylim([0, 11])
    xlim([datetime(1980,1,1), datetime(2020,1,1)])
    xticks([datetime(1980,1,1), datetime(1990,1,1), datetime(2000,1,1), datetime(2010,1,1), datetime(2020,1,1)])
    ylabel('\%', 'Interpreter', 'latex')
    title({'Return on capital'; ''}, 'Interpreter', 'latex')


    subplot(1,3,3) 
    p1 = plot(plotdata.month(plotdata.series2 == "Excess return"), plotdata.real_yield(plotdata.series2 == "Excess return"), ...
            '-k', 'LineWidth', 2); 
    ylabel('\%', 'Interpreter', 'latex')
    xlim([datetime(1980,1,1), datetime(2020,1,1)])
    xticks([datetime(1980,1,1), datetime(1990,1,1), datetime(2000,1,1), datetime(2010,1,1), datetime(2020,1,1)])
    title({'Excess return'; ''}, 'Interpreter', 'latex')


    saveas(gcf, "Output/Figures_Paper/data_IA.eps", 'epsc')
    saveas(gcf, "Output/Figures_Paper/data_IA.png")

end
