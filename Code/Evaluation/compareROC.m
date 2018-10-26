function [f, g, pval]  = compareROC( FPR, TPR, precision, names, outputPath, f, g, color )
    outputPlots = fullfile(outputPath, 'Compare');
    if(~exist(outputPlots, 'dir'))
        mkdir(outputPlots);
    end
    
    setColor = (exist('color', 'var') && ~isempty(color));
    
    AUC = zeros(2,1);
    maxF1 = zeros(2,1);
    F1 = cell(2,1);
    for ii = 1:2
        AUC(ii) = trapz(FPR{ii}, TPR{ii});
        F1{ii} = 2*(TPR{ii}.*precision{ii})./(TPR{ii} + precision{ii});
        maxF1(ii) = max(F1{ii});
    end
    
    if(~exist('f', 'var') || isempty(f))
        f = figure;
    else
        set(f, 'Position', get(0, 'Screensize'));
        h = findobj(f, 'Type', 'Axes');
        axes(h);
        hold on;
    end
    
    if(~setColor)
        plot(TPR{1}, precision{1}, 'LineWidth', 3);
    else
        plot(TPR{1}, precision{1}, ['--', color], 'LineWidth', 3);
    end
    
    xlabel('TPR');
    ylabel('Precision');
    hold on;
    
    if(~setColor)
        plot(TPR{2}, precision{2}, 'LineWidth', 3);
    else
        plot(TPR{2}, precision{2}, color, 'LineWidth', 3);
    end

    l = findobj(gcf, 'Type', 'Legend');
    newEnt = {[names{1}, ', max F1=', num2str(maxF1(1), 3)], [names{2}, ', max F1=', num2str(maxF1(2), 3)]};
    if(isempty(l))
        l = legend(newEnt);
        set(l, 'FontSize', 20);
        set(gca, 'FontSize', 24);
    else
        str = get(l, 'String');
        str{end+1} = newEnt{1};
        str{end+1} = newEnt{end};
        legend('off');
        l = legend(str);
        set(l, 'FontSize', 28);
        legend('show');
        set(gca, 'FontSize', 34);
    end
    set(l, 'Location', 'Best');
    
    xlim([0,1]); ylim([0,1]);
    saveas(f, fullfile(outputPlots, 'comparePR.jpg'));
    saveas(f, fullfile(outputPlots, 'comparePR.fig'));
    
    if(~exist('g', 'var') || isempty(g))
        g = figure;
    else
        set(g, 'Position', get(0, 'Screensize'));
        h = findobj(g, 'Type', 'Axes');
        axes(h);
        hold on;
    end
    
    if(~setColor)
        plot(FPR{1}, TPR{1}, 'LineWidth', 3);
    else
        plot(FPR{1}, TPR{1}, ['--', color], 'LineWidth', 3);
    end
    hold on;
    
    if(~setColor)
        plot(FPR{2}, TPR{2}, 'LineWidth', 3);
    else
        plot(FPR{2}, TPR{2}, color, 'LineWidth', 3);
    end
    hold off;
    xlabel('FPR'); ylabel('TPR');
    
    l = findobj(gcf, 'Type', 'Legend');
    newEnt = {[names{1}, ', AUC=', num2str(AUC(1), 3)], [names{2}, ', AUC=', num2str(AUC(2), 3)]};
    if(isempty(l))
        l = legend(newEnt);
        set(l, 'FontSize', 20);
        set(gca, 'FontSize', 24);
    else
        str = get(l, 'String');
        str{end+1} = newEnt{1};
        str{end+1} = newEnt{end};
        legend('off');
        l = legend(str);
        set(l, 'FontSize', 30);
        legend('show');
        set(gca, 'FontSize', 32);
    end
    set(l, 'Location', 'Best');
    
	xlim([0,0.5]); ylim([0,1]);
    saveas(g, fullfile(outputPlots, 'compareROCzoom.jpg'));
    saveas(g, fullfile(outputPlots, 'compareROCzoom.fig'));
	
    xlim([0,1]); ylim([0,1]);
    saveas(g, fullfile(outputPlots, 'compareROC.jpg'));
    saveas(g, fullfile(outputPlots, 'compareROC.fig'));
	
    [~, pval] = ttest2(F1{1}, F1{2}, 'Vartype', 'unequal');
    save(fullfile(outputPlots, 'compare.mat'), 'FPR', 'TPR', 'precision', 'F1', 'names', 'AUC', 'maxF1', 'pval', '-v7.3');
end

% load('E:\Topaz\Google Drive\Academy\Thesis\Code\Data\Outputs\GalitLahav\Jose\Jose3\M20150416_p21_20X\s102\Results_2018-02-03_22-22-51-561\Compare\compare.mat');
% names = {'SVD RPE1'; 'Ours stg1 RPE1'};
% [f, g, pval]  = compareROC( FPR, TPR, precision, names, '..\Data\Outputs', [], [], 'm' );
% load('E:\Topaz\Google Drive\Academy\Thesis\Code\Data\Outputs\Broad\Results_2018-02-17_16-11-08-766\Compare\compare.mat')
% names = {'SVD MCF-10A'; 'Ours stg1 MCF-10A'};
% [f, g]  = compareROC( FPR, TPR, precision, names, '..\Data\Outputs', f, g, 'r' );
% load('E:\Topaz\Google Drive\Academy\Thesis\Code\Data\Outputs\ISBI\Fluo-N2DL-HeLa\01\Results_2018-02-04_21-28-03-757\Compare\compare.mat')
% names = {'SVD N2DL-HeLa01'; 'Ours stg1 N2DL-HeLa01'};
% [f, g]  = compareROC( FPR, TPR, precision, names, '..\Data\Outputs', f, g, 'b' );
% load('E:\Topaz\Google Drive\Academy\Thesis\Code\Data\Outputs\ISBI\Fluo-N2DL-HeLa\02\Results_2018-02-03_09-25-32-522\Compare\compare.mat')
% names = {'SVD N2DL-HeLa02'; 'Ours stg1 N2DL-HeLa02'};
% [f, g]  = compareROC( FPR, TPR, precision, names, '..\Data\Outputs', f, g, 'g' );