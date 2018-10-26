function mitoticFeaturesSpaceAnalysis( mitoticFeatures, gm, c, outputPath, filtered_GT, candidates, radi, timeDiff, features_names, d )
% Analyze mitotic features space. 
% Signature:
%   mitoticFeaturesSpaceAnalysis( mitoticFeatures, gm, c, outputPath, filtered_GT, candidates, radi, timeDiff )
%   mitoticFeaturesSpaceAnalysis( mitoticFeatures, gm, c, outputPath, filtered_GT, candidates )
% Inputs:
%   mitoticFeatures: Nx2 features matrix.
%   plotFlag: 0/1 flag (Default=0).
%   outputPath: full path for saving outputs (relevant for plotFlag=1).
%   allTPs: Nx1 binary vector where 1 indicate true positives (Used for visualization only, plotFlag=1).
%%
    if(isempty(gm) || isempty(mitoticFeatures))
        return;
    end
    
    if(~exist('d', 'var') || isempty(d))
        d = min(3, size(mitoticFeatures, 2));
    end
    
    circSize = 30;

    outputPath = fullfile(outputPath, 'cluster');
    if(~exist(outputPath, 'dir'))
        mkdir(outputPath);
    end
    plotPath = fullfile(outputPath, 'plots');
    if(~exist(plotPath, 'dir'))
        mkdir(plotPath);
    end

    detections = checkDetections(filtered_GT, candidates, -mitoticFeatures(:,2), radi, timeDiff);
    allTPs = detections(detections ~= 0);
    
    group = zeros(size(mitoticFeatures, 1), 1);
    group(allTPs) = 1;
    
    if(d==3)
        f1 = figure;
        set(f1, 'Position', get(0, 'Screensize'));
        
        h = scatter3(mitoticFeatures(group==0, 1), mitoticFeatures(group==0, 2), mitoticFeatures(group==0, 3), 8, '.b'); hold on;
        scatter3(mitoticFeatures(group==1, 1), mitoticFeatures(group==1, 2), mitoticFeatures(group==1, 3), circSize, '.r');
        xlim([eps 1-eps]); ylim([eps 1-eps]); zlim([eps 1-eps]);
        xlabel(features_names{1});
        ylabel(features_names{2});
        zlabel('feature3');

        set(gca, 'FontSize', 36);
        l = legend(h, {'FP', 'TP'});
        set(l, 'FontSize', 36, 'Location', 'Best');
        saveas(f1, fullfile(plotPath, 'scatter3.fig'));
        saveas(f1, fullfile(plotPath, 'scatter3.jpg'));
        close(f1);
    end
    
    f1 = figure;
    set(f1, 'Position', get(0, 'Screensize'));
    cmap = flipud(colormap('gray'));
    
    [X,Y,Z] = plotPDF( gm, [0 1], 1000, c ); hold on;
    %mesh(X, Y, -Z);view(2);colormap(cmap);
    
    h = gscatter(mitoticFeatures(:, 1), mitoticFeatures(:, 2), group, ['b','r'], '.', [8, circSize]); hold on;

    xlim([min(mitoticFeatures(:,1)) max(mitoticFeatures(:,1))]); ylim([min(mitoticFeatures(:,2)) max(mitoticFeatures(:,2))]);
    xlabel(features_names{1});
    ylabel(features_names{2});
    
    set(gca, 'FontSize', 36);
    l = legend(h, {'FP', 'TP'});
    set(l, 'FontSize', 36, 'Location', 'Best');
    saveas(f1, fullfile(plotPath, 'scatter.fig'));
    saveas(f1, fullfile(plotPath, 'scatter.jpg'));

%     contour(X,Y,Z,[0.5, 0.5],'LineWidth',5);
    contour(X,Y,Z,'LineWidth',2); 
    colormap(cmap); colorbar; hold off;
    xlim([min(mitoticFeatures(:,1)) max(mitoticFeatures(:,1))]); ylim([min(mitoticFeatures(:,2)) max(mitoticFeatures(:,2))]);
    saveas(f1, fullfile(plotPath, 'scatter_cluster.fig'));
    saveas(f1, fullfile(plotPath, 'scatter_cluster.jpg'));
    
    [detection_labels, labels_connections] = labelTriplets(filtered_GT, candidates, radi, timeDiff );
    f2 = figure;
    set(f2, 'Position', get(0, 'Screensize'));
    circSize = 30;
    cmap = [
            0   0   1
            0.2081 0.1663 0.5292
            0.0282 0.6663 0.7574
            0.4783 0.7489 0.4877
            0.9763 0.9831 0.0538
            1   0   0
            ];
    C = zeros(size(mitoticFeatures, 1), 3);
    for cc = 1:size(mitoticFeatures, 1)
        C(cc, :) = cmap(detection_labels(cc)+1, :);
    end
    legend_handles=[];
    h=scatter(mitoticFeatures(:,1), mitoticFeatures(:,2), circSize, cmap(1,:), '.', 'LineWidth', 5, 'DisplayName', 'non valid'); hold on;
    if(any(detection_labels==0))
        legend_handles = [legend_handles, h];
    end
%     h=scatter(mitoticFeatures(detection_labels==1,1), mitoticFeatures(detection_labels==1,2), circSize, '*', 'MarkerFaceColor', cmap(2,:), 'LineWidth', 5, 'DisplayName', 'one daughter valid');
%     if(any(detection_labels==1))
%         legend_handles = [legend_handles, h];
%     end
%     h=scatter(mitoticFeatures(detection_labels==2,1), mitoticFeatures(detection_labels==2,2), circSize, 'x', 'MarkerFaceColor',cmap(3,:), 'LineWidth', 5, 'DisplayName', 'only mother valid');
%     if(any(detection_labels==2))
%         legend_handles = [legend_handles, h];
%     end
    h=scatter(mitoticFeatures(detection_labels==3,1), mitoticFeatures(detection_labels==3,2), 2*circSize, cmap(4,:), 'd', 'filled');
    set(h, 'LineWidth', 5, 'DisplayName', 'one daughter invalid');
    if(any(detection_labels==3))
        legend_handles = [legend_handles, h];
    end
    h=scatter(mitoticFeatures(detection_labels==4,1), mitoticFeatures(detection_labels==4,2), 2*circSize, cmap(5,:), 's', 'filled');
    set(h, 'MarkerEdgeColor',[0 0 0], 'DisplayName', 'only daughters valid');%
    if(any(detection_labels==4))
        legend_handles = [legend_handles, h];
    end
    h=scatter(mitoticFeatures(detection_labels==5,1), mitoticFeatures(detection_labels==5,2), circSize, cmap(6,:), '*');
    set(h, 'LineWidth', 5, 'DisplayName', 'all valid');
    if(any(detection_labels==5))
        legend_handles = [legend_handles, h];
    end
    xlim([min(mitoticFeatures(:,1)) max(mitoticFeatures(:,1))]); ylim([min(mitoticFeatures(:,2)) max(mitoticFeatures(:,2))]);
    xlabel('mother-daughters similarity');
    ylabel('daughters similarity');
    set(gca, 'FontSize', 36);
    l = legend(legend_handles);
    set(l, 'FontSize', 30, 'Location', 'Best');
    hold off;
    saveas(f2, fullfile(plotPath, 'scatter_labels.fig'));
    saveas(f2, fullfile(plotPath, 'scatter_labels.jpg'));

    [detection_labels_dt, labels_connections_dt] = labelTripletsTimeGap(filtered_GT, candidates, radi );
    f3 = figure;
    set(f3, 'Position', get(0, 'Screensize'));
    
    legend_handles=[];
    h=scatter(mitoticFeatures(:,1), mitoticFeatures(:,2), circSize, cmap(1,:), '.', 'LineWidth', 5, 'DisplayName', 'non valid'); hold on;
    if(any(detection_labels_dt==-1))
        legend_handles = [legend_handles, h];
    end
    
    for dt = 4:-1:0
        h=scatter(mitoticFeatures(detection_labels_dt==dt,1), mitoticFeatures(detection_labels_dt==dt,2), 2*circSize, cmap(end-dt,:), 'd', 'filled');
        set(h, 'LineWidth', 5, 'DisplayName', ['\Delta{t}=', num2str(dt)]);
        if(any(detection_labels_dt==dt))
            legend_handles = [legend_handles, h];
        end
    end
    xlim([min(mitoticFeatures(:,1)) max(mitoticFeatures(:,1))]); ylim([min(mitoticFeatures(:,2)) max(mitoticFeatures(:,2))]);
    xlabel(features_names{1});
    ylabel(features_names{2});
    set(gca, 'FontSize', 36);
    l = legend(legend_handles);
    set(l, 'FontSize', 30, 'Location', 'Best');
    hold off;
    saveas(f3, fullfile(plotPath, 'scatter_DT_labels.fig'));
    saveas(f3, fullfile(plotPath, 'scatter_DT_labels.jpg'));
    
    close(f1); close(f2); close(f3);
    save(fullfile(outputPath, 'clustering.mat'), 'mitoticFeatures', 'detection_labels', 'labels_connections', 'detection_labels_dt', 'labels_connections_dt', 'gm', 'c', '-v7.3');
end

function [X, Y, Z] = plotPDF( gm, range, n, c )
    dt = (range(end) - range(1))/n;
    x = range(1):dt:range(end);
    
    [X, Y] = meshgrid(x, x);
%     post = posterior(gm, [X(:), Y(:)]);
%     post = post(:, c);
%     Z = reshape(post, size(X));
    Z = reshape(ExponentialDistance([X(:), Y(:)], gm, c), size(X));
    
%     C = size(gm.mu, 1);
%     dists = cell(C,1);
%     Z = 0;
%     for i=1:C
%         MU = gm.mu(i,:);
%         Sigma = gm.Sigma(:,:,i);
%         pdf_val = mvnpdf([X(:), Y(:)], MU, Sigma);
%         dists{i} = reshape(pdf_val, size(X));
%         Z = Z + gm.ComponentProportion(i)*dists{i};
%     end
%     
%     if(exist('c', 'var') && ~isempty(c))
%         Z = dists{c};
%     end
end

% function dist = exponential_distance(x, gm, c)
%     P = gm.ComponentProportion(c);
%     MU = gm.mu(c,:);
%     Sigma = gm.Sigma(:,:,c);
%     dist = (1/P)*sqrt(det(Sigma))*exp(0.5*(x-MU)'*inv(Sigma)*(x-MU));
% end