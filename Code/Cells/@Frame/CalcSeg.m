function seg = CalcSeg(image, labeled)
% Re-calculate image segmentation
%   seg = CalcSeg(image, labeled)
% Inputs:
%   image = gray-level image.
%   labeled = image initial segmentation.
% Ouputs:
%   seg = output segmentation image.
    try
        close all;
        addpath(fullfile('Algorithms','FastMarching'));
        addpath(fullfile('Algorithms','FastMarching', 'functions'));
        [h,w] = size(image);
        s = regionprops(labeled, 'Centroid', 'MajorAxisLength');
        lengths = cat(1, s.MajorAxisLength);
        labels = find(lengths > 0);
        COMs = floor(cat(1, s(labels).Centroid));

        regNo = 2;
        image = double(image);
        GMMParams = fitgmdist(image(:),regNo,'RegularizationValue', 0.0005);%,'Options',statset('Display','final')); % 2 soft clusters: objects (cells) & background

        % calc Image Liklihood
        P = posterior(GMMParams,image(:));
        softLabel = zeros(h,w,regNo);
        for L = 1:regNo
           softLabel(:,:,L) = reshape(P(:,L),[h,w]);
           %figure;imagesc(softLabel(:,:,L));colorbar
        end

        %figure;imagesc(softLabel(:,:,1));colorbar
        %figure;imagesc(softLabel(:,:,2));colorbar

        [~,hardLabel] = max(softLabel,[],3);
        if (sum(hardLabel(:) ==1)) > (sum(hardLabel(:) ==2)),
             bckIdx =1; cellIdx =2;
        else
             bckIdx =2; cellIdx =1;
        end

        gradMap =  abs(gradient(image/max(image(:))))+0.001;

        ptsNo = size(COMs,1);
        PostProb= zeros(h,w,ptsNo+1);
        probGDist = zeros(h,w,ptsNo+1);

        for k=1:ptsNo %msfm2d -- mex file; Matlab Central: fast marching 
            GDist= msfm2d(1./(gradMap).^(0.2), COMs(k,:)', false, false); 
            %figure;imagesc(GDist);colorbar; title(['GDist-',int2str(k)]);
            probGDist(:,:,k) = 1./(1+ GDist);
            %figure;imagesc(probGDist(:,:,k));colorbar;title(['Prob-GDist-',int2str(k)]);
            PostProb(:,:,k)= probGDist(:,:,k).*softLabel(:,:,cellIdx);
            %figure;imagesc(PostProb(:,:,k));colorbar;title(['Prob-Prob-',int2str(k)]);
        end
        probGDist(:,:,ptsNo+1) = 1./(ptsNo+1).*ones(size(image));% for demonstration only
        [~,seg] = max(probGDist,[],3);

        %figure;imagesc(seg);colorbar % for demonstration only


        PostProb(:,:,ptsNo+1) = softLabel(:,:,bckIdx);

        [~,seg] = max(PostProb,[],3);

        % Keep the original label numbers:
        new_labeled = labeled;
        for n=1:length(labels)
            old_blob_idx = (labeled == label(n));
            overlap_idx = seg(old_blob_idx);
            if(any(overlap_idx))
                new_labeled(old_blob_idx) = 0;
                new_labeled(overlap_idx) = label(n);
            end
        end

        %figure;imagesc(seg);colorbar
        figure;imshow(image/max(image(:)));hold
        colArr = {'r','m','k','y','g','b','c'}; 
        for j=ptsNo:-1:1,
           contour(seg,[j+.1,j+.1],colArr{rem(j,7)+1});
        end
    catch EX
        throw(EX);
    end
end