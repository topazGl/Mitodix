function [ image, labeled, Z_dim ] = Read2D3D( im_path, seg_path, voxel_idx, rescaleFlag, downscaleFactor )
% Read intensity and labeled images.
%   [ image, labeled, Z_dim ] = Read2D3D( im_path, seg_path, voxel_idx, downscaleFactor )
% Inputs:
%   im_path = intensity image full path.
%   seg_path = segmentation image full path.
%   voxel_idx = sub image index for 3D (optional).
%   downscaleFactor = 0<downscaleFactor<=1 downscale factor (Default: 1, None).
    if(~exist('voxel_idx', 'var') || isempty(voxel_idx))
        voxel_idx = 0;
    end
    if(~exist('rescaleFlag', 'var') || isempty(rescaleFlag))
        rescaleFlag = 1;
    end
    if(~exist('downscaleFactor', 'var') || isempty(downscaleFactor))
        downscaleFactor = 1;
    end
    
    info = imfinfo(im_path);
    Z_dim = length(info);
    if(Z_dim == 1) % 2D
        image = imread(im_path);
        [labeled, status] = read_labeled(seg_path);
		[image, labeled] = preprocess(image, labeled, rescaleFlag, downscaleFactor);
    else % 3D
        if(not(voxel_idx))
            for z=1:Z_dim
                [labled_z, status] = read_labeled(seg_path, z);
                image_z = double(imread(im_path, 'Index', z));
				[image_z, labled_z] = preprocess(image_z, labled_z, rescaleFlag, downscaleFactor);
                if(z==1)
                    [Y, X] = size(image_z);
                    image = zeros(Y, X, Z_dim);
                    labeled = zeros(Y, X, Z_dim);
                end
                image(:, :, z) = image_z;
                labeled(:, :, z) = labled_z;
            end
        else
			Z_dim = 1;
			image = imread(im_path, 'Index', voxel_idx);
			[labeled, status] = read_labeled(seg_path, voxel_idx);
			[image, labeled] = preprocess(image, labeled, rescaleFlag, downscaleFactor);
        end        
    end
    if(~status)
        seg_path = fileparts(seg_path);
        if (~exist(seg_path, 'dir'))
            mkdir(seg_path);
        end
        [~, im_name, ext] = fileparts(im_path);
        imwrite(labeled, fullfile(seg_path, [im_name, ext]));
    end
end

function [labeled, status] = read_labeled(seg_path, voxel_idx)
    if(~exist(seg_path, 'file'))
        labeled = [];
        status = 0;
        return;
    end
    
    if(exist('voxel_idx', 'var') && ~isempty(voxel_idx))
        labeled = imread(seg_path, 'Index', voxel_idx);
    else
        labeled = imread(seg_path);
    end
    status = 1;
end

function [image, labeled, rescaleFlag] = preprocess(image, labeled, rescaleFlag, downscaleFactor)
	if(~exist('rescaleFlag', 'var') || isempty(rescaleFlag))
        rescaleFlag = 1;
    end
    
    if(~exist('downscaleFactor', 'var') || isempty(downscaleFactor))
        downscaleFactor = 1;
    end
    
    % Quantization:
    if(rescaleFlag && ~isa(image, 'uint8'))
        % Work in uint8:
        if(isa(image, 'uint8'))
            rescaleFlag = 8;
        else
%                 warning('8 bit input expected');
            if (isa(image, 'uint16'))
                rescaleFlag = 16;
            elseif(isa(image, 'uint32'))
                rescaleFlag = 32;
            else %uint64
                rescaleFlag = 64;
            end
        end
    end

    image = double(image);
    if(rescaleFlag)
        % map to 0-1:
        image = image/(2^rescaleFlag-1);
    end
    
    % clean salt-pepper noise:
    noisy_image = image;
    cImage = medfilt2(noisy_image);
    hImage = max(cImage, noisy_image);
    
    % segmentation based on intensities:
    if(isempty(labeled))
        dc = mean(cImage(:));
        fg_pixels = cImage(cImage(:)>dc);
        thr = AutoThr( fg_pixels, '', 0, 4);
%         thr = max(thr, 0.0015);
%         thr = min(thr, 0.0018);
        labeled = (cImage > thr);
        labeled = medfilt2(labeled);
        labeled = SplitBigClusters(labeled);
    end

    % Set BG label to 0:
    labeled = labeled - min(reshape(labeled,[],1));
    lblVec = reshape(labeled,[],1);
    if(sum(lblVec ~= 0) > sum(lblVec == 0))
        labeled = not(labeled);
    end
            
	% Fill holes:
    labeled = imfill(labeled,'holes');
    
    image(labeled == 0) = cImage(labeled == 0);
    image(labeled ~= 0) = hImage(labeled ~= 0);
    diff = abs(noisy_image - image);
    idx = diff > 0.01*max(reshape(noisy_image,[],1));
    image(idx) = noisy_image(idx);
    
    if(rescaleFlag)
        % Set BG to zero:
        bgPixels = image(labeled(:) == 0);
        image = image - mean(bgPixels);
%         image(image < 0) = 0;
        bg_val =  min(0, min(image(labeled(:) > 0)));
        if(isempty(bg_val))
            bg_val = 0;
        end
        image((labeled == 0) & (image<std(bgPixels))) = bg_val;
    end
            
	% downscale:	
	if(abs(downscaleFactor) > 1)
        downScaleLabeledFlag = (downscaleFactor>0);
        downscaleFactor = abs(downscaleFactor);
        if(~downScaleLabeledFlag)
            image = imgaussfilt(image,downscaleFactor);
%             image = imresize(image, 1/downscaleFactor, 'bicubic');
%             image = imresize(image, size(labeled), 'bicubic');
        else
            image = imresize(image, 1/downscaleFactor, 'bicubic');
            if(max(labeled(:)==1))
                CC = bwconncomp(labeled~=0);
                labeled = labelmatrix(CC);
            end
            labeled = imresize(labeled, 1/downscaleFactor, 'nearest');
            %labeled = uint16(labeled);
        end
    end
end