function SaveImages(frames, outputPath, com, R)
% Save frames patch images.
%   SaveImages(frames, outputPath, com, R)
% Inputs:
%   frames = array of Frame objects.
%   outputPath = output folder path.
%   COM = patch center.
%   R = roi around patch.
    N = length(frames);
    if(~N)
        return;
    end
    com = round(com);
    h = figure(); colormap(h,gray); 
    for n=1:N
        frame = frames(n);
        if(n==1)
            if(nargin < 4)
                R = ceil(0.5*frame.DistBetweenNeigbors.Mean);
            else
                R = ceil(R);
            end
            start = com - [R, R];
            start = max(start, [1, 1]);
            finish = com + [R, R];
            finish = min(finish, [frame.Size(2), frame.Size(1)]);
            p = round(com - start);
        end
        a = subplot(1,N,n, 'Parent', h);
        imagesc(frame.Image(start(2):finish(2), start(1):finish(1)), 'Parent', a);
        hold on; plot(p(1), p(2), 'r.', 'Parent', a); hold off;
        title(['Frame ', num2str(frame.Number)], 'Parent', a);
    end
    saveas(h,[outputPath,'_im.jpg']);
    close(h);

    h = figure(); colormap(h,jet); 
    for n=1:N
        frame = frames(n);
        a = subplot(1,N,n, 'Parent', h);
        imagesc(frame.Labeled(start(2):finish(2), start(1):finish(1)), 'Parent', a);%subimage(frame.Labeled(start(1):finish(1), start(2):finish(2)));
        %hold on; plot(p(1), p(2), 'r.', 'Parent', a); hold off;
        title(['Frame ', num2str(frame.Number)], 'Parent', a);
    end
    saveas(h,[outputPath,'_seg.jpg']);
    close(h);
end

