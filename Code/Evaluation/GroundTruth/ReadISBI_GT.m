        function [ mitosisGT ] = ReadISBI_GT(InputPath, OutputPath)
        % Each line: L B E P
        % Where L = cell label, B = birth frame, E = end frame,
        % P = parent label.
            mitosis = [];
            frames = [];
            file_name = fullfile(InputPath, 'man_track.txt');
            if(~exist(file_name, 'file'));
                file_name = fullfile(InputPath, 'res_track_4c.txt');
                if(~exist(file_name, 'file'))
                   file_name = fullfile(InputPath, 'res_track.txt');
                end
            end
            fileID = fopen(file_name,'r');
            
            size_data = [4 Inf];
            data = fscanf(fileID,'%d %d %d %d\n', size_data);

            fclose(fileID);
            N = size(data, 2);
            candidate_index = 0;
            for i=1:(N-1)
                if (~data(4, i))
                    continue;
                end
                second = 0;
                for j=i+1:N
                    % same parent and same birth frame:
                    if ((data(4, i) == data(4, j)) && (data(2, i) == data(2, j)))
                        second = j;
                        break;
                    end
                end
                if(~second)
                    continue;
                end
                
                candidate_index = candidate_index + 1;
                frameIndex = data(2, i);
                frames = [frames, frameIndex + 1];
                
                % Daughter1:
                cell.id = data(1, i);
                cell.frameIndex = frameIndex + 1;
                mitosis{candidate_index,1} = cell;


                
                % Daughter2:
                cell.id = data(1, second);
                cell.frameIndex = frameIndex + 1;
                mitosis{candidate_index,2} = cell;
                
                % Mother:
                cell.id = data(4, second);
                cell.frameIndex = frameIndex;
                mitosis{candidate_index,3} = cell;
            end

            max_frame = max(frames);
            min_index = 0;
            %sort:
            for i=1:candidate_index
                min_val = min(frames);
                if(min_val > max_frame)
                    break;
                end

                min_index = find(frames == min_val);
                if length(min_index) > 1
                    min_index = min_index(1);
                end

                %copy values:
                mitosisGT{i,1} = mitosis{min_index,1};
                mitosisGT{i,2} = mitosis{min_index,2};
                mitosisGT{i,3} = mitosis{min_index,3};
                
                frames(min_index) = max_frame+1;
            end
            
            if(~exist(OutputPath, 'dir'))
                mkdir(OutputPath);
            end
            save(fullfile(OutputPath, 'mitosisGT.mat'), 'mitosisGT');
        end

