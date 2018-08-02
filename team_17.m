function [] = dsaaproject(trainDir, fileDir, outputDir)
    
    %disp(char(trainDir));
    trainDir = char(trainDir);
    %disp(char(fileDir));
    fileDir = char(fileDir);
    %disp(char(outputDir));
    outputDir = char(outputDir);

    dirName = strcat(char(fileDir), '*.mat');
    files = dir(dirName);
    for file = files'
        % start of the directory's for loop

        clearvars -except file files dirName trainDir fileDir outputDir sig;

        filename = strcat(char(fileDir), char(file.name));

        % Loaded data sig
        load(filename);
    
        % Fixing parameters
        samplingRate    = 125;
        bandPassRange   = [24 240];
        [bw1 bw2]       = butter(4, bandPassRange/(samplingRate*60/2), 'bandpass');
        newSamplingRate = 25;
        factor          = samplingRate/newSamplingRate;
    
        iterations = floor((length(sig) - 8*samplingRate)/(2*samplingRate)) + 1;
        %disp(iterations);
        pred       = zeros(1, 125);
        FFTMat     = zeros(4, 1024);
        curRange   = [];
        iterations = min(iterations, 125);
    
        for i = (1 : iterations)     
            %% Preprocessing
        
            curRan = (i-1)*2*samplingRate + 1 : ((i-1)*2*samplingRate + 8*samplingRate);
            %disp(size(curRange));
            initRangeData = sig(:, curRan);
            %disp(i);
        
            for j = (1 : 5)
                initRangeData(j, :) = filter(bw1, bw2, initRangeData(j, :));
            end
        
            initRangeData(2, :) = (zscore(initRangeData(1, :), 0, 2) + zscore(initRangeData(2, :), 0, 2)) / 2;
        
            inRangeData = zeros(5, length(initRangeData(1,:))/5);
        
            for j = (2 : 5)
                inRangeData(j, :) = downsample(initRangeData(j, :), factor);
                FFTMat(j-1, :) = fft(inRangeData(j, :), 1024);
            end
        
            % explain and might be removable
            mapFreq = linspace(0, newSamplingRate, 1024);
            [~, lwBound] = min(abs(mapFreq - 1));
            [~, upBound] = min(abs(mapFreq - 3));
            mapFreq = mapFreq(lwBound : upBound);
            
            FFTMatRan = zeros(4, length(mapFreq));
        
            for j = (1 : 4)
                FFTMatRan(j, :) = FFTMat(j, lwBound : upBound);
            end
        
            %% Filtering
        
            len = 15;

            wnrFFT(i, :) = normalize(FFTMatRan(1, :));
            normFFTMat = zeros(size(FFTMatRan));
        
            if i == 1
                normFFTMat(1, :) = normalize(wnrFFT(i,:));
            else
                normFFTMat(1, :) = normalize(mean(wnrFFT(max(1, i-len):i, :), 1));
            end
        
            % normFFTMat(1, :) = normalize(normFFTMat(1, :));
        
            for j = (2 : 4)
                normFFTMat(j, :) = normalize(FFTMatRan(j, :));
            end
        
            Fin1 = 1 - 1/3 * sum(normFFTMat(2:4, :))./(normFFTMat(1, :));
            Fin1 (Fin1 < 0) = -1; % limiting value between -inf and -1
            %disp(sum(Fin1));
            Denoised1 = abs(FFTMatRan(1, :)).*Fin1;
            %Denoised1 = Denoised1/std(Denoised1);

            Fin2 = normFFTMat(1, :)./(normFFTMat(1, :) + sum(normFFTMat(2:4, :))./3);
            %disp(sum(Fin));
            Denoised2 = abs(FFTMatRan(1, :)).*Fin2;
            %Denoised2 = Denoised2/std(Denoised2);

            Denoised = Denoised1 + Denoised2;
        
            %% Post-processing
        
            beatRange = 25;
            if i > 15
                beatRange = max(abs(diff(pred(1:i-1)))) + 5;
            end
            
            [~, ind] = max(Denoised(1,:));
            rangeVal = ind(1);

            if i>1
                [~, ind] = max(Denoised(1, curRange));
                rangeVal = curRange(ind(1));
            end

            pred(i) = mapFreq(rangeVal) * 60;
            curRange = rangeVal - round(beatRange/((mapFreq(2)-mapFreq(1))*60)) : rangeVal + round(beatRange/((mapFreq(2)-mapFreq(1))*60));
        
            curRange(curRange<1) = []; curRange(curRange>length(mapFreq)) = [];
        
        end

        opFileName = strcat(outputDir, strcat('output_team_17_', char(file.name) ) );
        %opFileName = strcat('output_team_17_', filename(end-17:end));
        %disp(pred);
        %disp(opFileName);
        pred = pred.';
        save(opFileName, 'pred');
    end % end of the directory for loop'


end

function [output] = normalize(input)
    output = abs(input)/max(abs(input));
end
