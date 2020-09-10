function [ t1f, V1f ] = filterMaster( V1, t, USE_LAB_FILTER, SAMPLE_RATE )

    % Filtering parameters
    %BREAK = 500;
    Fline = 60; % power line frequency
    if (USE_LAB_FILTER ~= 1)
        %disp(size(V1'))
        sbar = removePLI_multichan([V1',V1'],SAMPLE_RATE,7,[70,0.01,4],[0.01,4,4],1,Fline,1);
        % Good for first data set:
        %sbar = removePLI_multichan(V1',fs,4,[70,0.01,4],[0.09,4,4],1,60,1); 
        %makeSpectralPlot(t(BREAK:end), sbar(BREAK:end), e1, '_filtered_Keshtkaran')
        %fprintf('[!]chck1[!]\n')
        V1f = sbar((length(V1)+1):end)';
        %V1f = sbar(1:length(V1))';
        %V1f = sbar;
        %fprintf('[!]chck2[!]\n')
        t1f = t(1:end)';
        %fprintf('[!]chck3[!]\n')
        %disp(size(t1f))
        %disp(size(V1f))
        %disp(size(V1))
        %fprintf('EXISTFROMFILTERMASTER\n')
    else
        %BREAK = 1;
        %Vpad = normrnd(mean(V1),sqrt(var(V1)),10*SAMPLE_RATE,1); % 10 seconds worth of padding V1(1)
        Vpad = V1;
        %Vpad = [];
        % Apply lab filter
        % notch filter at 60 Hz and harmonics
        V1s = [Vpad;V1];
        for f=Fline:Fline:(SAMPLE_RATE/2+Fline)
            V1s = FilterData(V1s,SAMPLE_RATE,'notch',f);
        end
        %makeSpectralPlot(t(BREAK:end), V1s(BREAK:end), e1, '_filtered_labFilter')
        V1f = V1s((length(Vpad)+1):end);
        %t1f = t((length(Vpad)+1):end)';
        t1f = t;
    end
    %fprintf('SLIM SHADY\n')
    %disp(size(V1f));
    %disp(size(t1f));
end

