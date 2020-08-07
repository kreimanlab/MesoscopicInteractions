close all;
clear;

%fn_o = 'opti_broadband2_p';
fn_o = 'out/oo6';

fo = fopen(fn_o,'r');
Subs = {};
Metrics = {};
Coh_t = [];
Ct_t = [];
Cp_t = [];
F_x = [];
F_x_n = [];
F_y = [];
F_y_n = [];
count = 1;
while ~feof(fo)
    line = fgets(fo); %# read line by line
    line = strsplit(line);
    
    if (~strcmp(line{1},'Est.'))
    
        % [TEMP] Trim first element
        line = line(2:end);

        Metrics{count} = line{2};
        if (startsWith(line{1},'m00'))
            Subs{count} = line{1};
            Coh_t(count) = str2double(line{3});
            Ct_t(count) = str2double(line{4});
            F_x(count) = str2double(line{5});
            F_x_n(count) = str2double(line{6});
        elseif (startsWith(line{1},'xsub'))
            Cp_t(count) = str2double(line{3});
            F_y(count) = str2double(line{4});
            F_y_n(count) = str2double(line{5});
        end
        %fprintf(line);
        count = count + 1;
    end
end
fclose(fo);

% final variables
CohT = [];
CtT = [];
CpT = [];
Fx = {};
Fxn = {};
Fy = [];
Fyn = [];

% init
f_x = [];
f_x_n = [];
count = 1;
Subs = {Subs{:},[]};
for i = 1:length(Subs)
    if (isempty(Subs{i}))
        CohT = [CohT; coh_t];
        CtT = [CtT; ct_t];
        CpT = [CpT; Cp_t(i)];
        Fx = {Fx{:}, sprintf('[%.3f]%.3f[%.3f]',min(f_x),mean(f_x),max(f_x)) };
        Fxn = {Fxn{:}, sprintf('[%.3f]%.3f[%.3f]',min(f_x_n),mean(f_x_n),max(f_x_n)) };
        Fy = [Fy; F_y(i)];
        Fyn = [Fyn; F_y_n(i)];
        %return
        % reset
        f_x = [];
        f_x_n = [];
        count = 1;
    else
        f_x(count) = F_x(i);
        f_x_n(count) = F_x_n(i);
        coh_t = Coh_t(i);
        ct_t = Ct_t(i);
    end
    
    count = count + 1;
end

% Print output
UNK = '[-----]-----[-----]';
uCohT = unique(CohT);
uCtT = unique(CtT);
uCpT = unique(CpT);
fprintf('%s\n\n',Metrics{1});
for i1 = 1:length(uCohT)
    % Table head
    head = sprintf('Coh=%.3f\t|\t',uCohT(i1));
    for i2 = 1:length(uCtT)
        head = [head,sprintf('CT:%.3f\t\t\t',uCtT(i2))];
    end
    fprintf('%s\n',head);
    fprintf('---------------------------------------------------------------------------------------------------------------------------------------------------------\n')
    for i3 = 1:length(uCpT)
        
        % First line
        line = sprintf('CP:%.3f\t|\t',uCpT(i3));
        for i2 = 1:length(uCtT)
            % Find index
            idx = find( ( CpT == uCpT(i3) ) & ( CtT == uCtT(i2) & ( CohT == uCohT(i1) )) ,1);
            if (~isempty(idx))
                addon = sprintf('%s/%.3f\t',Fx{idx},Fy(idx));
            else
                addon = sprintf('%s/-----\t',UNK);
            end
            line = [line,addon];
        end
        fprintf('%s\n',line)
        
        % Second line
        line = sprintf('        \t|\t');
        for i2 = 1:length(uCtT)
            % Find index
            idx = find( ( CpT == uCpT(i3) ) & ( CtT == uCtT(i2) & ( CohT == uCohT(i1) )) ,1);
            if (~isempty(idx))
                addon = sprintf('%s/%.3f\t',Fxn{idx},Fyn(idx));
            else
                addon = sprintf('%s/-----\t',UNK);
            end
            line = [line,addon];
        end
        fprintf('%s\n',line)
        
        % Third line
        fprintf('        \t|\t\n')
        
    end
    fprintf('\n');
end