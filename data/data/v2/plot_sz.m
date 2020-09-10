clear; close all;

% Patients
Subjects = {'example'};

dir_sub = '../h5_notch20';
dir_raw = '../txt';
dir_raw1 = '../txt';
dir_raw2 = '../txt';
dir_raw3 = '../txt';

annotation_paths = {};
sids = {};
prefix = {};
a_idx = 1;
for i = 1:length(Subjects)
    sid = Subjects{i};
    subf = h5readatt([dir_sub,'/',sid,'.h5'],'/h5eeg','files');
    fprintf('%s\n',sid)
    for j = 1:length(subf)
        % Look for .ent.txt
        [~,fn_ent] = system(sprintf('find %s/%s/ | grep -i %s | grep ".ent.txt"',...
            dir_raw,sid,subf{j}));
        if (isempty(fn_ent) || contains(fn_ent,'No such file'))
            [~,fn_ent] = system(sprintf('find %s/%s/ | grep -i %s | grep ".ent.txt"',...
                dir_raw1,sid,subf{j}));
        end
        if (isempty(fn_ent) || contains(fn_ent,'No such file'))
        %if (isempty(fn_ent) && ~contains(fn_ent,'No such file'))
            [~,fn_ent] = system(sprintf('find %s/%s/ | grep -i %s | grep ".ent.txt"',...
                dir_raw2,sid,subf{j}));
        end
        if (isempty(fn_ent) || contains(fn_ent,'No such file'))
            [~,fn_ent] = system(sprintf('find %s/%s/ | grep -i %s | grep ".ent.txt"',...
                dir_raw3,sid,subf{j}));
        end

        % Report
        if (isempty(fn_ent) || contains(fn_ent,'No such file'))
            fprintf('\t%s\n',subf{j});
        else
            fn_ent = fn_ent(1:(end-1));
            fprintf('\t%s - found\n',subf{j});
            fprintf('\t\t%s\n',fn_ent);
            annotation_paths{a_idx} = fn_ent;
            sids{a_idx} = sid;
            prefix{a_idx} = subf{j};
            a_idx = a_idx + 1;
        end
    end
end

save('annotation_paths','annotation_paths','sids','prefix');
