clear; close all;

% Patients
Subjects = {'m00001','m00003','m00005','m00006','m00019','m00021','m00022','m00023',...
    'm00024','m00025','m00026','m00027','m00028','m00030','m00032','m00033',...
    'm00035','m00037','m00038','m00039','m00043','m00044','m00045','m00047',...
    'm00048','m00049','m00052','m00053','m00055','m00056','m00058','m00059',...
    'm00060','m00061','m00068','m00071','m00073','m00075','m00079','m00083',...
    'm00084','m00095','m00096','m00097','m00100','m00107','m00122','m00124',...
    };

dir_sub = '../h5_notch20';
dir_raw = '/mnt/symphony_1/RawData';
dir_raw1 = '/mnt/symphony_2/RawData';
dir_raw2 = '/mnt/synology_1/rawdata';
dir_raw3 = '/mnt/synology_2/RawData';

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
