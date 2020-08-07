clear; close all;
% run on cuebuntu

wordlist = {'cortical stim','cort stim','cort. stim'};
wordlist_neg = {'ended','no','thanks','ready','done','between','over','end','ends','prepare','given','touching','testing','last','that was 7 ma','no ad for all','5 sec','during','out of romm','post','complete','rate','7 words','train','oct_25','oct_26','oct_28','2nd','without','duration','recall','Start Stimpres','sample'};

load('annotation_paths.mat');

% h5 dir
dir_h5 = '/nas_share/cuenap/data/h5_notch20';
dir_stamp = '/nas_share/cuenap/data/stamps';

% remove vim tmp files and expand multiples
annotation_paths2 = {};
sids2 = {};
prefix2 = {};
idx = 1;
for i = 1:length(annotation_paths)
    ap = annotation_paths{i};

    n_paths = count(ap,'/mnt/');
    if (n_paths > 1)
        apL = strsplit(ap,'/mnt/');
        apL = apL(2:end);
        for j = 1:n_paths
            apC{j} = ['/mnt/',apL{j}];
        end
    elseif (n_paths == 1)
        apC = {ap};
    end
 
    for k = 1:length(apC)
        ap = apC{k};
        if (endsWith(ap,'.txt'))
            %annotation_paths2{idx} = annotation_paths{i};
            annotation_paths2{idx} = ap;
            sids2{idx} = sids{i};
            prefix2{idx} = prefix{i};
            idx = idx + 1;
        end
    end
end
annotation_paths = annotation_paths2;
sids = sids2;
prefix = prefix2;

f_annot_idx = [];
f_annot_dt = {};
f_annot_str = {};
f_prefix = {};
f_sid = {};
f_annotation_paths = {};
i_c = 1;
for i = 1:length(annotation_paths)
    ap = annotation_paths{i};
    %fprintf('---\n')
    %fprintf('%s\n',ap)
    %fprintf('---\n')
   
    infile = fopen(ap,'r');
    D = textscan(infile,'%s','Delimiter','\n');
    fclose(infile);

    % parse creation date
    creation_date_raw = D{1}{4};
    creation_date_raw = strsplit(creation_date_raw);
    creation_day_raw = strjoin(creation_date_raw(4:end));
    creation_date_raw = strjoin(creation_date_raw(3:end));
    creation_day = datetime(creation_day_raw,'InputFormat','MMM dd, yyyy');
    creation_date = datetime(creation_date_raw,'InputFormat','HH:mm:ss MMM dd, yyyy');

    % build datetime for each annotation
    if (isnat(creation_date))
        fprintf('W: creation date could not be parsed. skipping.\n')
    else
        annotations_raw = D{1}(7:end);
        n_annot = length(annotations_raw);


        for i2 = 1:n_annot
            aL = strsplit(annotations_raw{i2});

            if (startsWith(aL{1},'d'))
                os_fmt = 0;
            else
                os_fmt = -1;
            end

            if (os_fmt == 0)
                os_d = str2double(aL{1 + os_fmt}(2:end)) - 1;
            else
                os_d = 0;
            end


            dstr = aL{2 + os_fmt};
            dstr = strsplit(dstr,'.');
            dstr = dstr{1};
            dstr = strsplit(dstr,':');
            os_h = str2double(dstr{1});
            os_m = str2double(dstr{2});
            os_s = str2double(dstr{3});
            annot_str = strjoin(aL((3 + os_fmt):end));
            annot_dt = creation_date + days(os_d) + hours(os_h) + minutes(os_m) + seconds(os_s);
            %f_annot_idx(i_c) = nan;



            if (contains(lower(annot_str),wordlist)) && (~contains(lower(annot_str),wordlist_neg))
                %fprintf('%s\t%s\t%s\t%s\n',sids{i},prefix{i},datestr(annot_dt),annot_str)

                % --- build h5file info ---
                h5fn = [dir_h5,'/',sids{i},'.h5'];
                h5_files = h5readatt(h5fn,'/h5eeg','files');
                h5_datetime = h5readatt(h5fn,'/h5eeg','datetime');
                h5_timestamp = h5readatt(h5fn,'/h5eeg','timestamp');
                h5_n_samples = h5readatt(h5fn,'/h5eeg','n_samples');
                h5_ds_factors = h5readatt(h5fn,'/h5eeg','ds_factors');
                h5_n_chan = h5readatt(h5fn,'/h5eeg/eeg','n_chan');
                h5_fs = h5readatt(h5fn,'/h5eeg/eeg','rate');

                % find file
                fidx = find(strcmp(h5_files,prefix{i}),1);
                if (isempty(fidx))
                    fprintf(2,'W: Could not find file %s: %s\n',sids{i},prefix{i})
                end
                f_start_dt = datetime(replace(h5_datetime{fidx},'T',' '));
                cs_samples = cumsum([0; h5_n_samples]) + 1;
                f_start_index = cs_samples(fidx);

                % find index where annotation exists
                ecog = H5eeg(h5fn);
                f_annot_idx_t0 = ecog.getSampleFromTime(dir_stamp, datestr(annot_dt));

                if (isempty(f_annot_idx_t0))
                    fprintf(2,'[*] Could not map "%s" for %s: %s, default to estimate\n',datestr(annot_dt),sids{i},prefix{i})
                    f_annot_idx_t = f_start_index + round(seconds(annot_dt - f_start_dt)*h5_fs);
                    isEst = true;
                else
                    f_annot_idx_t = f_annot_idx_t0;
                    isEst = false;
                end

                if ((f_annot_idx_t >= cs_samples(fidx + 1)) || isEst)
                    fprintf(2,'\tW: annotation datetime not found, skipping..\n');
                else
                    fprintf('%s\t%s\t%s\t%s\n',sids{i},prefix{i},datestr(annot_dt),annot_str)
                    f_annot_dt{i_c} = annot_dt;
                    f_annot_str{i_c} = annot_str;
                    f_prefix{i_c} = prefix{i};
                    f_sid{i_c} = sids{i};
                    f_annotation_paths{i_c} = annotation_paths{i};
                    f_annot_idx(i_c) = f_annot_idx_t;
                    i_c = i_c + 1;
                    %return
                end

                %working_dt = f_start_dt;
                %isfound = false;
                %stride = 250;
                %for i3 = 1:stride:h5_n_samples(fidx)
                %    if (annot_dt < working_dt)
                %        f_annot_idx = f_start_index + i3 - 1;
                %        isfound = true;
                %        break;
                %    end

                %    if (mod((i3 - 1),round(h5_n_samples(fidx)/50)) == 0)
                %        fprintf('\t%.2f%%\n',100*(i3 - 1)/h5_n_samples(fidx));
                %    end
                %    working_dt = working_dt + seconds( stride/(h5_fs) );%*h5_ds_factors(fidx)) );
                %end

                %if (~isfound)
                %    fprintf(2,'W: annotation datetime not found.\n');
                %end
                
                %return

            end

        end
    end

    %return
end


% Clean up
%nidx = isnan(f_annot_idx);
%annotation_paths = annotation_paths(~nidx);
%sids = sids(~nidx);
%prefix = prefix(~nidx);
%f_annot_idx = f_annot_idx(~nidx);
%f_annot_dt = f_annot_dt(~nidx);
%f_annot_str = f_annot_str(~nidx);

save('annotation_paths2_stim.mat','f_annotation_paths','f_sid','f_prefix','f_annot_idx','f_annot_dt','f_annot_str');
