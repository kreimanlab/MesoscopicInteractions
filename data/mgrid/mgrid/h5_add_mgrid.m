close all;
clear;

%addpath('/media/jerry/internal/data/reverse')

%subjects_dir = '/media/jerry/internal/data/coreg';
subjects_dir = '../coregistration';
fns_h5 = carveF('../h5_notch20','.h5');
fns_h5 = fns_h5(~contains(fns_h5,'mSu'));
elec_space = 10; % mm

% % ----------------- BYPASS -------------------------------
% fns_h5 = {'h5_notch20/sub48.h5','h5_notch20/sub34.h5'};

for i = 1:length(fns_h5)
    fn_h5 = fns_h5{i};
    ecog = H5eeg(fn_h5);

    % extract patient name
    fn_h5_tmp = strsplit(fn_h5,'/');
    fn_h5_tmp = strsplit(fn_h5_tmp{end},'.h5');
    sid = fn_h5_tmp{1};
    %sid = fn_h5((end-8):(end-3));

    fprintf('[*] %s\n',sid)
    fn_mgrid = sprintf('%s/%s/elec_recon/%s.mgrid',subjects_dir,sid,sid);
    if (~exist(fn_mgrid,'file'))
        fprintf(2,'[!] Skip: did not find mgrid file for: %s.\n',sid);
    else
        f_mgrid = fopen(fn_mgrid,'r');
        D = textscan(f_mgrid,'%s','Delimiter','\n');
        D = D{1};
        fclose(f_mgrid);
        
        % parse
        g_names = {};
        g_dims = [];
        n_skip_description = 1;
        for j = 1:length(D)
            l = D{j};
            if (contains(l,'Description'))
                if (n_skip_description > 0)
                    n_skip_description = n_skip_description - 1;
                else
                    g_names = [g_names; D{j+1}];
                end
            end
            if (contains(l,'Number of Grids'))
                n_grids = str2double(D{j+1});
            end
            if (contains(l,'Dimensions'))
                dims = strsplit(D{j+1});
                dims = [str2double(dims{1}) str2double(dims{2})];
                g_dims = [g_dims; dims];
            end
            
        end
        
        
        % match to existing channel labels
        chan_labels = h5readatt(fn_h5,'/h5eeg/eeg','labels');
        if (strcmp(sid,'sub28'))
            chan_labels{59} = 'PF59';
        end
        %chan_hemi = h5readatt(fn_h5,'/h5eeg/eeg','hemi');
        
        % build mgrid channel and hemisphere labels
        mgrid_labels = {};
        mgrid_hemi = {};
        mgrid_labels_i = 1;
        mgrid_gnames = {};
        mgrid_gnums = [];
        for ii = 1:n_grids
            n_grid_chans = prod(g_dims(ii,:));
            mgrid_gname = strsplit(g_names{ii},'_');
            for iii = 1:n_grid_chans
                mgrid_labels{mgrid_labels_i,1} = sprintf('%s%i',mgrid_gname{2},iii);
                mgrid_hemi{mgrid_labels_i,1} = mgrid_gname{1}(1);
                mgrid_gnums(mgrid_labels_i) = iii;
                mgrid_labels_i = mgrid_labels_i + 1;
            end
            mgrid_gnames{ii,1} = mgrid_gname{2};
        end
        
        % check they are the same size
        if (length(mgrid_labels) ~= length(chan_labels))
            fprintf('[ERROR] Subject %s chan_labels does not match .mgrid file.\n',sid)
        else
            %h5writeatt(fn_h5,'/h5eeg/eeg','isRH',uint8(strcmp(mgrid_hemi,'R')));
        end
        
        if (all(strcmp(mgrid_labels,chan_labels)) && false)
            fprintf('[*] .mgrid electrode names match .h5, no mapping needed.\n')
        else
            fprintf('[!] Channel labels do not match. mapping needed.\n')
            
            chan_gnames = cell(length(chan_labels),1);
            chan_gnums = zeros(length(chan_labels),1);
            for jj = 1:length(chan_labels)
                mlab = chan_labels{jj};
                m = isstrprop(mlab,'alpha');
                chan_gnames{jj,1} = mlab(m);
                chan_gnums(jj,1) = str2double(mlab(~m));
            end
            
            u_chan_gnames = unique(chan_gnames);
            n_grids_equal = (length(u_chan_gnames) == length(mgrid_gnames));
            mgrid_gnames = sort(mgrid_gnames);
            u_chan_gnames = sort(u_chan_gnames);
            n_elec_per_grid_equal = false(length(mgrid_gnames),1);
            for jj = 1:length(mgrid_gnames)
                n_elec_per_grid_mgrid = sum(contains(mgrid_labels,mgrid_gnames{jj}));
                n_elec_per_grid_chan = sum(contains(chan_labels,u_chan_gnames{jj}));
                n_elec_per_grid_equal(jj) = (n_elec_per_grid_mgrid == n_elec_per_grid_chan);
            end
            
            if (~all(n_elec_per_grid_equal))
                fprintf(2,'[!] Number of electrodes per grid not equal.\n')
            end
            %return
            
            if (n_grids_equal)
            else
                fprintf(2,'[!] Number of grids in .mgrid not the same as .h5\n')
            end
%             chan_gnames = {};
%             for jj = 1:length(chan_labels)
%             end
            %sreturn
            
            % Build theoretical distance matrix
            n_chan = length(chan_labels);
            n_bchan = ecog.n_bchan;
            Dmat = zeros(n_bchan,n_bchan);
            for i2 = 1:(n_bchan-1)
                for j2 = (i2+1):n_bchan
                    b1c1 = ecog.bip(i2,1);
                    b1c2 = ecog.bip(i2,2);
                    b2c1 = ecog.bip(j2,1);
                    b2c2 = ecog.bip(j2,2);
                    
                    pt_b1c1 = ecog.bip(i2,4:6);
                    pt_b1c2 = ecog.bip(i2,7:9);
                    pt_b2c1 = ecog.bip(j2,4:6);
                    pt_b2c2 = ecog.bip(j2,7:9);
                    dists = zeros(2,2);
                    dists(1,1) = sqrt(sum((pt_b1c1 - pt_b2c1).^2));
                    dists(1,2) = sqrt(sum((pt_b1c1 - pt_b2c2).^2));
                    dists(2,1) = sqrt(sum((pt_b1c2 - pt_b2c1).^2));
                    dists(2,2) = sqrt(sum((pt_b1c2 - pt_b2c2).^2));
                    
                    bIdx = [b1c1; b1c2; b2c1; b2c2];
                    
                    Hemi = mgrid_hemi(bIdx);
                    hemi_b1c1 = mgrid_hemi{b1c1};
                    hemi_b1c2 = mgrid_hemi{b1c2};
                    hemi_b2c1 = mgrid_hemi{b2c1};
                    hemi_b2c2 = mgrid_hemi{b2c2};
                    
                    Grid = chan_gnames(bIdx);
                    grid_b1c1 = chan_gnames{b1c1};
                    grid_b1c2 = chan_gnames{b1c2};
                    grid_b2c1 = chan_gnames{b2c1};
                    grid_b2c2 = chan_gnames{b2c2};
                    
                    Gnum = mgrid_gnums(bIdx);
                    gnum_b1c1 = mgrid_gnums(b1c1);
                    gnum_b1c2 = mgrid_gnums(b1c2);
                    gnum_b2c1 = mgrid_gnums(b2c1);
                    gnum_b2c2 = mgrid_gnums(b2c2);
                    
                    % Get grid dimensions
                    gdim_b1c1 = g_dims(contains(mgrid_gnames,grid_b1c1),:);
                    gdim_b1c2 = g_dims(contains(mgrid_gnames,grid_b1c2),:);
                    gdim_b2c1 = g_dims(contains(mgrid_gnames,grid_b2c1),:);
                    gdim_b2c2 = g_dims(contains(mgrid_gnames,grid_b2c2),:);
                    Gdim = [gdim_b1c1; gdim_b1c2; gdim_b2c1; gdim_b2c2;];
                    
                    sameHemi = [[strcmp(hemi_b1c1,hemi_b2c1), strcmp(hemi_b1c1,hemi_b2c2)]; ...
                                [strcmp(hemi_b1c2,hemi_b2c1), strcmp(hemi_b1c2,hemi_b2c1)]];
                    
                    sameGrid = [[strcmp(grid_b1c1,grid_b2c1), strcmp(grid_b1c1,grid_b2c2)]; ...
                                [strcmp(grid_b1c2,grid_b2c1), strcmp(grid_b1c2,grid_b2c1)]];
                            
                    d = min(dists(:));   
                    if (not(all(sameHemi(:))))
                        d = (-1) * Inf;
                    elseif (all(sameGrid))
                        dim = gdim_b1c1;
                        C1_num = [gnum_b1c1 gnum_b1c2];
                        C2_num = [gnum_b2c1 gnum_b2c2];
                        Is_n = zeros(length(C1_num),length(C2_num));
                        for ic1 = 1:length(C1_num)
                            c1_num = C1_num(ic1);
                            for ic2 = 1:length(C2_num)
                                c2_num =  C2_num(ic2);
                                % determine if diagonal neighbors
                                % Are they immediate neighbors?
                                % neighbor if: electrodes are consecutive, and the
                                % preceding electrode [min] does is not the last electrode
                                % in the long dimension [dim(2)]
                                is_n1 = (abs(c1_num - c2_num) == 1) & (mod(min([c1_num,c2_num]),dim(2)) ~= 0);
                                %is_n1 = (abs(c1_num - c2_num) == 1);
                                %is_n1 = false;
                                % or:
                                
                                % Are they diagonal neighbors?
                                if (dim(1) > 1)
                                    is_n1 = is_n1 | ((max([c1_num,c2_num]) - dim(2)) == min([c1_num,c2_num]));

                                    % check top-right, bottom-left
                                    is_n2 = ((c1_num + dim(2) + 1) == c2_num);
                                    is_n2 = is_n2 | ((c2_num + dim(2) + 1) == c1_num);
                                    % check top-left, bottom-right
                                    is_n2 = is_n2 | ((c1_num + dim(2) - 1) == c2_num);
                                    is_n2 = is_n2 | ((c2_num + dim(2) - 1) == c1_num);
                                else
                                    is_n2 = false;
                                end
                                
                                Is_n(ic1,ic2) = (is_n1 | is_n2);
                            end
                        end
                        
                        if (any(Is_n(:)))
                            d = -100;

                            %fprintf('[!] is any n: %i %i; %i %i\n',b1c1,b1c2,b2c1,b2c2)

                        end
                    end
                    
                    Dmat(i2,j2) = d;
                    Dmat(j2,i2) = d;
                end
            end
            
            % DEBUG
%             [~,sIdx] = sort(ecog.bip(:,1));
%             imagesc(Dmat(sIdx,sIdx))
%             set(gca,'YTick',1:length(chan_labels(ecog.bip(sIdx,1))))
%             set(gca,'YTickLabel',chan_labels(ecog.bip(sIdx,1)))
            
            % Save Dmat to file
            fprintf('\t(*) Saving Dmat to file..\n');
            h5create(fn_h5,'/h5eeg/elec_recon/Dmat',size(Dmat));
            h5write(fn_h5,'/h5eeg/elec_recon/Dmat',Dmat);
            h5create(fn_h5,'/h5eeg/elec_recon/g_dims',size(g_dims));
            h5write(fn_h5,'/h5eeg/elec_recon/g_dims',g_dims);
            fprintf('\t\tDone.\n');
            %return
        end
        
        %return
    end
end
