close all;
clear;

art_idx_dir = 'verify/art_idx';
art_idx_dir_out = 'art';

outL = carveF(art_idx_dir, '.mat');

system(['mkdir ',art_idx_dir_out]);

% for i = 1:length(outL)
for i = length(outL)
    af = outL{i};
    sid = strsplit(af,'.mat');
    sid = strsplit(sid{1},'/');
    sid = sid{end};
    load(af);
    
    if (strcmp(sid,'mSu'))
        % remove artifacts for mSu
        Art(Art == 3) = 0;
    end
    
%     h = figure;
%     AxesH = axes;
%     drawnow;
%     InSet = get(AxesH, 'TightInset');
%     InSet = InSet * 1.5;
%     InSet(1) = InSet(1) * 1.1; %ytick label offset
%     set(AxesH, 'Position', [InSet(1:2), 1-InSet(1)-InSet(3), 1-InSet(2)-InSet(4)])
    for j = 1:length(bchan_labels)
        Aj = Art(j,:);
        Ajp = (j-1) * ones(size(Aj));
        T = linspace(0,length(Aj)*((1)/(3600)),length(Aj)); % hours
        [~,n_arts] = size(Aj);
        
        % Split by time segment
%         sec_time_seg = 60;
%         [n_bchan,n_art] = size(Aj);
%         Starts = 1:sec_time_seg:n_art;

        
%         for k = 0:3
%             plot(T(Aj==k),Ajp(Aj==k),'.','MarkerSize',20,'Color',Y_col(k+1,:));
%             hold on;
%         end
        
    end
%     axis tight
%     xlabel('Hours')
%     yticks(1:n_bchan)
%     yticklabels(bchan_labels)
%     return

    % Horizontal and vertical artifact fractions
    horiz = sum(Art~=0,2)/length(Aj);
    vert = sum(Art~=0)/n_bchan;
    [horiz_s,horiz_si] = sort(horiz,'Descend');
    [vert_s,vert_si] = sort(vert,'Descend');
    
    fprintf('%s\n',sid)
    fprintf('\thoriz > 0.5: %i of %i\n',sum(horiz>0.5),n_bchan)
    fprintf('\tvert > 0.5: %.2f %%\n',100*sum(vert>0.5)/length(vert))
    
    % Apply horiz and vert to artifact files
    Art(horiz>0.5,:) = 4;
    Art(:,vert>0.5) = 5;
%     if (strcmp(sid,'m00001'))
%         imagesc(Art); colorbar;
%         return
%     end
    
    if (true)
        % hdf5 save
        % art_width
        outfn = sprintf('%s/%s_art.h5',art_idx_dir_out,sid);
        if (exist(outfn,'file'))
            system(['rm -f ',outfn]);
        end
        h5create(outfn,'/artifacts',size(Art),'Datatype','uint8')
        h5writeatt(outfn,'/artifacts','width',art_width);
        h5writeatt(outfn,'/artifacts','n_bchan',n_bchan);
        h5writeatt(outfn,'/artifacts','n_arts',n_arts);
        h5writeatt(outfn,'/artifacts','version',3);
        h5write(outfn,'/artifacts',Art,[1 1],size(Art));

        % compute legacy artifact labels (for perm.py)
        arts = zeros(2,n_arts);
        % copy starting indices
        arts_o = h5read(['/mnt/cuenap/data/h5_notch20/',sid,'.h5'],'/h5eeg/artifacts');
        arts(2,:) = arts_o(2,1:n_arts);
        % compute if artifact
        arts(1,:) = n_bchan*(sum(Art(:,1:n_arts) ~=0) > 0);
        h5create(outfn,'/artifacts_v1',size(arts),'Datatype','single')
        h5write(outfn,'/artifacts_v1',single(arts),[1 1],size(arts))
        h5writeatt(outfn,'/artifacts_v1','width',art_width)
        h5writeatt(outfn,'/artifacts_v1','n_arts',n_arts)
        h5writeatt(outfn,'/artifacts_v1','n_bchan',n_bchan)
        h5writeatt(outfn,'/artifacts_v1','version',1)
    %     for ii = 1:n_arts
    %         arts(1,ii) = ;
    %     end



        % --- re-compute art_idx ----------------------------------------------

        % Calculate usable electrodes
        n_comb = nchoosek(n_bchan,2);
        n_w = round(n_arts/w) - 1;
        art_idx = false(n_comb,n_w);
        for l = 1:n_w
            l_start = (l-1)*w + 1;
            l_end = l_start + w - 1;

            % Precompute conditions
            condp = zeros(n_bchan,1);
            for k4 = 1:n_bchan
                condp(k4) = sum(Art(k4,l_start:l_end) ~= 0) > 0;
            end

            k3 = 1;
            for k1 = 1:(n_bchan-1)

                for k2 = (k1+1):n_bchan

                    art_idx(k3,l) = (condp(k1) || condp(k2));

                    %art_idx(k3,l) = sum([ak1,Art(k2,l_start:l_end)]) > 0;
                    k3 = k3 + 1;
                end
            end
        end
        frac_art = sum(sum(art_idx ~= 0))/numel(art_idx);

        fprintf('%s\tw: %.1f\t%.2f %% artifacts\n',...
            sid,w,100*frac_art)

        h5create(outfn,'/art_idx',size(art_idx),'Datatype','uint8');
        h5write(outfn,'/art_idx',uint8(art_idx),[1 1],size(art_idx));
        h5writeatt(outfn,'/art_idx','w',w);
        h5writeatt(outfn,'/art_idx','frac_art',frac_art);
    end
    
    
%     save(sprintf('%s/%s.mat',art_idx_dir_out,sid),'Art','-append')
    
    
    %disp(horiz_s(1:20)')
    %disp(vert_s(1:20))
    
end