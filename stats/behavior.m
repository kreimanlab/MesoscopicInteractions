clear; close all;

% Process behavior .index file
behaviors = {'BODY MOVEMENT','CONTACT','EATING','HEAD MOVEMENT',...
    'PATIENT IS TALKING','QUIET','SLEEP','SOMEONE IS TALKING','VIDEO GAMES',...
    'WATCH TV'};
f_in = fopen('sub21.index','r');
D = textscan(f_in,'%s','Delimiter','\n');
D = D{1};
beh_etime = zeros(1,length(D));
beh_str = cell(1,length(D));
for i = 1:length(D)
    d_raw = strsplit(D{i},',');
    beh_etime(i) = str2double(d_raw{1}(1:(end-1)));
    beh_str{i} = d_raw{2}(2:end);
end

% Parse behaviors
fn_h5 = '/media/klab/KLAB101/h5_notch20/sub21.h5';
startEtime = h5read(fn_h5,'/h5eeg/aux',[1 1],[1 1]);
N_COH = 53934; % number of coherence time points computed
thresh_warn_annot_hrs = 1;
w = 10;
width = 250;
n_behaviors = length(behaviors);
coh_behavior_idx = false(n_behaviors,N_COH);
for i = 1:n_behaviors
    beh = behaviors{i};
    state_beh_on = false;
    for j = 1:length(beh_str)
        cond_beh_on = contains(beh_str{j},beh) && endsWith(beh_str{j},'ON');
        cond_beh_off = contains(beh_str{j},beh) && endsWith(beh_str{j},'OFF');
        if ( cond_beh_on && (~state_beh_on) )
            beh_start_idx = beh_etime(j) - startEtime + 1;
            beh_start_idx = ceil(beh_start_idx/(w*width));
            state_beh_on = true;
        elseif ( cond_beh_on && (state_beh_on) )
            beh_start_idx = beh_etime(j) - startEtime + 1;
            beh_start_idx = ceil(beh_start_idx/(w*width));
        elseif ( cond_beh_off && (state_beh_on) )
            beh_end_idx = beh_etime(j) - startEtime + 1;
            beh_end_idx = ceil(beh_end_idx/(w*width));
            state_beh_on = false;
            if (beh_end_idx >= beh_start_idx)
                coh_behavior_idx(i,beh_start_idx:beh_end_idx) = true;
                %fprintf('start: %i\tend: %i\n',beh_start_idx,beh_end_idx);
                if ((beh_end_idx-beh_start_idx)*w > thresh_warn_annot_hrs*60*60)
                    fprintf(2,'\tW: %s annot longer than %.2f hrs\n',...
                        beh_str{j},thresh_warn_annot_hrs);
                end
            end
        end
    end
end

%return;

%fn_h5 = '/media/klab/KLAB101/sub21_4ce6814b.hdf5';
h5_info = h5info(fn_h5);

% H5eeg object
ecog = H5eeg(fn_h5);

% constants
etime_min = 96660604;
etime_max = 133239262;

% read subfiles
file_idx = 3; % the annotations are referenced to the third file: 4ce6814b
files = h5readatt(fn_h5,'/h5eeg','files');
n_samples = h5readatt(fn_h5,'/h5eeg','n_samples');
ds_factors = h5readatt(fn_h5,'/h5eeg','ds_factors');
n_samples_offset = sum(n_samples(1:(file_idx-1)));

% correct for downsampling, then the 2 files before annotated file
sidx = etime_min/ds_factors(file_idx) + n_samples_offset;
manual_shift = -47443; % offset found using trial-and-error
start_idx = round(sidx + manual_shift);
end_idx = round((sidx + manual_shift+(etime_max-etime_min+10)/2));
start_etime = h5read(fn_h5,'/h5eeg/aux',[1 start_idx],[1 1]);
end_etime = h5read(fn_h5,'/h5eeg/aux',[1 end_idx],[1 1]);
fprintf('start elapsed time index: %i\n',start_etime);
fprintf('end elapsed time index: %i\n',end_etime);
fprintf('start h5 index: %i\n',start_idx);
fprintf('end h5 index: %i\n',end_idx);

start_pos = [1 start_idx];
read_block = [ecog.n_chan (end_idx-start_idx+1)];
bytes = 4;
fprintf('[*] Reading %.2f GBytes..\n',prod(read_block)*bytes*(1e-9))
data = h5read(fn_h5,'/h5eeg/eeg',start_pos,read_block);
aux = h5read(fn_h5,'/h5eeg/aux',start_pos,[1 read_block(2)]);

fprintf('[*] Saving..\n')
save('sub21_24hr','-v7.3','-nocompression')
fprintf('[!] Done.\n')



