close all;
%clear;

if ismac
    resultsDir = '/Volumes/RawData/data/results';
    h5Dir = '/Volumes/RawData/scripts/synth/out';
elseif isunix
    [~,hname] = system('hostname');
    if strcmp(strip(hname),'hopper')
        %resultsDir = '/media/klab/44/data/results';
        %resultsDir = '/media/klab/D0BA8713BA86F56E/data/results';
        %resultsDir = '~/data/h5eeg/artifact_removed/test/opencl/remote/results';
        resultsDir = '~/data/h5eeg/artifact_removed/test/opencl/remote/results/scipy';
        %resultsDir = '~/data/h5eeg/artifact_removed/test/opencl/remote/cpu';
        h5Dir = '/media/klab/44/h5';
    else
        resultsDir = '/mnt/cuenap2/data/results';
        h5Dir = '/mnt/cuenap2/scripts/synth/out';
    end
end

if (~exist('A','var'))
    A = Analysis(resultsDir,h5Dir);
end

sid = 'sub2';%'sub12';%'sub4';%'sub2';% 'sub1';
metric = 'pcBroadband';
bchan1_str = 'IF33-IF34';%'AT1-AT2';%'RMT4-RMT5';%'LA10-LA11';%'IF33-IF34';%'LF53-LF54';
bchan2_str = 'IF36-IF37';%'AT6-AT7';%'RMT7-RMT8';%'LP26-LP27';%'IF36-IF37';%'PT39-PT40';
bro_s = 0.5; % Hz
bro_e = 125; % Hz

fname_graph = sprintf('%s_graph-%s.h5',sid,metric);
file_i = find(strcmp(A.r.graph_f,fname_graph),1);
subj_i = find(strcmp(A.h5eeg.subjects,sid),1);
bchan1 = find(strcmp(A.h5eeg.bip_labels{subj_i},bchan1_str),1);
bchan2 = find(strcmp(A.h5eeg.bip_labels{subj_i},bchan2_str),1);
n_bchan = A.h5eeg.n_bchan{subj_i};

i_r = 1;
for i = 1:(n_bchan-1)
    for j = (i+1):n_bchan
        if ((i == bchan1) && (j == bchan2))
            %return
            % Load GPU result
            fname = [resultsDir,'/',fname_graph];
            h5inf = h5info(fname,'/R');
            n_w = h5inf.Dataspace.Size(2);
            r_gpu = sqrt(h5read(fname,'/R',[i_r 1],[1 n_w]));
            chan1 = h5read(fname,'/chan1',[i_r 1],[1 1])+1;
            chan2 = h5read(fname,'/chan2',[i_r 1],[1 1])+1;
            w = h5read(fname,'/w');
            
            % Recompute coherence
            h5fname = [h5Dir,'/',sid,'.h5'];
            fs = A.h5eeg.fs{subj_i};
            n_samples = A.h5eeg.n_samples{subj_i};
            w_samp = double(round(fs*w));
            bip = A.h5eeg.bip{subj_i};
            Starts = 1:w_samp:n_samples;
            r_cpu = nan(1,length(Starts));
%             fprintf('Reading..\n')
%             v_b1c1 = h5read(h5fname,'/h5eeg/eeg',[bip(chan1,1), 1],[1 n_samples]);
%             fprintf('Finished 1/4.\n')
%             v_b1c2 = h5read(h5fname,'/h5eeg/eeg',[bip(chan1,2), 1],[1 n_samples]);
%             fprintf('Finished 2/4.\n')
%             v_b2c1 = h5read(h5fname,'/h5eeg/eeg',[bip(chan2,1), 1],[1 n_samples]);
%             fprintf('Finished 3/4.\n')
%             v_b2c2 = h5read(h5fname,'/h5eeg/eeg',[bip(chan2,2), 1],[1 n_samples]);
%             fprintf('Finished 4/4.\n')
%             v1 = v_b1c1 - v_b1c2;
%             v2 = v_b2c1 - v_b2c2;
            
            for k = 1:length(Starts)
                tic;
                start_i = double(Starts(k));
                end_i = start_i + w_samp - 1;
                if (end_i > n_samples)
                    end_i = n_samples;
                end
                v_b1c1 = h5read(h5fname,'/h5eeg/eeg',[bip(chan1,1), start_i],[1 (end_i-start_i+1)]);
                v_b1c2 = h5read(h5fname,'/h5eeg/eeg',[bip(chan1,2), start_i],[1 (end_i-start_i+1)]);
                v_b2c1 = h5read(h5fname,'/h5eeg/eeg',[bip(chan2,1), start_i],[1 (end_i-start_i+1)]);
                v_b2c2 = h5read(h5fname,'/h5eeg/eeg',[bip(chan2,2), start_i],[1 (end_i-start_i+1)]);
                v1 = v_b1c1 - v_b1c2;
                v2 = v_b2c1 - v_b2c2;
                %return
                %r = corr(v1(start_i:end_i)',v2(start_i:end_i)');
                %[cxy,f] = mscohere(v1',v2',round(fs),round(0.5*fs),round(fs)*2,round(fs));
                %[cxy,f] = mscohere(v1',v2',2*round(fs),[],[],round(fs));
                [cxy,f] = mscohere(v1',v2',hamming(2*round(fs)),round(0.5*round(fs)),[],round(fs));
                
                
                r_cpu(k) = mean(sqrt(cxy((f >= bro_s) & (f < bro_e))));
                t_sing = toc;
                fprintf('(%i/%i) r = %.4f\tETA: %.2f min\n',k,length(Starts),r_cpu(k),(length(Starts)-k)*t_sing/60);
                
                if true%(r_cpu(k) > 0.3)
                    plot(f,sqrt(cxy),'black');
                    hold on
                    plot([min(f) max(f)],[r_cpu(k) r_cpu(k)],'red--')
                    xlabel(sprintf('Frequency (Hz), res: %.2f Hz, nfft: %i',f(2)-f(1),length(f)))
                    ylabel('Coherence')
                    axis([min(f) max(f) 0 1])
                    return;
                end
                
%                 if (k == 2000)
%                     break
%                 end
            end
            
            h = figure;
            subplot(2,1,1)
            plot(r_gpu,'black.')
            xlabel('Minutes')
            ylabel('Broadband Coherence GPU')
            axis([0 n_w -0.8 0.8])
            subplot(2,1,2)
            plot(r_cpu,'black.')
            xlabel('Minutes')
            ylabel('Broadband Coherence CPU')
            axis([0 n_w -0.8 0.8])
            print(h,sprintf('figures/check_coherence/%s_%s_%s_%s',sid,metric,bchan1_str,bchan2_str));
            
            return
        end
        
        i_r = i_r + 1;
    end
end

