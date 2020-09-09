close all;
clear;



dir_h5 = '/media/klab/internal/data/h5_notch20'; %'/media/klab/KLAB101/h5_notch20';
dir_art = sprintf('%s/art_nosz',dir_h5); %dir_artLp = '/media/klab/KLAB101/h5_notch20/art_nosz';
dir_res = '/home/jerry/data/results/coh_w10_phase'; %dir_resLp = '/media/klab/KLAB101/results/coh_w10'; 
dir_res_perm = '/home/jerry/data/opencl/results_res5hz'; 
dir_cor = '/media/klab/internal/data/coreg';
dir_cache = './cache';
dir_stamp = '/media/klab/internal/data/stamps';
dir_video = '/media/klab/internal/data/videos';

metrics = {'pcBroadband','pcTheta','pcAlpha','pcBeta','pcGamma'};

% Patients
Subjects = {'sub1','sub2','sub3','sub4','sub5','sub6','sub7','sub8',...
    'sub9','sub10','sub11','sub12','sub13','sub14','sub15','sub16',...
    'sub17','sub18','sub19','sub20','sub21','sub22','sub23','sub24',...
    'sub25','sub26','sub27','sub28','sub29','sub30','sub31','sub32',...
    'sub33','sub34','sub35','sub36','sub37','sub38','sub39','sub40',...
    'sub41','sub42','sub43','sub44','sub45','sub46','sub47','sub48',...
    'mSu'};

% Exclude monkey
Subjects = Subjects(1:(end-1));

Subjects = {'sub3'};


trig_sim = false;
if (trig_sim)
    n_sim = 1000;
    fs = 250;
    w = 10;
    
    pli = 60;
    PLI_S = pli-4;
    PLI_E = pli+4;
    PL2_S = (fs-(3*pli))-2;
    PL2_E = (fs-(3*pli))+2;
    PL3_S = (2*pli)-3;
    PL3_E = (2*pli)+3;
    GAM_S = 30;
    GAM_E = 100;
    
    %Coh = zeros(2,n_sim);
    
    Cohs = [];
    Freq_res = [0.5,5];
    for j = 1:length(Freq_res)
        freq_res = Freq_res(j);
        Coh = [];
        window = hamming(round((1/freq_res)*fs));
        parfor i = 1:n_sim
            v1 = normrnd(0,1,[1 fs*w]);
            v2 = normrnd(0,1,[1 fs*w]);
            %v1 = rand([1 fs*w]);
            %v2 = rand([1 fs*w]);
            %[C12,f] = mscohere(v1,v2,[],[],[],fs);
            [P12,~] = cpsd(v1,v2,window,[],[],fs);
            [P11,~] = cpsd(v1,v1,window,[],[],fs);
            [P22,f] = cpsd(v2,v2,window,[],[],fs);
            mask_gam = ((f >= GAM_S) & (f < PLI_S)) | ((f > PLI_E) & (f < PL2_S)) | ((f > PL2_E) & (f < GAM_E));
            CSD12 = P12./sqrt(P11 .* P22);
            ph = atan2(imag(CSD12),real(CSD12));
            r = abs(CSD12);
            %Coh(1,i) = mean(r(mask_gam));
            %Coh(2,i) = mean(ph(mask_gam));
            Coh = [Coh, [mean(r(mask_gam)); mean(ph(mask_gam))]];
        end
        Cohs = [Cohs; Coh];
    end
    
    
    
    h = figure('visible','off','PaperUnits','inches','PaperPosition',[0 0 4 4]);
    hold all;
    
    for k = 1:length(Freq_res)
        rang = (k*2 - 1):(k*2);
        Coh = Cohs(rang,:);
        Rg = Coh(1,:);
        PHg = Coh(2,:);
        x = Rg .* cos(PHg);
        y = Rg .* sin(PHg);
        x = x(:);
        y = y(:);
        
        if k ~= 1
            plot(x,y,'.','MarkerSize',1,'Color',[0 0.6 0]);
        else
            plot(x,y,'black.','MarkerSize',1);
        end
    
    end
    
    % center marker
    plot([0],[0],'black+','MarkerSize',10);

    axis(0.5*[-1 1 -1 1]);
    daspect([1 1 1]);
    xlabel('Re(Coherence)');
    ylabel('Im(Coherence)');
    set(gca,'TickDir','out');
    box off;
    xticks(linspace(-0.5,0.5,5));
    yticks(linspace(-0.5,0.5,5));
    
    sid = Subjects{1};
    metric = metrics{1};
    print(h,sprintf('figures/coherence_phase_p5_5_%s_%s',sid,metric),'-depsc');
    close(h);
        
end


for iM = 1
    for i = 1:length(Subjects)
        sid = Subjects{i};
        metric = metrics{iM};
        fn_graph = sprintf('%s/%s_graph-%s.h5',dir_res,sid,metric);
        fn_perm = sprintf('%s/5hz/%s_perm-%s-1000.h5',dir_res_perm,sid,metric);
        PHg = h5read(fn_graph,'/PH');
        Rg = h5read(fn_graph,'/R');
        PHp = h5read(fn_perm,'/PH');
        Rp = h5read(fn_perm,'/R');
        
        h = figure('visible','off','PaperUnits','inches','PaperPosition',[0 0 4 4]);
        x = Rg .* cos(PHg);
        y = Rg .* sin(PHg);
        x = x(:);
        y = y(:);
        x = x((Rg(:))<0.99);
        y = y((Rg(:))<0.99);
        n_cut = 10000; %numel(x);
        plot(x(1:n_cut),y(1:n_cut),'black.','MarkerSize',1);
        hold all;
        x = Rp .* cos(PHp);
        y = Rp .* sin(PHp);
        x = x(:);
        y = y(:);
        clear Rp;
        clear PHp;
        plot(x,y,'.','Color',[0 0 0.6],'MarkerSize',1);
        % center marker
        plot([0],[0],'black+','MarkerSize',10);
        
        axis([-1 1 -1 1]);
        daspect([1 1 1]);
        xlabel('Re(Coherence)');
        ylabel('Im(Coherence)');
        set(gca,'TickDir','out');
        box off;
        xticks(linspace(-1,1,5));
        yticks(linspace(-1,1,5));
        
        print(h,sprintf('figures/coherence_phase_%s_%s',sid,metric),'-depsc');
        close(h);
        
        %p = polarhistogram(PHg(:),1000,'BinLimits');
        return;
    end
end