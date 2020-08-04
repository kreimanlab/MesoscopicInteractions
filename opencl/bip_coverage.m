close all;

DIST_THRESH_MM = 15;
%Subjects = {'m00006','m00019','m00023','m00024','m00026','m00030','m00037','m00068'};
SubjectsAll = {'m00006','m00019','m00023','m00024','m00026','m00030','m00037','m00038','m00043','m00060','m00068','m00083'};
hcc = jet(length(SubjectsAll));
h = figure;
PctCovs = zeros(1,length(SubjectsAll));

for k = 1:length(SubjectsAll)

    Subjects = SubjectsAll(1:k);
    
    Adj_cov = zeros(180,180);
    Adj_sub = cell(180,180);

    for i = 1:length(Subjects)
        sid = Subjects{i};
        elecLabel = read_label(sprintf('%',sid),sprintf(...
            '/home/klab/coregistration/%s/label/all_surf_ielvis.label',sid));
        load(sprintf('%s/label/all_surf_ielvis_mmp.mat',sid));
        eBip = bip(sid);
        [n_bchan,~] = size(eBip);
        eBipA = zeros(n_bchan,2);
        for j = 1:n_bchan
            e = eBip(j,1:2);
            ss1 = strsplit(AtlasLabels{1}{e(1)},'_');
            ss2 = strsplit(AtlasLabels{1}{e(2)},'_');

            % unknown electrodes
            if (length(ss1) == 1)
                atIdx1 = 0;
            else
                atIdx1 = find(contains(cLH.struct_names,sprintf('L_%s_ROI',ss1{2})));
            end
            if (length(ss2) == 1)
                atIdx2 = 0;
            else
                atIdx2 = find(contains(cLH.struct_names,sprintf('L_%s_ROI',ss2{2})));
            end
            eBipA(j,:) = [atIdx1,atIdx2];
        end

        % Map to coverage matrix
        for m = 1:n_bchan
            for n = 1:n_bchan
                for a1 = 1:2
                    for a2 = 1:2
                        % Distance thresholding
                        cBip = [eBip(m,4:end),eBip(n,4:end)];
                        cBip = reshape(cBip,4,3);
                        Adist = zeros(4,4);
                        for d1 = 1:4
                            for d2 = 1:4
                                Adist(d1,d2) = sqrt(sum((cBip(d1,:) - cBip(d2,:)).^2));
                            end
                        end
                        skip = (min(Adist(~(Adist==0))) > DIST_THRESH_MM);
                        if (~skip)
                            try
                                Adj_cov(eBipA(m,a1),eBipA(n,a2)) = Adj_cov(eBipA(m,a1),eBipA(n,a2)) + 1;
                                Adj_sub{eBipA(m,a1),eBipA(n,a2)} = [Adj_sub{eBipA(m,a1),eBipA(n,a2)},i];
                            catch
                                %fprintf('Skip. \n')
                            end
                        end
                    end
                end
            end
        end
    end

    % Number of subjects per ROI pair    
    Adj_nsub = zeros(size(Adj_sub));
    for m = 1:length(Adj_sub)
        for n = 1:length(Adj_sub)
            Adj_nsub(m,n) = length(unique(Adj_sub{m,n}));
        end
    end
    [f,x] = hist(Adj_nsub(Adj_nsub ~= 0),length(Subjects));
    plot(1:length(Subjects),f,'-o','Color',hcc(k,:)); hold on;

    pctCov = 100*(sum(sum(Adj_cov ~= 0))/(numel(Adj_cov)));
    PctCovs(k) = pctCov;
    fprintf('Pct coverage: %.2f\n',pctCov)

end
xlabel('Number of subjects')
ylabel('Number of ROI pairs')
legStr = cell(1,length(SubjectsAll));
for i = 1:length(SubjectsAll)
    legStr{i} = sprintf('n = %i [%.1f]',i,PctCovs(i));
end
legend(legStr)
print(h,sprintf('mmp_Adj_nsubhist-%i',length(Subjects)),'-depsc')


h = figure;
rIdx = ~(sum((Adj_cov == 0)) == 180);
imagesc(Adj_cov(rIdx,rIdx))
cc = parula(10*64);
cc = [[1 1 1];cc];
colormap(cc);
fsize = 3;
set(gca,'xtick',[1:length(cLH.struct_names(rIdx))],'Ticklength',[2e-3 1e-5],...
    'xticklabel',cLH.struct_names(rIdx),'fontsize',fsize,'TickDir','out');
xtickangle(90);
set(gca,'ytick',[1:length(cLH.struct_names(rIdx))],'Ticklength',[2e-3 1e-5],...
    'yticklabel',cLH.struct_names(rIdx),'fontsize',fsize,'TickDir','out');
title(sprintf('HCP-MMP parcelation coverage (%.2f %%, n = %i)',pctCov,length(Subjects)),'fontsize',8)
c = colorbar; c.FontSize = 8;
print(h,sprintf('mmp_Adj_bipcov-%i',length(Subjects)),'-depsc')

h = figure;
rIdx = ~(sum((Adj_nsub == 0)) == 180);
imagesc(Adj_nsub(rIdx,rIdx))
cc = parula(10*64);
cc = [[1 1 1];cc];
colormap(cc);
fsize = 3;
set(gca,'xtick',[1:length(cLH.struct_names(rIdx))],'Ticklength',[2e-3 1e-5],...
    'xticklabel',cLH.struct_names(rIdx),'fontsize',fsize,'TickDir','out');
xtickangle(90);
set(gca,'ytick',[1:length(cLH.struct_names(rIdx))],'Ticklength',[2e-3 1e-5],...
    'yticklabel',cLH.struct_names(rIdx),'fontsize',fsize,'TickDir','out');
title(sprintf('HCP-MMP parcelation coverage (%.2f %%, n = %i)',pctCov,length(Subjects)),'fontsize',8)
c = colorbar; c.FontSize = 8;
print(h,sprintf('mmp_Adj_bipnsub-%i',length(Subjects)),'-depsc')
