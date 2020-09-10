function [ bipChan ] = bip( sid )
% output columns: e1Chan e2Chan edist e1c(3) e2c(3)
%   e1Chan - channel number of first electrode
%   e2Chan - channel number of second electrode
%   edist - distance (mm) between electrodes
%   e1c - RAS coordinates of first electrode on pial surface (iELvis)
%   e2c - RAS coordinates of second electrode on pial surface (iELvis)

% [2018-02-19 19:19] Release 1.1 - Jerry Wang - jwang04@g.harvard.edu
%   * now takes grids of all dimensions
% [2018-02-04 00:13] Release 1.0 - Jerry Wang - jwang04@g.harvard.edu

%sid = 'sub4';

if (strcmp(sid,'mSu'))
    bipChan = bip_su(sid);
    return
end

if strcmp(sid, 'mChibi')
    bipChan = bip_chibi(sid);
    return;
end

% Debug
SHOW_RESULT = true;
DEBUG = false;

% Constant
WARN_THRESHOLD_MM = 20;

% Load mgrid
mgfn = sprintf('%s/elec_recon/%s.mgrid',sid,sid);
MG = textscan(fopen(mgfn,'r'),'%s');
MG = MG{1};

dim_idx = find(contains(MG,'Dimensions','IgnoreCase',true));
Gnames = cell(size(dim_idx));
Gcols = zeros(size(dim_idx));
Grows = zeros(size(dim_idx));
for i = 1:length(dim_idx)
    Gnames{i} = MG{dim_idx(i)-1};
    col = strsplit(MG{dim_idx(i)+1});
    row = strsplit(MG{dim_idx(i)+2});
    Gcols(i) = str2double(col{1});
    Grows(i) = str2double(row{1});
end

Gnames4sort = cell(size(Gnames));
for i2 = 1:length(Gnames4sort)
    st = strsplit(Gnames{i2},'_');
    Gnames4sort{i2} = st{2};
end
[~,sortI] = sort(Gnames4sort);
Gnames = Gnames(sortI);
Gcols = Gcols(sortI);
Grows = Grows(sortI);

% Load electrode names
infn = sprintf('%s/elec_recon/%s.electrodeNames',sid,sid);
D = textscan(fopen(infn,'r'),'%s %s %s');

% Load ielvis coordinates
infn2 = sprintf('%s/elec_recon/%s.PIAL',sid,sid);
D2 = textscan(fopen(infn2,'r'),'%s %s %s');

% Load ielvis labels
elecLabel = read_label(sprintf('%',sid),sprintf(...
        './%s/label/all_surf_ielvis.label',sid));

% Process 
n_header_lines = 2;
eNames = D{1}((n_header_lines+1):end);
eGridName = cell(size(eNames));
eGridNum = nan(size(eNames));
for i = 1:length(eNames)
    isL = isletter(eNames{i});
    eGridName{i} = eNames{i}(isL);
    eGridNum(i) = str2double(eNames{i}(~isL));
end
gridNames = sort(unique(eGridName));
%gridNames = gridNames(sortI);

if (DEBUG)
    for i = 1:length(gridNames)
        fprintf('[%s] %s: %s (%i X %i)\n',sid,Gnames{i},gridNames{i},Gcols(i),Grows(i))
    end
end

% Process coordinates
n_chan = length(D2{1})-n_header_lines;
elvisCoord = nan(n_chan,3);
for i = 1:3
    cTmp = nan(n_chan,1);
    for j = 1:n_chan
        cTmp(j) = str2double(D2{i}{n_header_lines+j});
    end
    elvisCoord(:,i) = cTmp;
    %elvisCoord(:,i) = D2{i}((n_header_lines+1):end);
end

n_warn = 0;
bipChan = [];
bcIdx = 1;
for i = 1:length(gridNames)
    isGrid = strcmp(eGridName,gridNames{i});
    gName = eGridName(isGrid);
    gNum = eGridNum(isGrid);
    for j = 1:length(gName)
        if (mod(j,Grows(i)) ~= 0)
            %fprintf('---\n')
            e1 = sprintf('%s%i',gridNames{i},j);
            e2 = sprintf('%s%i',gridNames{i},j+1);
            %fprintf('e1: %s, e2: %s\n',e1,e2)
            e1Idx = find(strcmp(eNames,e1));
            e2Idx = find(strcmp(eNames,e2));
            %fprintf('e1Idx: %i, e2Idx: %i\n',e1Idx,e2Idx)
            e1Chan = elecLabel(e1Idx,end);
            e2Chan = elecLabel(e2Idx,end);
            %fprintf('e1Chan: %i, e2Chan: %i\n',e1Chan,e2Chan)
            e1c = elvisCoord(e1Idx,:);%elecLabel(e1Idx,2:4);
            e2c = elvisCoord(e2Idx,:);%elecLabel(e2Idx,2:4);
            edist = sqrt(sum((e1c - e2c).^2));
            if (edist < WARN_THRESHOLD_MM)
                %fprintf('%s[%i]-%s[%i] %.2f mm\n',e1,e1Chan,e2,e2Chan,edist);
            else
                fprintf(2,'[*] %s[%i]-%s[%i] %.2f mm\n',e1,e1Chan,e2,e2Chan,edist);
                n_warn = n_warn + 1;
            end

            if (DEBUG)
                fprintf('\t%s[%i]-%s[%i] %.2f mm\n',e1,e1Chan,e2,e2Chan,edist);
            end
            bipChan(bcIdx,:) = [e1Chan e2Chan edist e1c e2c];
            bcIdx = bcIdx + 1;
        end
    end
end

if (SHOW_RESULT)
    fprintf((1 + (n_warn>0)),...
        '[%s] BIP montage %i of %i channels exceed %.2f mm.\n',...
        sid,n_warn,bcIdx-1,WARN_THRESHOLD_MM);
end


end
