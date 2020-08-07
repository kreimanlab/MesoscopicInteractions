close all;
clc;

ART_INDEX_MANUAL = 6;
% see: http://wiki.neurotycho.org/Anesthesia_Task_Details
% 20, 20, 30, 30, 15 minutes (total: 115)
cutoff_minutes_after_antagonist = 15;

msu_art = h5read('/media/klab/internal/data/h5_notch20/art_nosz/mSu_art.h5','/artifacts');
[~, n_art_segs] = size(msu_art);
se = round(linspace(1,n_art_segs,4));
sessions_start = se(1:3);
sessions_end = se(2:4)-1;
sessions_end(end) = n_art_segs;
n_session_segs = diff(sessions_start);
n_session_segs = round(n_session_segs(1));
for i = 1:length(sessions_start)
    start_idx = sessions_start(i) + (40/115) * n_session_segs;
    end_idx = sessions_start(i) + ((70 + cutoff_minutes_after_antagonist)/115) * n_session_segs;
    msu_art(:,round(start_idx):round(end_idx)) = ART_INDEX_MANUAL;
end

msu_art_idx = h5read('/media/klab/internal/data/h5_notch20/art_nosz/mSu_art.h5','/art_idx');
[~, n_art_segs] = size(msu_art_idx);
se = round(linspace(1,n_art_segs,4));
sessions_start = se(1:3);
sessions_end = se(2:4)-1;
sessions_end(end) = n_art_segs;
n_session_segs = diff(sessions_start);
n_session_segs = round(n_session_segs(1));
for i = 1:length(sessions_start)
    start_idx = sessions_start(i) + (40/115) * n_session_segs;
    end_idx = sessions_start(i) + ((70 + cutoff_minutes_after_antagonist)/115) * n_session_segs;
    msu_art_idx(:,round(start_idx):round(end_idx)) = ART_INDEX_MANUAL;
end

save('./cache/mSu_artifact_trim','msu_art','msu_art_idx');

% h5write('/media/klab/internal/data/h5_notch20/art/mSu_art.h5','/artifacts',msu_art);
% fprintf('[*] wrote to artifacts.\n')
% h5write('/media/klab/internal/data/h5_notch20/art/mSu_art.h5','/art_idx',msu_art_idx);
% fprintf('[*] wrote to art_idx.\n')


