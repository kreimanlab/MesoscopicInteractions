function [ bipChan ] = bip_chibi( sid )
% output columns: e1Chan e2Chan edist e1c(3) e2c(3)
%   e1Chan - channel number of first electrode
%   e2Chan - channel number of second electrode
%   edist - distance (mm) between electrodes
%   e1c - RAS coordinates of first electrode on pial surface (iELvis)
%   e2c - RAS coordinates of second electrode on pial surface (iELvis)

% [2018-02-19 19:19] Release 1.1 - Jerry Wang - jwang04@g.harvard.edu
%   * now takes grids of all dimensions
% [2018-03-23 23:30] bip_su is for macaque subject Su

%sid = 'sub4';
DEBUG = true;

if (~ strcmp(sid,'mChibi'))
    fprintf('[*] Only works with mChibi\n')
end

load('/Volumes/RawData/scripts/synth/mChibi/ChibiMAP.mat');
ecoord = [X, Y, zeros(size(X))];
ecoord = ecoord/(2*sqrt(2));

breaks = [4,13,22,33,44,51,55,62,64,79,91,104,114,118,120];
bipChan = [];
iB = 1;
for i = 1:(length(X)-1)
    if (~any(i == breaks))
        e1 = i;
        e2 = e1 + 1;
        e1c = ecoord(e1,:);
        e2c = ecoord(e2,:);
        edist = sqrt(sum((e1c-e2c).^2));
        bipChan(iB,:) = [e1 e2 edist e1c e2c];
        if (DEBUG)
            fprintf('> %i %i %.4f\n',e1,e2,edist)
        end
        iB = iB + 1;
    end
end



end
