function [ bipChan ] = bip_su( sid )
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
DEBUG = false;

if (~ strcmp(sid,'mSu'))
    fprintf('[*] Only works with mSu\n')
end

load([sid,'/SuMAP.mat']);
ecoord = [Su.X, Su.Y, zeros(size(Su.X))];
ecoord = ecoord/(2*sqrt(2));

breaks = [8,16,25,35,46,55,64,73,82,86,90,94,106,121,125];
bipChan = [];
iB = 1;
for i = 1:(length(Su.X)-1)
    if (sum(i == breaks) == 0)
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
