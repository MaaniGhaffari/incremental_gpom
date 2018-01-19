function [fp, yp] = beamSampling(laserScan)
% This function samples along scanned beams.
% Current implementation is a simple uniform sampling but adapted to the
% range.

nStates = length(laserScan);
yp = cell(size(laserScan));
fp = cell(size(laserScan));
Nmin = 1; % if laserScan is empty
for i = 1:nStates
    for j = 1:size(laserScan{i},2)
        Nmin = ceil(6.25 * sqrt(laserScan{i}(1,j).^2 + laserScan{i}(2,j).^2));
        Nmin = min(Nmin,15);
        fp{i}(:,j) = [laserScan{i}(1,j); laserScan{i}(2,j)] * 1/Nmin;
        yp{i}(j,1) = 1; % target
    end
    
    for m = 1:Nmin-2
        k = m * size(laserScan{i},2)+1:(m+1)*size(laserScan{i},2); % index tracker
        for j = 1:size(laserScan{i},2)
            fp{i}(:,k(j)) = [laserScan{i}(1,j); laserScan{i}(2,j)] * (m+1)/Nmin;
            yp{i}(k(j),1) = 1; % target
        end
    end
end