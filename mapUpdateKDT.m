function [fmu, fcov, fmap, ids] = mapUpdateKDT(mu, s2, maps, tdata, currentPose)
% This function updates GP continuous occupancy maps by using BCM.

% Incremental Gaussian processes occupancy mapping using range-finder sensors:
% Ghaffari Jadidi, M., Valls Miro, J. & Dissanayake, G. Auton Robot (2017). https://doi.org/10.1007/s10514-017-9668-3
% Pre-print available on: https://arxiv.org/abs/1605.00335

%{  
    Copyright (C) 2018  Maani Ghaffari Jadidi
    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details. 
%}

b = maps.param.boundaries; % valid regions for map points

fmu = maps.Mu;  % mean map
fcov = maps.C;  % variance map
ids = [];       % updated indices

% query points in global coordinates
q = repmat(currentPose(1:2),size(tdata.t,1),1)' + R2d(currentPose(3)) * tdata.t';

for j = 1:size(tdata.t,1)
    % check if the current query point is inside the map
    xi = q(:,j);
    if (xi(1) > b(2) || xi(1) < b(1) || xi(2) > b(4) || xi(2) < b(3))
        continue;
    end
    
    % find the corresponding global index
    idx = knnsearch(maps.mdl, [xi(1), xi(2)]);
    ids = [ids; idx];
    
    % point-wise map fusion using BCM
    [muf, sf] = fusionBCM(mu(j), maps.Mu(idx), s2(j), maps.C(idx));
    
    % update the corresponding map mean and variance values
    fmu(idx) = muf;
    fcov(idx) = sf;
end

% occupancy probabilities using Logistic regression
% remove points with zero variance to avoid division by zero
idx = fcov ~= 0;
lambda = zeros(size(fcov));
lambda(idx) = min(fcov(idx))./fcov(idx); % bounded information (inverse variance)

% Logistic regression to get probablity map in [0, 1] range for each point
fmap = logisticReg(fmu, lambda, maps.param.gamma);
