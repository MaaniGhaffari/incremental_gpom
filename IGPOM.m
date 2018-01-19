function output = IGPOM(robotPose, laserScan, Parameters, gpom_store)
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

% test data is a grid centred at the sensor
x = -Parameters.gp.param.testAreaSize:Parameters.gp.param.testPointDist:Parameters.gp.param.testAreaSize;
[X,Y] = meshgrid(x,x);
tdata = [];
tdata.t = [X(:), Y(:)];
tdata.size = size(X);

% sample free area points along each beam
[fp, yp] = beamSampling(laserScan);

meanfunc = Parameters.gp.meanfunc;
covfunc = Parameters.gp.covfunc;
likfunc = Parameters.gp.likfunc;
inffunc = Parameters.gp.inffunc;

% Computing GPOM
x_free = cell(size(laserScan));
x_occupied = cell(size(laserScan));

numSteps = length(laserScan); % number of steps/scans
fusedMaps = gpom_store.maps;

if Parameters.plot
    figure; hold on
end

for i = 1:numSteps
    currentPose = [robotPose.x(i), robotPose.y(i), robotPose.h(i)];
    
    % training data
    x_occupied{i} = [laserScan{i}(1,:)', laserScan{i}(2,:)']; % occupied points returned by laser
    x_free{i} = [fp{i}(1,:)', fp{i}(2,:)']; % sampled points from free area along each laser beam
    X_training = [x_free{i}; x_occupied{i}]; % design matrix
    y_training = [-yp{i}(:); ones(size(x_occupied{i},1),1)]; % target values
    
    % GP regression
    if i == 1 && Parameters.gp.opthp
        gpom_store.hyp = minimize(gpom_store.hyp, @gp, Parameters.gp.iterMaxf, inffunc, meanfunc, covfunc, likfunc, X_training, y_training);        
    end
        
    [mu, s2] = gp(gpom_store.hyp, inffunc, meanfunc, covfunc, likfunc, X_training, y_training, tdata.t);
    
    % Updating the maps
    [fmu, fcov, fmap, ids] = mapUpdateKDT(mu, s2, fusedMaps, tdata, currentPose);
    
    fusedMaps.P = fmap;
    fusedMaps.C = fcov;
    fusedMaps.Mu = fmu;
    fusedMaps.ids = ids;
    
    if Parameters.plot
        pcolor(Parameters.gp.testPoint.X, Parameters.gp.testPoint.Y, reshape(fusedMaps.P,fusedMaps.size)), caxis([0, 1])
        colormap jet; colorbar, shading interp, grid off
        axis equal, axis(Parameters.plotArea)
        drawnow
    end

end

output = gpom_store;
output.maps = fusedMaps;

end
