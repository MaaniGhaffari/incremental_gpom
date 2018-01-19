function p = logisticReg(mu,lambda, gamma)
% This function given all parameters computes probabilities using a Logistic
% regression classifier
p = 1 ./ (1 + exp(-gamma * mu .* sqrt(lambda)));
