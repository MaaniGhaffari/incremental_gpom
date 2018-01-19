function R = R2d(theta)
% 2D rotation matrix
% Input angle is in radian.
R = [cos(theta) -sin(theta);
     sin(theta) cos(theta)];
end