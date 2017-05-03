% central_diff.m
% Fahim Mannan (fmannan@cim.mcgill.ca)
function [dx, dy] = central_diff(u)
% ignore boundary
dx = u;
dy = u;

dx(:,2:end-1) = 0.5 * (u(:,3:end) - u(:,1:end-2));
dy(2:end-1,:) = 0.5 * (u(3:end,:) - u(1:end-2,:));
