% Fahim Mannan (fmannan@cim.mcgill.ca)
% compute 2nd derivative using central difference
% ignore the boundary in the computation
% call the function with appropriate padding
function [dxx, dyy, dxy] = central_diff2(u, padtype)
if(~exist('padtype', 'var'))
    padtype = 'replicate';
end
% initialize the border cells
[dxx, dyy] = central_diff(u);
[t, dxy1] = central_diff(padarray(dxx, [1,1], padtype));
[dyx1, t] = central_diff(padarray(dyy, [1,1], padtype));
dxy = 0.5 * (dxy1(2:end-1,2:end-1) + dyx1(2:end-1,2:end-1));

% compute using 2nd central difference
% (http://en.wikipedia.org/wiki/Finite_difference)
dxx(:,2:end-1) = u(:,3:end) + u(:,1:end-2) - 2*u(:,2:end-1);
dyy(2:end-1,:) = u(3:end,:) + u(1:end-2,:) - 2*u(2:end-1,:);
dxy(2:end-1,2:end-1) = 0.25 * -(u(1:end-2,3:end) - u(3:end,3:end) ...
                                - u(1:end-2,1:end-2) + u(3:end, 1:end-2));
                            