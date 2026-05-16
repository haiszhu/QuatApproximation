function [x] = qaq_paraboloid_uv_mex(npts, H, u, v, x)
% Elliptic paraboloid x(3,npts) = [u; v; (H/2)(u^2+v^2)].
npts = double(npts);
H    = double(H);
u    = double(u(:)');
v    = double(v(:)');
if nargin < 5 || isempty(x), x = zeros(3, npts); end
x    = double(reshape(x, 3, npts));
mex_id_ = 'qaq_paraboloid_uv_mex(c i int64_t[x], c i double[x], c i double[x], c i double[x], c io double[xx])';
[x] = QuatApproximation_mex(mex_id_, npts, H, u, v, x, 1, 1, npts, npts, 3, npts);
end

% --------------------------------------------------------------------------
