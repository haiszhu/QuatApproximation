function [x, nx, sp, rts, rps, rtts, rpps, rtps] = qaq_paraboloidparam_mex(npts, h, a, t, p, x, nx, sp, rts, rps, rtts, rpps, rtps)
% Downward-opening paraboloid eval: r(t,p) = [t; p; h*(1 - (t^2+p^2)/a^2)],
% with first/second partials and the (cross(rp, rt))-based outward normal.
npts = double(npts);
h    = double(h);
a    = double(a);
t    = double(t(:)');
p    = double(p(:)');
if nargin < 6  || isempty(x),    x    = zeros(3, npts); end
if nargin < 7  || isempty(nx),   nx   = zeros(3, npts); end
if nargin < 8  || isempty(sp),   sp   = zeros(1, npts); end
if nargin < 9  || isempty(rts),  rts  = zeros(3, npts); end
if nargin < 10 || isempty(rps),  rps  = zeros(3, npts); end
if nargin < 11 || isempty(rtts), rtts = zeros(3, npts); end
if nargin < 12 || isempty(rpps), rpps = zeros(3, npts); end
if nargin < 13 || isempty(rtps), rtps = zeros(3, npts); end
x    = double(reshape(x,    3, npts));
nx   = double(reshape(nx,   3, npts));
sp   = double(sp(:)');
rts  = double(reshape(rts,  3, npts));
rps  = double(reshape(rps,  3, npts));
rtts = double(reshape(rtts, 3, npts));
rpps = double(reshape(rpps, 3, npts));
rtps = double(reshape(rtps, 3, npts));
mex_id_ = 'qaq_paraboloidparam_mex(c i int64_t[x], c i double[x], c i double[x], c i double[x], c i double[x], c io double[xx], c io double[xx], c io double[x], c io double[xx], c io double[xx], c io double[xx], c io double[xx], c io double[xx])';
[x, nx, sp, rts, rps, rtts, rpps, rtps] = QuatApproximation_mex(mex_id_, npts, h, a, t, p, x, nx, sp, rts, rps, rtts, rpps, rtps, 1, 1, 1, npts, npts, 3, npts, 3, npts, npts, 3, npts, 3, npts, 3, npts, 3, npts, 3, npts);
end

% --------------------------------------------------------------------------
