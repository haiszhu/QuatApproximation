function [fx, fy, fz, f, ijidx] = qah_evaltensorproductharmonicgrad_mex(nt, r, order, fx, fy, fz, f, ijidx)
% Tensor-product harmonic polynomial basis H_{i,j} and Cartesian gradient
% evaluated at nt 3D points r(3, nt). Outputs are (nt, order^2) real
% arrays; ijidx is (2, order^2) integer pairs.
nt     = double(nt);
order  = double(order);
order2 = order*order;
r      = double(reshape(r, 3, nt));
if nargin < 4 || isempty(fx),    fx    = zeros(nt, order2); end
if nargin < 5 || isempty(fy),    fy    = zeros(nt, order2); end
if nargin < 6 || isempty(fz),    fz    = zeros(nt, order2); end
if nargin < 7 || isempty(f),     f     = zeros(nt, order2); end
if nargin < 8 || isempty(ijidx), ijidx = zeros(2, order2);  end
fx    = double(reshape(fx,    nt, order2));
fy    = double(reshape(fy,    nt, order2));
fz    = double(reshape(fz,    nt, order2));
f     = double(reshape(f,     nt, order2));
ijidx = double(reshape(ijidx, 2,  order2));
mex_id_ = 'qah_evaltensorproductharmonicgrad_mex(c i int64_t[x], c i double[xx], c i int64_t[x], c io double[xx], c io double[xx], c io double[xx], c io double[xx], c io int64_t[xx])';
[fx, fy, fz, f, ijidx] = QuatApproximation_mex(mex_id_, nt, r, order, fx, fy, fz, f, ijidx, 1, 3, nt, 1, nt, order2, nt, order2, nt, order2, nt, order2, 2, order2);
end

% --------------------------------------------------------------------------
