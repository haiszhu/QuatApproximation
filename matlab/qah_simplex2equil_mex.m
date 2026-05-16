function [uvs_eq] = qah_simplex2equil_mex(npts, h_side, uvs, uvs_eq)
% Map npts simplex points uvs(2,npts) to an equilateral triangle of
% side length h_side centered at the origin.
npts    = double(npts);
h_side  = double(h_side);
uvs     = double(reshape(uvs, 2, npts));
if nargin < 4 || isempty(uvs_eq), uvs_eq = zeros(2, npts); end
uvs_eq  = double(reshape(uvs_eq, 2, npts));
mex_id_ = 'qah_simplex2equil_mex(c i int64_t[x], c i double[x], c i double[xx], c io double[xx])';
[uvs_eq] = QuatApproximation_mex(mex_id_, npts, h_side, uvs, uvs_eq, 1, 1, 2, npts, 2, npts);
end

% --------------------------------------------------------------------------
