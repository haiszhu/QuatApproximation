function [uvs_eq, detJ] = qaq_simplex2equil_with_detj_mex(npts, h_side, uvs, uvs_eq, detJ)
% simplex2equil with Jacobian determinant (detJ = sqrt(3)/2 * h_side^2).
npts    = double(npts);
h_side  = double(h_side);
uvs     = double(reshape(uvs, 2, npts));
if nargin < 4 || isempty(uvs_eq), uvs_eq = zeros(2, npts); end
if nargin < 5 || isempty(detJ),   detJ   = 0;              end
uvs_eq  = double(reshape(uvs_eq, 2, npts));
detJ    = double(detJ);
mex_id_ = 'qaq_simplex2equil_with_detj_mex(c i int64_t[x], c i double[x], c i double[xx], c io double[xx], c io double[x])';
[uvs_eq, detJ] = QuatApproximation_mex(mex_id_, npts, h_side, uvs, uvs_eq, detJ, 1, 1, 2, npts, 2, npts, 1);
end

% --------------------------------------------------------------------------
