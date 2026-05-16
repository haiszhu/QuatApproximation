function [Omega] = qao_omegaall_mex(r0, M_all, order, onm_0, onm_1, onm_2, onm_3, h_dim, ijIdx)
order   = double(order);
h_dim   = double(h_dim);
r0      = double(reshape(r0, 3, []));
m       = size(r0, 2);
M_all   = double(reshape(M_all, [], m));
dim1    = size(M_all, 1);
n       = size(onm_0, 1) / h_dim;
morder  = 2*order + 2;
h_dim4  = 4 * h_dim;
nh_dim  = n * h_dim;
onm_0   = double(reshape(onm_0, nh_dim, 4));
onm_1   = double(reshape(onm_1, nh_dim, 4));
onm_2   = double(reshape(onm_2, nh_dim, 4));
onm_3   = double(reshape(onm_3, nh_dim, 4));
ijIdx   = double(reshape(ijIdx, 2, h_dim));
Omega   = zeros(m, h_dim4);
mex_id_ = 'qao_omegaall_mex(c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i double[xx], c i double[xx], c i double[xx], c i double[xx], c i double[xx], c i double[xx], c i int64_t[xx], c io double[xx])';
[Omega] = QuatApproximation_mex(mex_id_, m, dim1, n, h_dim, morder, r0, M_all, onm_0, onm_1, onm_2, onm_3, ijIdx, Omega, 1, 1, 1, 1, 1, 3, m, dim1, m, nh_dim, 4, nh_dim, 4, nh_dim, 4, nh_dim, 4, 2, h_dim, m, h_dim4);
end

% --------------------------------------------------------------------------
% qatg_line3quadr_3dline_T_mex
% Tensor-product analog of qakg_line3quadr_3dline_mex: GL boundary nodes,
% weights, tangents, speeds for a square source patch.
%
% Inputs (MATLAB shapes match the qotential 'T' branch of quadr_3dline):
%   x_uvs (3, ordert^2)  source patch values at GL tensor grid
%   ordert (scalar)      1D tensor order (inferred from sqrt(size(x_uvs,2)))
%   nquad  (scalar)      panel order
%   tgl, wgl (nquad)     GL nodes/weights
%   Dgl    (nquad,nquad) GL differentiation matrix
%   tpan   (sbdnp+1)     panel breakpoints in [0, 2pi]
% Outputs:
%   sxbd   (3, nbd)      boundary node positions
%   swbd   (nbd)         speed-weighted GL weights
%   stangbd(3, nbd)      unit tangent vectors
%   sspbd  (nbd)         speeds
% --------------------------------------------------------------------------
