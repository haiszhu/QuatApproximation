function [omega_slp, omega] = qao_omegasdlp_mex(m, nterms, ncoeff, h_dim, r0, Ichi, Ialpha, omega_slp, omega)
m = double(m);
nterms = double(nterms);
ncoeff = double(ncoeff);
h_dim = double(h_dim);
nc4 = 4*ncoeff;
m4 = 4*m;
r0 = double(reshape(r0, 3, m));
Ichi = double(reshape(Ichi, m, nc4));
Ialpha = double(reshape(Ialpha, m, nc4));
if nargin < 8 || isempty(omega_slp), omega_slp = zeros(h_dim, m); end
if nargin < 9 || isempty(omega), omega = zeros(h_dim, m, 4); end
Ichi = reshape(Ichi, m, nc4);
Ialpha = reshape(Ialpha, m, nc4);
omega = reshape(omega, h_dim, m4);
mex_id_ = 'qao_omegasdlp_mex(c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i double[xx], c i dcomplex[xx], c i dcomplex[xx], c io double[xx], c io double[xx])';
[omega_slp, omega] = QuatApproximation_mex(mex_id_, m, nterms, ncoeff, h_dim, r0, Ichi, Ialpha, omega_slp, omega, 1, 1, 1, 1, 3, m, m, nc4, m, nc4, h_dim, m, h_dim, m4);
omega = reshape(omega, h_dim, m, 4);
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
