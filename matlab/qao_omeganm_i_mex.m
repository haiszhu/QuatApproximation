function [onm] = qao_omeganm_i_mex(r, dr, q_i, q_j, q_k)
r       = double(reshape(r,  3, []));
N       = size(r, 2);
dr      = double(reshape(dr, 3, N));
ncoeff  = size(q_i, 2);
nslots  = size(q_i, 3);
Ncoeff  = N * ncoeff;
q_i     = double(reshape(q_i, Ncoeff, nslots));
q_j     = double(reshape(q_j, Ncoeff, nslots));
q_k     = double(reshape(q_k, Ncoeff, nslots));
onm     = zeros(Ncoeff, nslots);
mex_id_ = 'qao_omeganm_i_mex(c i int64_t[x], c i int64_t[x], c i int64_t[x], c i double[xx], c i double[xx], c i double[xx], c i double[xx], c i double[xx], c io double[xx])';
[onm] = QuatApproximation_mex(mex_id_, N, ncoeff, nslots, r, dr, q_i, q_j, q_k, onm, 1, 1, 1, 3, N, 3, N, Ncoeff, nslots, Ncoeff, nslots, Ncoeff, nslots, Ncoeff, nslots);
onm     = reshape(onm, N, ncoeff, nslots);
end

% --------------------------------------------------------------------------
% qao_omegaall_mex
% Target-dependent Omega assembly across all four kernel coefficient
% channels (omega0..omega3).  Mirrors qotential's omegaall_mex.
%
% Inputs (MATLAB shapes match qotential's omegaall_mex.m signature):
%   r0   (3, m)        target points
%   M_all (n*morder, m) or (n, morder, m); reshaped internally to (dim1, m)
%   order (scalar)     -> morder = 2*order + 2
%   onm_0..onm_3 (n*h_dim, 4)  flattened (n, h_dim, 4)
%   h_dim (scalar)     -> h_dim = order^2
%   ijIdx (2, h_dim)
% Output:
%   Omega (m, 4*h_dim)
% --------------------------------------------------------------------------
