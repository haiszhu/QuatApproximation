function [q_i, q_j, q_k] = qak_qnm_i_mex(idx, sx, lptype_id, F, F1, F2, F3, gradxyz)
idx       = double(idx);
lptype_id = double(lptype_id);
sx        = double(reshape(sx, 3, []));
N         = size(sx, 2);
F1        = double(reshape(F1, N, []));
ncoeff    = size(F1, 2);
F         = double(reshape(F,  N, ncoeff));
F2        = double(reshape(F2, N, ncoeff));
F3        = double(reshape(F3, N, ncoeff));
if nargin < 8 || isempty(gradxyz), gradxyz = zeros(3, N); end
gradxyz   = double(reshape(gradxyz, 3, N));
Ncoeff    = N * ncoeff;
q_i = zeros(Ncoeff, 5);
q_j = zeros(Ncoeff, 5);
q_k = zeros(Ncoeff, 5);
mex_id_ = 'qak_qnm_i_mex(c i int64_t[x], c i int64_t[x], c i int64_t[x], c i double[xx], c i int64_t[x], c i double[xx], c i double[xx], c i double[xx], c i double[xx], c i double[xx], c io double[xx], c io double[xx], c io double[xx])';
[q_i, q_j, q_k] = QuatApproximation_mex(mex_id_, idx, N, ncoeff, sx, lptype_id, F, F1, F2, F3, gradxyz, q_i, q_j, q_k, 1, 1, 1, 3, N, 1, N, ncoeff, N, ncoeff, N, ncoeff, N, ncoeff, 3, N, Ncoeff, 5, Ncoeff, 5, Ncoeff, 5);
q_i = reshape(q_i, N, ncoeff, 5);
q_j = reshape(q_j, N, ncoeff, 5);
q_k = reshape(q_k, N, ncoeff, 5);
end

% --------------------------------------------------------------------------
% qao_omeganm_i_mex
% omega^{(n,m)}_i boundary-integrand assembly: (dr x r) . q.
% nslots is inferred from the q_* inputs (size along dim 3).
% --------------------------------------------------------------------------
