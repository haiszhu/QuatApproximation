function [F, Fx, Fy, Fz, ier] = qah_l3dtavecevalmat_r128_mex(ztargs, nt, nterms, F, Fx, Fy, Fz, ier, flag)
% r128 evaluation with the design-spec §10 trailing-flag contract.
%   flag = 0 : cast r64 inputs to r128, run harmonic_mod::l3dtavecevalmat_r128,
%              cast outputs back to r64.
%   flag = 1 : reserved for HDF5 round-trip with a downstream r128 consumer
%              (layout deferred; currently aborts).
nt     = double(nt);
nterms = double(nterms);
nc     = (nterms+1)*(nterms+1);
ztargs = double(reshape(ztargs, 3, nt));
if nargin < 4 || isempty(F),    F    = complex(zeros(nt, nc)); end
if nargin < 5 || isempty(Fx),   Fx   = complex(zeros(nt, nc)); end
if nargin < 6 || isempty(Fy),   Fy   = complex(zeros(nt, nc)); end
if nargin < 7 || isempty(Fz),   Fz   = complex(zeros(nt, nc)); end
if nargin < 8 || isempty(ier),  ier  = 0;                       end
if nargin < 9 || isempty(flag), flag = 0;                       end
F    = complex(reshape(F,  nt, nc));
Fx   = complex(reshape(Fx, nt, nc));
Fy   = complex(reshape(Fy, nt, nc));
Fz   = complex(reshape(Fz, nt, nc));
ier  = double(ier);
flag = double(flag);
mex_id_ = 'qah_l3dtavecevalmat_r128_mex(c i double[xx], c i int64_t[x], c i int64_t[x], c io dcomplex[xx], c io dcomplex[xx], c io dcomplex[xx], c io dcomplex[xx], c io int64_t[x], c i int64_t[x])';
[F, Fx, Fy, Fz, ier] = QuatApproximation_mex(mex_id_, ztargs, nt, nterms, F, Fx, Fy, Fz, ier, flag, 3, nt, 1, 1, nt, nc, nt, nc, nt, nc, nt, nc, 1, 1);
end

% --------------------------------------------------------------------------
% qak_qnm_i_mex
% q^{(n,m)}_i kernel-coefficient builder. idx in {0,1,2,3}.
% lptype_id: 0='d' (DLP/SLPn, 4 active slots), 1='s' (SLP, 5 slots),
%            2='T' (DLPn, NOT YET IMPLEMENTED -- error_stop).
% Outputs are sized at the max rank (5 slots); inactive slots are zero.
% --------------------------------------------------------------------------
