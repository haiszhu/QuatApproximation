function [F, Fx, Fy, Fz, ier] = qah_l3dtavecevalmat_mex(ztargs, nt, nterms, F, Fx, Fy, Fz, ier)
% 3D Laplace solid spherical harmonics: values F and Cartesian gradient
% (Fx, Fy, Fz) at target points ztargs(3, nt). Each output is shaped
% (nt, (nterms+1)^2), one column per basis index (l, m) for
% 0 <= l <= nterms, -l <= m <= l.
%
% Thin r64 wrapper for harmonic_mod::l3dtavecevalmat_r64.
nt     = double(nt);
nterms = double(nterms);
nc     = (nterms+1)*(nterms+1);
ztargs = double(reshape(ztargs, 3, nt));
if nargin < 4 || isempty(F),   F   = complex(zeros(nt, nc)); end
if nargin < 5 || isempty(Fx),  Fx  = complex(zeros(nt, nc)); end
if nargin < 6 || isempty(Fy),  Fy  = complex(zeros(nt, nc)); end
if nargin < 7 || isempty(Fz),  Fz  = complex(zeros(nt, nc)); end
if nargin < 8 || isempty(ier), ier = 0;                       end
F   = complex(reshape(F,  nt, nc));
Fx  = complex(reshape(Fx, nt, nc));
Fy  = complex(reshape(Fy, nt, nc));
Fz  = complex(reshape(Fz, nt, nc));
ier = double(ier);
mex_id_ = 'qah_l3dtavecevalmat_mex(c i double[xx], c i int64_t[x], c i int64_t[x], c io dcomplex[xx], c io dcomplex[xx], c io dcomplex[xx], c io dcomplex[xx], c io int64_t[x])';
[F, Fx, Fy, Fz, ier] = QuatApproximation_mex(mex_id_, ztargs, nt, nterms, F, Fx, Fy, Fz, ier, 3, nt, 1, 1, nt, nc, nt, nc, nt, nc, nt, nc, 1);
end

% --------------------------------------------------------------------------
