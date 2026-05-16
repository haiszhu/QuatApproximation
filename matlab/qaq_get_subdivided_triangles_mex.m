function [all_tris] = qaq_get_subdivided_triangles_mex(C0, s0, o0, L, all_tris)
% Recursively apply subdivide_simplex L times; output is (4^L, 4).
C0 = double(C0(:));
s0 = double(s0);
o0 = double(o0);
L  = double(L);
nL = 4^L;
if nargin < 5 || isempty(all_tris), all_tris = zeros(nL, 4); end
all_tris = double(reshape(all_tris, nL, 4));
mex_id_ = 'qaq_get_subdivided_triangles_mex(c i double[x], c i double[x], c i double[x], c i int64_t[x], c io double[xx])';
[all_tris] = QuatApproximation_mex(mex_id_, C0, s0, o0, L, all_tris, 2, 1, 1, 1, nL, 4);
end

% --------------------------------------------------------------------------
