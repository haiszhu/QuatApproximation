function [uvs_all] = qaq_assemble_subdivided_nodes_mex(N, M, uvs_simplex, sub_tris, uvs_all)
% Lay out N reference simplex VR nodes inside each of M subdivided
% triangles described by sub_tris(M, 4).
N    = double(N);
M    = double(M);
NM   = N * M;
uvs_simplex = double(reshape(uvs_simplex, 2, N));
sub_tris    = double(reshape(sub_tris,    M, 4));
if nargin < 5 || isempty(uvs_all), uvs_all = zeros(2, NM); end
uvs_all     = double(reshape(uvs_all, 2, NM));
mex_id_ = 'qaq_assemble_subdivided_nodes_mex(c i int64_t[x], c i int64_t[x], c i double[xx], c i double[xx], c io double[xx])';
[uvs_all] = QuatApproximation_mex(mex_id_, N, M, uvs_simplex, sub_tris, uvs_all, 1, 1, 2, N, M, 4, 2, NM);
end

% --------------------------------------------------------------------------
