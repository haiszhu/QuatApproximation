function [uvs] = qakg_get_vioreanu_nodes_mex(norder, npols, uvs)
% Vioreanu-Rokhlin reference-triangle nodes of order norder
% (npols = (norder+1)*(norder+2)/2 columns).
norder = double(norder);
npols  = double(npols);
if nargin < 3 || isempty(uvs), uvs = zeros(2, npols); end
uvs = double(reshape(uvs, 2, npols));
mex_id_ = 'qakg_get_vioreanu_nodes_mex(c i int64_t[x], c i int64_t[x], c io double[xx])';
[uvs] = QuatApproximation_mex(mex_id_, norder, npols, uvs, 1, 1, 2, npols);
end

% --------------------------------------------------------------------------
