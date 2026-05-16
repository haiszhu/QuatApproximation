function [wts] = qakg_get_vioreanu_wts_mex(norder, npols, wts)
% Vioreanu-Rokhlin reference-triangle quadrature weights of order norder.
norder = double(norder);
npols  = double(npols);
if nargin < 3 || isempty(wts), wts = zeros(npols, 1); end
wts = double(wts(:));
mex_id_ = 'qakg_get_vioreanu_wts_mex(c i int64_t[x], c i int64_t[x], c io double[x])';
[wts] = QuatApproximation_mex(mex_id_, norder, npols, wts, 1, 1, npols);
end

% --------------------------------------------------------------------------
