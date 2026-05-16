function [w_bclag] = qaq_bclaginterpweights_mex(n, tgl, w_bclag)
% Barycentric Lagrange interp weights at the given nodes tgl.
n   = double(n);
tgl = double(tgl(:));
if nargin < 3 || isempty(w_bclag), w_bclag = zeros(n, 1); end
w_bclag = double(w_bclag(:));
mex_id_ = 'qaq_bclaginterpweights_mex(c i int64_t[x], c i double[x], c io double[x])';
[w_bclag] = QuatApproximation_mex(mex_id_, n, tgl, w_bclag, 1, n, n);
end

% --------------------------------------------------------------------------
