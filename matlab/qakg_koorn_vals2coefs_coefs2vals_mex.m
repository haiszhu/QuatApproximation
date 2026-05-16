function [umatr, vmatr] = qakg_koorn_vals2coefs_coefs2vals_mex(korder, kpols, umatr, vmatr)
% Koornwinder vals<->coefs matrices at the VR nodes of order korder.
korder = double(korder);
kpols  = double(kpols);
if nargin < 3 || isempty(umatr), umatr = zeros(kpols, kpols); end
if nargin < 4 || isempty(vmatr), vmatr = zeros(kpols, kpols); end
umatr = double(reshape(umatr, kpols, kpols));
vmatr = double(reshape(vmatr, kpols, kpols));
mex_id_ = 'qakg_koorn_vals2coefs_coefs2vals_mex(c i int64_t[x], c i int64_t[x], c io double[xx], c io double[xx])';
[umatr, vmatr] = QuatApproximation_mex(mex_id_, korder, kpols, umatr, vmatr, 1, 1, kpols, kpols, kpols, kpols);
end

% --------------------------------------------------------------------------
