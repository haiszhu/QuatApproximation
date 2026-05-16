function [tgl, wgl, Dgl] = qaq_gauss_mex(n, tgl, wgl, Dgl)
% Gauss-Legendre nodes (tgl), weights (wgl), and differentiation matrix (Dgl)
% on [-1, 1] for the given order n.
n = double(n);
if nargin < 2 || isempty(tgl), tgl = zeros(n, 1);    end
if nargin < 3 || isempty(wgl), wgl = zeros(n, 1);    end
if nargin < 4 || isempty(Dgl), Dgl = zeros(n, n);    end
tgl = double(tgl(:));
wgl = double(wgl(:));
Dgl = double(reshape(Dgl, n, n));
mex_id_ = 'qaq_gauss_mex(c i int64_t[x], c io double[x], c io double[x], c io double[xx])';
[tgl, wgl, Dgl] = QuatApproximation_mex(mex_id_, n, tgl, wgl, Dgl, 1, n, n, n, n);
end

% --------------------------------------------------------------------------
