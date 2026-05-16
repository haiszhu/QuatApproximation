function [sxbd, swbd, stangbd, sspbd, r_vert] = qakg_line3quadr_3dline_mex(x_uvs, korder, kpols, umatr, nquad, tgl, wgl, Dgl, sbdnp, tpan, nbd, sxbd, swbd, stangbd, sspbd, r_vert)
% GL boundary quadrature (nodes/weights/tangents/speeds) along the three
% edges of a curved triangle described by Koornwinder coeffs.
korder = double(korder);
kpols  = double(kpols);
nquad  = double(nquad);
sbdnp    = double(sbdnp);
nbd      = double(nbd);
sbdnpp1  = sbdnp + 1;
x_uvs    = double(reshape(x_uvs, 3, kpols));
umatr    = double(reshape(umatr, kpols, kpols));
tgl      = double(tgl(:));
wgl      = double(wgl(:));
Dgl      = double(reshape(Dgl, nquad, nquad));
tpan     = double(tpan(:));
if nargin < 12 || isempty(sxbd),    sxbd    = zeros(3, nbd); end
if nargin < 13 || isempty(swbd),    swbd    = zeros(1, nbd); end
if nargin < 14 || isempty(stangbd), stangbd = zeros(3, nbd); end
if nargin < 15 || isempty(sspbd),   sspbd   = zeros(1, nbd); end
if nargin < 16 || isempty(r_vert),  r_vert  = zeros(3, 3);   end
sxbd    = double(reshape(sxbd,    3, nbd));
swbd    = double(swbd(:)');
stangbd = double(reshape(stangbd, 3, nbd));
sspbd   = double(sspbd(:)');
r_vert  = double(reshape(r_vert, 3, 3));
mex_id_ = 'qakg_line3quadr_3dline_mex(c i double[xx], c i int64_t[x], c i int64_t[x], c i double[xx], c i int64_t[x], c i double[x], c i double[x], c i double[xx], c i int64_t[x], c i double[x], c i int64_t[x], c io double[xx], c io double[x], c io double[xx], c io double[x], c io double[xx])';
[sxbd, swbd, stangbd, sspbd, r_vert] = QuatApproximation_mex(mex_id_, x_uvs, korder, kpols, umatr, nquad, tgl, wgl, Dgl, sbdnp, tpan, nbd, sxbd, swbd, stangbd, sspbd, r_vert, 3, kpols, 1, 1, kpols, kpols, 1, nquad, nquad, nquad, nquad, 1, sbdnpp1, 1, 3, nbd, nbd, 3, nbd, nbd, 3, 3);
end

% ============================================================
% harmonic_mod wrappers (qah_*)
% ============================================================

% ============================================================
% quatapproximation_mod wrappers (qaq_*)
% ============================================================

% --------------------------------------------------------------------------
