function [R, c, alpha] = qakg_circumcircle_transform_3d_r128_mex(r_vert, R, c, alpha, ktri, flag)
r_vert = double(reshape(r_vert, 3, 3));
if nargin < 2 || isempty(R),     R     = zeros(3, 3); end
if nargin < 3 || isempty(c),     c     = zeros(3, 1); end
if nargin < 4 || isempty(alpha), alpha = 0;           end
if nargin < 5 || isempty(ktri),  ktri  = 0;           end
if nargin < 6 || isempty(flag),  flag  = 0;           end
R     = double(reshape(R, 3, 3));
c     = double(c(:));
alpha = double(alpha);
ktri  = double(ktri);
flag  = double(flag);
mex_id_ = 'qakg_circumcircle_transform_3d_r128_mex(c i double[xx], c io double[xx], c io double[x], c io double[x], c i int64_t[x], c i int64_t[x])';
[R, c, alpha] = QuatApproximation_mex(mex_id_, r_vert, R, c, alpha, ktri, flag, 3, 3, 3, 3, 3, 1, 1, 1);
end

