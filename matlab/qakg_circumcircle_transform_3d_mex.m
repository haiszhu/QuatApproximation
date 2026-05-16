function [R, c, alpha] = qakg_circumcircle_transform_3d_mex(r_vert, R, c, alpha)
r_vert = double(reshape(r_vert, 3, 3));
R = double(reshape(R, 3, 3));
c = double(c(:));
alpha = double(alpha);
mex_id_ = 'qakg_circumcircle_transform_3d_mex(c i double[xx], c io double[xx], c io double[x], c io double[x])';
[R, c, alpha] = QuatApproximation_mex(mex_id_, r_vert, R, c, alpha, 3, 3, 3, 3, 3, 1);
end

