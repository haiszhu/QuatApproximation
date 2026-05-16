function [sub_tris] = qaq_subdivide_simplex_mex(Cx, Cy, scale, orientation, sub_tris)
% One uniform-refinement step on a simplex triangle.
Cx          = double(Cx);
Cy          = double(Cy);
scale       = double(scale);
orientation = double(orientation);
if nargin < 5 || isempty(sub_tris), sub_tris = zeros(4, 4); end
sub_tris    = double(reshape(sub_tris, 4, 4));
mex_id_ = 'qaq_subdivide_simplex_mex(c i double[x], c i double[x], c i double[x], c i double[x], c io double[xx])';
[sub_tris] = QuatApproximation_mex(mex_id_, Cx, Cy, scale, orientation, sub_tris, 1, 1, 1, 1, 4, 4);
end

% --------------------------------------------------------------------------
