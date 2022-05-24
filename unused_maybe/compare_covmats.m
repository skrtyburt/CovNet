function [p_mat, z_obs_mat, z1_mat, z2_mat] = compare_covmats(mat1, mat2, n1, n2)

narginchk(4,4)

if ~isequal(size(mat1),size(mat2))
    fprinf(2,'Uequal matrix sizes.\n')
    return
end

[sz1, sz2] = size(mat1);
if sz1 ~= sz2
    fprintf(2,'diiferent number of row and column in matrices\n')
    return
end

for i = 1: sz1
    for j= 1: sz1
        [p_mat(i,j), z_obs_mat(i,j), z1_mat(i,j), z2_mat(i,j)] = fisher_r2z_compare(mat1(i,j), n1, mat2(i,j), n2);
    end
end



        