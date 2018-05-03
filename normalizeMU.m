function [A B C] = normalizeMU(MU0, MU1, MU2)

[n m] =size(MU0);

for i=1:n
    for j=1:m
        vect = [MU0(i,j) MU1(i,j) MU2(i,j)];
        norm_vect = vect/max(vect);
        A(i,j) = norm_vect(1);
        B(i,j) = norm_vect(2);
        C(i,j) = norm_vect(3);
        %input('')
    end
end


end