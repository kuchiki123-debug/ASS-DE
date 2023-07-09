function vi=boundConstraint(vi, lu);

[NP, D] = size(vi);  % the population size and the problem's dimension

for i= 1:NP
    for j=1:D
        if vi(i,j)<lu(1,j)||vi(i,j)>lu(2,j)
            vi(i,j) = lu(1,j)+(lu(2,j)-lu(1,j))*rand;
        end
    end
end

vi(:,5) = round(vi(:,5));
vi(:,6) = round(vi(:,6));