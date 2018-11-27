function [indic] = zerone(a,row,first,last);
indic = 0;
for j = first:last
    rj = a{row}(j);
    if rj == '1'
        indic = j;
    end
end