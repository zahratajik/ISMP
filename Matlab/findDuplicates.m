function [repeated] = findDuplicates(array)
    repeated = []; 
    for i = 1:length(array)
        if sum(array == array(i)) > 1 && ~ismember(array(i), repeated)
            repeated = [repeated, array(i)];
        end
    end
end
