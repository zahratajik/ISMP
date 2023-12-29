function [repeated, missing] = findDuplicatesAndMissing(array)
    repeated = []; 
    missing = []; 

    for i = 1:length(array)
        if sum(array == array(i)) > 1 && ~ismember(array(i), repeated)
            repeated = [repeated, array(i)];
        end
        
        if ~ismember(i, array)
            missing = [missing, i];
        end
    end
end
