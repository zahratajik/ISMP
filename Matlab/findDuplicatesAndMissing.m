function [repeated, missing] = findDuplicatesAndMissing(array)
    repeated = []; % آرایه‌ای برای ذخیره اعداد تکراری
    missing = []; % آرایه‌ای برای ذخیره اعداد گم شده

    for i = 1:length(array)
        % بررسی تکراری بودن عدد
        if sum(array == array(i)) > 1 && ~ismember(array(i), repeated)
            repeated = [repeated, array(i)];
        end
        
        % بررسی عدد گم شده
        if ~ismember(i, array)
            missing = [missing, i];
        end
    end
end
