function Lp_2 = Merge(repeated_2, PL, Lp_1, Lp_2, n)
    if n==2
        for t = 1:length(repeated_2)
            % Find duplicate indexes
            ww(:,1) = find(PL(:, 2) == repeated_2(t));
            merge(:,1) = PL(ww, 1);
            
            % Calculate the area to find the maximum value
            Area = zeros(length(merge), 1);
            for j2 = 1:length(merge)
                Area(j2) = Lp_1(merge(j2)).area;
            end
            
            w_max = find(Area == max(Area)); %biggest Area
            
            % Perform tagging for each duplicate element
            Lp_2(repeated_2(t)).label = "M-" + num2str(Lp_1(merge(w_max(1))).label);
            
            clear w_max ww merge Area
        end
    else
        for t=1:length(repeated_2)
            ww(:,1) = find( PL(:,2) == repeated_2(t)); % Find duplicate numbers in step new
            merge(:,1) = PL(ww,1); % equivalent to the counters in step old
            for j2=1:length(merge(:,1))
                Area(j2,1)= double(Lp_1(merge(j2,1)).area); % find Area of each patch
            end
            w_max = find( Area(:,1) == max(Area(:,1))); % maximum Area
            
            numbers = regexp ((Lp_1(merge(w_max(1))).label), '[0-9]+', 'match');
            
            if length(numbers(1,:))==1
                Lp_2(repeated_2(t)).label = ("M-" + "" + numbers);
            else
                
                Lp_2(repeated_2(t)).label = ("M-"+""+numbers(1,2));
            end
            clear w_max ww merge Area
        end
    end
end
