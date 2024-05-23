function Lp_2 = Fragmentation(repeated_1, PL, Lp_2, old_label, ind,qs)

for t = 1:length(repeated_1)
    
    % Find duplicate indexes
    w1(:,1) = find(PL(:, 1) == repeated_1(t));
    fragment(:,1) = PL(w1, 2);
    
    % Perform tagging for each duplicate element
    for j1 = 1:length(fragment)
        Lp_2(fragment(j1,1)).label = sprintf("F%d-%d", j1, (max(old_label) + length(ind) + qs));
        qs = qs + 1;
    end
    
    
    clear fragment w1 numbers
end

end
