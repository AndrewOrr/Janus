function [Area] = Gini(Original)
% clc
% Original = [1; 1; 1];
% c = [ones(length(Original),1) Original]

% L = 10;
% Original = zeros(L,1);
% for i = 1:L
%     Original(i) = f(i/L);
% end
% Original = sortrows(Original,-1)

    if sum(Original) == 0
        Area = 0;
    else   
        Normalized = zeros(length(Original)+1,2);
        Normalized(1,2) = sum(Original);
        Normalized(length(Normalized),1) = length(Normalized)-1;
        for i = 1:length(Original)
            Normalized(i+1,2) = Normalized(i,2)-Original(i);
            Normalized(i+1,1) = i;
        end
    %     Normalized;

        Normalized(:,2) = Normalized(:,2)/Normalized(1,2);
        Normalized(:,1) = Normalized(:,1)/Normalized(length(Normalized),1);

        Output = polyarea(Normalized(:,1),Normalized(:,2));

        Area = zeros(length(Normalized)-1,1);
        for j = 2:length(Normalized)
            Area(j-1) = (Normalized(j,2) + Normalized(j-1,2)) * 0.5 * (Normalized(j,1) - Normalized(j-1,1));
        end

        Area = sum(Area);
        Area = (0.5-Area)/0.5;
    end
end


% (0.5-trapz(Normalized(:,1),Normalized(:,2)))/0.5