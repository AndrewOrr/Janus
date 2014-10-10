function [Graph] = RGER(Size,p)

    Graph = zeros(Size);

    for i = 1:Size
        for j = 1:Size
            if rand(1) <= p
                Graph(i,j) = 1;
            end
        end
    end

end