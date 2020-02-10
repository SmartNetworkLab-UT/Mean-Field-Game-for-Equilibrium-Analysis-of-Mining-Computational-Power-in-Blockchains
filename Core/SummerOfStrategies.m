function Sum = SummerOfStrategies(ListOfAgents)
    Sum = 0;
    for i = 1:length(ListOfAgents)
        Sum = Sum + ListOfAgents(i).x;
    end
end