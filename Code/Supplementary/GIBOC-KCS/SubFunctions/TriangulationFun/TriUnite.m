function [ TR3 ] = TriUnite( TR1, TR2 )
%UNITE two unconnected triangulation objects , TR1, TR2

ConnectivityList = vertcat( TR1.ConnectivityList ,...
    TR2.ConnectivityList + length(TR1.Points));
Points = vertcat( TR1.Points , TR2.Points);

TR3 =  triangulation(ConnectivityList,Points);
end

