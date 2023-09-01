function [ TRout ] = TriReduceMesh( TR, ElmtsKept, NodesKept )
%Remove unnecessary Node and Renumber the elements accordingly 
%   Detailed explanation goes here
% TR : A triangulation object with n Elements
% ElmtsKept : A nx1 vectors of index of the rows of kept elments, or a
% binary indicating kept Elements
% NodesKept : ID of kept nodes of corresponding to TR connectibity list OR 
% list of nodes coordinates

if nargin >2
    
    if sum(sum(mod(NodesKept,1))) == 0 % NodesID given
        NodesIDKept = NodesKept;
    else % Nodes Coordinates given
        NodesIDKept = TR.nearestNeighbor(NodesKept);
        NodesIDKept = unique(NodesIDKept);
    end
    
    ElmtsKept = TR.vertexAttachments(NodesIDKept);
    ElmtsKept = unique(transpose(horzcat(ElmtsKept{:})));
end
     
    ElmtsNodesIDKept = TR.ConnectivityList(ElmtsKept,:);
    NodesKept = unique(ElmtsNodesIDKept(:));
    PointsKept = TR.Points(NodesKept,:);
    
    IndexTrsfrm = zeros(length(TR.Points),1);
    IndexTrsfrm(NodesKept) = 1 : length(NodesKept);
    
    TRout = triangulation(IndexTrsfrm(ElmtsNodesIDKept),PointsKept);
end

