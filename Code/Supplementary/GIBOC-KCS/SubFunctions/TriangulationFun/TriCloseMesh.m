function [ TRout ] = TriCloseMesh( TRsup , TRin , nbElmts )
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here
[ TR ] = TriDilateMesh( TRsup, TRin, nbElmts );
[ TRout ] = TriErodeMesh( TR, nbElmts );

end

