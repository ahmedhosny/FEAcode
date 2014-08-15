
%clean up
clear all; clc;
%Define File Name
% truss1, truss2, truss3 
% twoD1, twoD1a, twoD2, twoD3
% twoDQ1, twoDQ2, twoDQ3, twoDQ4, twoDQ5
NameTextFile='twoDQ2.inp'; 
%%Function to Read all files
[ElemConnectivityBar,ElemConnectivityTriPS,ElemConnectivityTriPE,ElemConnectivityQuadPS,ElemConnectivityQuadPE,myDistributedLoad,NodeForceBC,...
    NodeCoor,A,EnNU,NodeDispBC]=Readinput(NameTextFile);

% IF Bar and quad are empty,  run TRI element function
if isempty(ElemConnectivityBar) && isempty(ElemConnectivityQuadPS) && isempty(ElemConnectivityQuadPE)
    [DF,RE]=ContTri2d(ElemConnectivityTriPS,ElemConnectivityTriPE,myDistributedLoad,NodeForceBC,NodeCoor,A,EnNU,NodeDispBC);
    
    % IF Quad and tri are empty , run Bar 2d function
elseif isempty(ElemConnectivityQuadPS) && isempty(ElemConnectivityQuadPE) && isempty(ElemConnectivityTriPS) && isempty(ElemConnectivityTriPE)
    [RE,DF,myStress,myStrain]=Bar2d(ElemConnectivityBar,myDistributedLoad,NodeForceBC,NodeCoor,A,EnNU,NodeDispBC);
    
    % Else, run quad function
else
    
    [DF,RE]=ContQuad2d(ElemConnectivityQuadPS,ElemConnectivityQuadPE,myDistributedLoad,NodeForceBC,NodeCoor,A,EnNU,NodeDispBC);
end

% NOTE: Change all inv(K) to K\d.. More efficient

    
    
    
    
    
    
    
