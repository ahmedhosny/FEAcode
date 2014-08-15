function [ElemConnectivityBar,ElemConnectivityTriPS,ElemConnectivityTriPE,ElemConnectivityQuadPS,ElemConnectivityQuadPE,myDistributedLoad,NodeForceBC,NodeCoor,A,EnNU,NodeDispBC]=Readinput(NameTextFile);

%Open up the file
myFileID=fopen(NameTextFile);
%initiate a variable - to allow the big if statement to run
ElemConnectivityBar=[];
ElemConnectivityTriPS=[];
ElemConnectivityTriPE=[];
ElemConnectivityQuadPS =[];
ElemConnectivityQuadPE =[];
myDistributedLoad =[];
NodeForceBC = [];
%Go through file
myStop=0;
while myStop==0
    myLine=fgetl(myFileID);
    myLine=lower(myLine); %to have it read lowercase
    if length(myLine)>4 && strcmp(myLine,'*node print')==0
        switch myLine(1:5); % read from 1 to 5 number of characters in a string array
            case '*node'
                NodeCoor=fscanf(myFileID, '%d, %f, %f' , [3,Inf])' % to transpose
            case '*elem'
                % If statement in switchcase to read different element types ,
                % now I have different variable names...Solved with big if
                % statement
                if strcmp(myLine(16:19),'t2d2')==1
                    ElemConnectivityBar=fscanf(myFileID, ' %d, %d, %d' , [3,Inf])'
                elseif  strcmp(myLine(16:19),'cps3')==1
                    ElemConnectivityTriPS=fscanf(myFileID, ' %d, %d, %d, %d' , [4,Inf])'
                elseif strcmp(myLine(16:19),'cpe3')==1
                    ElemConnectivityTriPE=fscanf(myFileID, ' %d, %d, %d, %d' , [4,Inf])'
                   elseif  strcmp(myLine(16:19),'cps4')==1
                    ElemConnectivityQuadPS=fscanf(myFileID, ' %d, %d, %d, %d, %d' , [5,Inf])'  
                  elseif  strcmp(myLine(16:19),'cpe4')==1
                    ElemConnectivityQuadPE=fscanf(myFileID, ' %d, %d, %d, %d, %d' , [5,Inf])'
                end
            case '*soli'
                A=fscanf(myFileID, '%f')
            case '*elas' % Continuum elements have E,nu
                EnNU=fscanf(myFileID, '%f,%f') %first row is E, second is nu - if it exists
            case '*boun'
                NodeDispBC=fscanf(myFileID, '%d, %d, %d, %f' , [4,Inf])'
            case '*dloa' % for distributed load, P1,2,3 are the edges of the triangle on which the load is acting
                myDistributedLoad=fscanf(myFileID, '%d, %*c%d,  %f' , [3,Inf])'
            case '*cloa'
                NodeForceBC=fscanf(myFileID, '%d, %d, %f' , [3,Inf])'
                %dummy last case
            case '*end '
                %Do Nothing
                myStop=1;
        end
        
    end

    
end