function [DF,RE]=ContTri2d(ElemConnectivityTriPS,ElemConnectivityTriPE,myDistributedLoad,NodeForceBC,NodeCoor,A,EnNU,NodeDispBC);



    
    
    %% Plotting (PS and PE)
    if isempty(ElemConnectivityTriPE)           % it is a plain stress problem
        T = ElemConnectivityTriPS;
    else                                        % it is a plain strain problem
        T = ElemConnectivityTriPE;
    end
    
    
    
    
    %% Calculate element stifness matrix
    % Calculate D for plane stress and Plane strain
    E = EnNU(1,1);
    NU = EnNU(2,1);
    if isempty(ElemConnectivityTriPE)  % it is plane stress
        
        D = ( E  / (1-(NU  ^ 2))) * [1,NU,0; ...
            NU,1,0; ...
            0,0, (1-NU)/2];
        
    else  % it is a Plane strain
        
        D = (E / ((1+NU) * (1-2*NU))) * [1-NU,NU,0;...
            NU,1-NU,0;...
            0,0,(1-2*NU)/2];
        
    end

    
% get the size of T
            [m,n] = size(T);

     
    %Declare 3dmatrix for all elements X and Y coordinates
    myElementsCoor = zeros (1,6,m);
    %Declare 3dMatrix for element index
    index=zeros(1,6,m);
    %Declare 3dMatrix for B 
    myElementB=zeros(3,6,m);
    %Global Stifness matrix
    [a,b] = size(NodeCoor);
    K = zeros(a*2,a*2);
    %loop through elements
    
    
    for i = 1:m
        %Get coordinates in X and Y

            X1 = NodeCoor(T(i,2),2);
            Y1 = NodeCoor(T(i,2),3);
            X2 = NodeCoor(T(i,3),2);
            Y2 = NodeCoor(T(i,3),3);
            X3 = NodeCoor(T(i,4),2);
            Y3 = NodeCoor(T(i,4),3);
            myElementsCoor(:,:,i)=[X1,Y1,X2,Y2,X3,Y3];

        %Get Area
        Area = abs ( (X1*(Y2 - Y3) + X2*(Y3 - Y1) + X3*(Y1 - Y2)) / 2 );
        %Calculate N (shape function)
        syms X;
        syms Y;
        N1 = 1/(2*Area) * (X2*Y3 - X3*Y2 + (Y2-Y3) * X + (X3 - X2) * Y) ;
        N2 = 1/(2*Area) * (X3*Y1 - X1*Y3 + (Y3-Y1) * X + (X1 - X3) * Y) ;
        N3 = 1/(2*Area) * (X1*Y2 - X2*Y1 + (Y1-Y2) * X + (X2 - X1) * Y) ;
        %Calculate Be = nabla * Ne
        N1X = diff(N1,X);
        N2X = diff(N2,X);
        N3X = diff(N3,X);
        N1Y = diff(N1,Y);
        N2Y = diff(N2,Y);
        N3Y = diff(N3,Y);
        B = [N1X,0,N2X,0,N3X,0;...
            0,N1Y,0,N2Y,0,N3Y;...
            N1Y,N1X,N2Y,N2X,N3Y,N3X];
        % Add B to a 3d matrix
        myElementB(:,:,i) = B;
        %Calculate Ke = B(T) * D * B * thickness (A) * Area
        KE = B' * D * B * A * Area;
        %Comibine into K
 
            n1= T(i,2);
            n2= T(i,3);
            n3= T(i,4);

        index(:,:,i)=[2*n1-1 2*n1 2*n2-1 2*n2 2*n3-1 2*n3];
        K(index(:,:,i),index(:,:,i))=K(index(:,:,i),index(:,:,i))+KE;
    end
    
    %Get all x and y DOFS of all node
    myDOF=1:(a*2);
    %Make a copy of my myDOF
    myFreeDOFVector = myDOF;
    myEssentialDOFVector = myDOF;
    %Sort NodeDispBC rows according to first column
    NodeDispBC = sortrows(NodeDispBC,1);
    %Get all restricted nodes and remove from myFreeNodesVector
    [r,t]=size(NodeDispBC);
    for i=r:-1:1
        if NodeDispBC(i,2) ~= NodeDispBC(i,3)
            myFreeDOFVector(:,(NodeDispBC(i,1)*2))=[];
            myFreeDOFVector(:,(NodeDispBC(i,1)*2-1))=[];
        elseif NodeDispBC(i,2)== 1 && NodeDispBC(i,3)==1
            myFreeDOFVector(:,(NodeDispBC(i,1)*2-1))=[];
        elseif NodeDispBC(i,2)== 2 && NodeDispBC(i,3)==2
            myFreeDOFVector(:,(NodeDispBC(i,1)*2))=[];
        end
    end
    %Get KF
    KF = K(myFreeDOFVector,myFreeDOFVector);
    %Get my Essential nodes vector by deleting free nodes from it
    myEssentialDOFVector(myFreeDOFVector) = [];
    %Get KEF
    KEF = K(myEssentialDOFVector,myFreeDOFVector);
    %Get KE
    KE = K(myEssentialDOFVector,myEssentialDOFVector);
    
    
    %Construct Global force Vector
    myGlobalForceVector = zeros([1,a*2]);
    % if both mydistributedLoad and NodeForceBC are empty ie no forces just
    % displacements
    if  isempty(myDistributedLoad) && isempty(NodeForceBC) %%%%%%%%%%%%%%% case 1
        %do nothing
   
    elseif isempty(myDistributedLoad) && ~isempty(NodeForceBC) % use NodeForceBC %%%%%%%%%%%%%%% case 2
            
            %Populate this with the forces
            [L,b]=size(NodeForceBC) ;
            for i=1:L;
                if NodeForceBC(i,2)==1
                    myGlobalForceVector(1,NodeForceBC(i,1)*2-1)=NodeForceBC(i,3);
                else
                    myGlobalForceVector(1,NodeForceBC(i,1)*2)=NodeForceBC(i,3);
                end
            end
        else % use myDistributedLoad   %%%%%%%%%%%%%%% case 3
            %Solution for dload
            
            %Get size
            [a,b] = size(myDistributedLoad);
            for i = 1:a
                %get element number
                myElement = myDistributedLoad(i,1);
                %get myElement index
                myElementI = index(:,:,myElement);
                
                %get edge number 1 or 2 or 3
                if myDistributedLoad(i,2)==1 %acting on edge one
                    %get angle
                    x1 = myElementsCoor(1,1,myElement);
                    y1 = myElementsCoor(1,2,myElement);
                    x2 = myElementsCoor(1,3,myElement);
                    y2 = myElementsCoor(1,4,myElement);
                    %get normal to edge
                    dx=x2-x1;
                    dy=y2-y1;
                    %get angle f normal in angles
                    angle = atan((-dx-dx) /(dy+dy)) * 180/pi +180; %?
                    q1=cosd(angle) * myDistributedLoad(i,3);
                    q2=sind(angle) * myDistributedLoad(i,3);
                    Length = sqrt ((x2-x1)^2 + (y2-y1)^2);
                    Fe = [q1;q2;q1;q2;0;0] * 1/2 * A * Length;
                    myGlobalForceVector (:,myElementI) = myGlobalForceVector (:,myElementI) + Fe';
                    
                elseif myDistributedLoad(i,2)==2 %acting on edge two
                    %get angle
                    x1 = myElementsCoor(1,3,myElement);
                    y1 = myElementsCoor(1,4,myElement);
                    x2 = myElementsCoor(1,5,myElement);
                    y2 = myElementsCoor(1,6,myElement);
                    %get normal to edge
                    dx=x2-x1 ;
                    dy=y2-y1;
                    %get angle f normal in angles
                    angle = atan((-dx-dx) /(dy+dy)) * 180/pi + 180; %?
                    q1=cosd(angle) * myDistributedLoad(i,3);
                    q2=sind(angle) * myDistributedLoad(i,3);
                    Length = sqrt ((x2-x1)^2 + (y2-y1)^2);
                    Fe = [0;0;q1;q2;q1;q2] * 1/2 * A * Length;
                    myGlobalForceVector (:,myElementI) = myGlobalForceVector (:,myElementI) + Fe';
                    
                else %acting on edge three
                    %get angle
                    x1 = myElementsCoor(1,5,myElement);
                    y1 = myElementsCoor(1,6,myElement);
                    x2 = myElementsCoor(1,1,myElement);
                    y2 = myElementsCoor(1,2,myElement);
                    %get normal to edge
                    dx=x2-x1 ;
                    dy=y2-y1;
                    %get angle f normal in angles
                    angle = atan((-dx-dx) /(dy+dy)) * 180/pi + 180; %?
                    q1=cosd(angle) * myDistributedLoad(i,3);
                    q2=sind(angle) * myDistributedLoad(i,3);
                    %Another calculation for angles
                    % myCosValue=(dy+dy)/(sqrt(power((-dx-dx),2)+power((dy+dy),2)));
                    %mySinValue=(-dx-dx)/(sqrt(power((-dx-dx),2)+power((dy+dy),2)));
                    % q1=myCosValue * myDistributedLoad(i,3);
                    % q2=mySinValue * myDistributedLoad(i,3);
                    Length = sqrt ((x2-x1)^2 + (y2-y1)^2);
                    Fe = [q1;q2;0;0;q1;q2] * 1/2 * A * Length;
                    myGlobalForceVector (:,myElementI) = myGlobalForceVector (:,myElementI) + Fe';
                    
                end
            end
    end
        
        
        %Calculate DE
        %fix myEssentialDOFVector
        myEssentialDOFVector = myEssentialDOFVector';
        DE = myEssentialDOFVector;
        if  ~isempty(NodeDispBC)
            [a,b]=size(NodeDispBC);
            for i=1:a
                odd = NodeDispBC(i,1)*2 -1;
                even = NodeDispBC(i,1)*2;
                if NodeDispBC(i,2) ~= NodeDispBC(i,3)
                    DE(DE==odd)= NodeDispBC(i,4);
                    DE(DE==even)= NodeDispBC(i,4);
                elseif NodeDispBC(i,2)== 1 && NodeDispBC(i,3)==1
                    DE(DE==odd)= NodeDispBC(i,4);
                else % both equla two
                    DE(DE==even)= NodeDispBC(i,4);
                end
            end
        end
        
        
        %%%
        %Get freeDOFS force vector
        FF = myGlobalForceVector(myFreeDOFVector);
        %Transpose FF to get a column vector
        FF = FF';
        % Now to Calculate DF
        DF = inv(KF) * (FF - KEF'*DE)
        % Now to Calculate RE
        RE = KE * DE + KEF * DF
        
        %%%%%%%%%%%%%%%Plotting of deformed structure%%%%%%%%%%%%%%%%%
        %transpose DF so it matches myFreeNodesVector
        DF=DF';
        %make a new displacedNodeCoor from NodeCoor
        myDispNodeCoor = NodeCoor;
        %Make a global matrix for displacments - will help to print later
        myNodeDisp = zeros (a,2);
       %Go grab X and Y displacments of free nodes and add them to displacedNodeCoor
       [aa,b]=size(myFreeDOFVector);
        for i=b:-1:1
           if mod(myFreeDOFVector(:,i),2)==0 % ie even (y axis)
            [k,l] = find(myDispNodeCoor(:,1)==myFreeDOFVector(:,i)/2);
               myDispNodeCoor(k,3) = myDispNodeCoor(k,3) + DF(1,i);
                %Extra line for printing
               myNodeDisp(myFreeDOFVector(1,i)/2, 2) = DF(1,i);
            else
                [k,l] = find(myDispNodeCoor(:,1)==(myFreeDOFVector(:,i)+1)/2);
                 myDispNodeCoor(k,2) = myDispNodeCoor(k,2) + DF(1,i);
                %Extra line for printing
                myNodeDisp((myFreeDOFVector(1,i)+1)/2, 1) = DF(1,i);
           end
        end
        N = myDispNodeCoor;
        N(:,1)=[];
        [a,b] = size(N);
        myZeros = zeros(a,1);
        N = [N myZeros];
      
        
        
        %%%%%%%%%%%%%%% Calculation of stress and strain %%%%%%%%%%%%%%%%%
        %Calculate de for each element
        allD = zeros(1,a*2);
        allD(myFreeDOFVector) = DF; % so elegant...
     
        
        % calculate strain = B of element * de
        % declare 3d matrix for element strain and stress
        myElementStrain=zeros(3,1,m);
        myElementStress=zeros(3,1,m);
        %Loop through elements
        for i = 1:m
            de = allD(index(:,:,i));
            myElementStrain (:,:,i) = myElementB(:,:,i)  * de';
            % now calculate stress = De * strain
        myElementStress (:,:,i) = D * myElementStrain (:,:,i);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
      S11 = zeros(m,1);
      S22 = zeros(m,1);
      S12 = zeros(m,1);
      E11 = zeros(m,1);
      E22 = zeros(m,1);
      G12 = zeros(m,1);
      for i = 1:m
         S11(i,1) =  myElementStress (1,1,i);
         S22(i,1) =  myElementStress (2,1,i);
         S12(i,1) =  myElementStress (3,1,i);
         E11(i,1) =  myElementStrain (1,1,i);
         E22(i,1) =  myElementStrain (2,1,i);
         G12(i,1) =  myElementStrain (3,1,i);
      end
      
 T(:,1)=[];
        trisurf(T,N(:,1),N(:,2),N(:,3),S11); %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Change plot value here
        view(2);
        
      
        %Convert allD into a matrix
        allD = vec2mat(allD,2);
        %Convert myGlobal force into a matrix
        myGlobalForceVector = vec2mat(myGlobalForceVector,2);
        %get the Essential nodes and put their reaction force in the correct
        %index in the myGlobalForce index (Compare myEssentialNodesVector to RE)
        [j,k] = size(RE);
        for i = 1:j
            if mod(myEssentialDOFVector(i,1),2)==0 %ie even , Y Coor
                myGlobalForceVector (myEssentialDOFVector(i,1)/2,2) = RE (i,1);
            else
                myGlobalForceVector ((myEssentialDOFVector(i,1)+1)/2,1) = RE (i,1) ;
            end
        end
        
        
          %Open file to print to
        fid = fopen('HW6.txt','w');
        
        
            %Print Node
        fprintf(fid,'Node\tU_1\t\t\t\t\tU_2\t\t\t\t\tF_1\t\t\t\t\tF_2\n');
        
     
        for i = 1:4
            fprintf(fid,'%d\t\t',i); %allD myGlobalForceVector
            fprintf(fid,'%f\t\t\t',allD(i,:));
            fprintf(fid,'%f\t\t\t',myGlobalForceVector(i,:));
            fprintf(fid,'\n');
        end
        

        fprintf(fid,'\n');
                
        %Print Stress & Strain
        fprintf(fid,'Element\tS11\t\t\t\t\tS22\t\t\t\t\tS12\t\t\t\t\tE11\t\t\t\t\tE22\t\t\t\t\tE12\n');
        
        for i = 1:size(S11,1)
            fprintf(fid,'%d\t\t',i);
            fprintf(fid,'%f\t\t\t',S11(i,:));
            fprintf(fid,'%f\t\t\t',S22(i,:));
            fprintf(fid,'%f\t\t\t',S12(i,:));
             fprintf(fid,'%f\t\t\t',E11(i,:));
            fprintf(fid,'%f\t\t\t',E22(i,:));
            fprintf(fid,'%f\t\t\t',G12(i,:));
            
            
            fprintf(fid,'\n');
        end
        
        fclose(fid)
        
        %%%%%%%%%%%%%%for quad inner if
    end


