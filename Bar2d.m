function [RE,DF,myStress,myStrain]=Bar2d(ElemConnectivityBar,myDistributedLoad,NodeForceBC,NodeCoor,A,EnNU,NodeDispBC);




        %Get size of Element Connectivity
        [m,n] = size(ElemConnectivityBar);
        %myCombinedMatrix shows the node coord values to be connected with elements. Only works for 2D
        for i1=1:m
            myCombinedMatrix(i1,1)= NodeCoor(ElemConnectivityBar(i1,2),2);
            myCombinedMatrix(i1,2)= NodeCoor(ElemConnectivityBar(i1,2),3);
            myCombinedMatrix(i1,3)= NodeCoor(ElemConnectivityBar(i1,3),2);
            myCombinedMatrix(i1,4)= NodeCoor(ElemConnectivityBar(i1,3),3);
        end
        
        %Figuring out angles of rotated elements -  Only works for 2D
        for i3 = 1:m
            x1=myCombinedMatrix(i3,1);
            y1=myCombinedMatrix(i3,2);
            x2=myCombinedMatrix(i3,3);
            y2=myCombinedMatrix(i3,4);
            myCosValue(i3)=(x2-x1)/(sqrt(power((y2-y1),2)+power((x2-x1),2)));
            mySinValue(i3)=(y2-y1)/(sqrt(power((y2-y1),2)+power((x2-x1),2)));
            k(i3)= A * EnNU(1,1) / (sqrt(power((y2-y1),2)+power((x2-x1),2)));
            myOriginalLength(i3,1) = sqrt(power((y2-y1),2)+power((x2-x1),2));
            
            % myAngles=atan2(y2-y1, x2-x1)*180/pi older technique - not accurate
        end
        
        %Get angles in Degrees - only frst row
        %myAngles=acosd(myCosValue) % dont do this - will make problems
        
        %Get the number of nodes
        [a,b]=size(NodeCoor);
        %Construct empty global stifness matrix
        myGlobalStifnessMatrix = zeros(a*2);
        for i=1:m
            % Rotation
            myRotationMatrix=[myCosValue(1,i),mySinValue(1,i),0,0;0,0,myCosValue(1,i),mySinValue(1,i)]
            % to get the seperate stifness matrix for each element
            myElementStifnessMatrix= k(1,i)*myRotationMatrix' *[1,-1;-1,1]*  myRotationMatrix;
            %Split local stifness matrix into four quarters - go to my nodes and get
            %node number
            n1= ElemConnectivityBar(i,2);
            n2= ElemConnectivityBar(i,3);
            index=[2*n1-1 2*n1 2*n2-1 2*n2];
            myGlobalStifnessMatrix(index,index)=myGlobalStifnessMatrix(index,index)+myElementStifnessMatrix;
            % the top code does that simply...
            % myGlobalStifnessMatrix((2*node1-1):(2*node1), (2*node1-1):(2*node1)) = myGlobalStifnessMatrix((2*node1-1):(2*node1), (2*node1-1):(2*node1)) + myElementStifnessMatrix(1:2,1:2);
            % myGlobalStifnessMatrix((2*node2-1):(2*node2), (2*node2-1):(2*node2)) = myGlobalStifnessMatrix((2*node2-1):(2*node2), (2*node2-1):(2*node2)) + myElementStifnessMatrix(3:4,3:4);
            % myGlobalStifnessMatrix((2*node1-1):(2*node1), (2*node2-1):(2*node2)) = myGlobalStifnessMatrix((2*node1-1):(2*node1), (2*node2-1):(2*node2)) + myElementStifnessMatrix(3:4,1:2);
            % myGlobalStifnessMatrix((2*node2-1):(2*node2), (2*node1-1):(2*node1)) = myGlobalStifnessMatrix((2*node1-1):(2*node1), (2*node2-1):(2*node2)) + myElementStifnessMatrix(1:2,3:4);
        end
        
        %Get all x and y positions of all node
        myNodesVector=1:(a*2);
        %Make a copy of my myFreeNodesVector
        myFreeNodesVector = myNodesVector;
        %Sort NodeDispBC rows according to first column
        NodeDispBC = sortrows(NodeDispBC,1);
        %Get all restricted nodes and remove from myFreeNodesVector
        [r,t]=size(NodeDispBC);
        for i5=r:-1:1
            if NodeDispBC(i5,2) ~= NodeDispBC(i5,3)
                myFreeNodesVector(:,(NodeDispBC(i5,1)*2))=[];
                myFreeNodesVector(:,(NodeDispBC(i5,1)*2-1))=[];
            elseif NodeDispBC(i5,2)== 1 && NodeDispBC(i5,3)==1
                myFreeNodesVector(:,(NodeDispBC(i5,1)*2-1))=[];
            elseif NodeDispBC(i5,2)== 2 && NodeDispBC(i5,3)==2
                myFreeNodesVector(:,(NodeDispBC(i5,1)*2))=[];
                
            end
        end
        
        %Copy myGlobalStifnessMatrix into..
        KF=myGlobalStifnessMatrix;
        KE=myGlobalStifnessMatrix;
        KFE=myGlobalStifnessMatrix;
        KEF=myGlobalStifnessMatrix;
        
        %Produce KF
        for i = 2*a:-1:1
            if ~max(i==myFreeNodesVector);
                KF(i,:)=[];
                KF(:,i)=[];
            end
        end
        %Produce KE
        for i = 2*a:-1:1
            if max(i==myFreeNodesVector);
                KE(i,:)=[];
                KE(:,i)=[];
            end
        end
        %Produce KEF
        for i = 2*a:-1:1
            if ~max(i==myFreeNodesVector);
                KEF(:,i)=[];
            end
            if max(i==myFreeNodesVector);
                KEF(i,:)=[];
            end
        end
        %Produce KFE
        for i = 2*a:-1:1
            if max(i==myFreeNodesVector);
                KFE(i,:)=[];
            end
            if max(i==myFreeNodesVector);
                KFE(:,i)=[];
            end
        end
        
        %Construct myGlobalForce, all zeros to start with...
        myGlobalForce = zeros(1,a*2) ;
        %Populate this with the forces
        [L,b]=size(NodeForceBC) ;
        for i=1:L;
            if NodeForceBC(i,2)==1
                myGlobalForce(1,NodeForceBC(i,1)*2-1)=NodeForceBC(i,3);
            else
                myGlobalForce(1,NodeForceBC(i,1)*2)=NodeForceBC(i,3);
            end
        end
        %Copy myGlobalForce to FF
        FF=myGlobalForce;
        %  Now to delete the restrained nodes from FF (They have 0 values)
        for i = a*2:-1:1
            if ~max(i==myFreeNodesVector);
                FF(:,i)=[];
            end
        end
        %Transpose FF to get a column vector
        FF = FF';
        % Now to Calculate DF
        DF = inv(KF) * FF
        % Now to Calculate RE
        RE = KEF* DF
        
        
        
        %transpose DF so it matches myFreeNodesVector
        DF=DF'
        %make a new displacedNodeCoor from NodeCoor
        myDispNodeCoor = NodeCoor;
        %Make a global matrix for displacments - will help to print later
        myNodeDisp = zeros (a,2);
        %Go grab X and Y displacments of free nodes and add them to displacedNodeCoor
        [aa,b]=size(myFreeNodesVector);
        for i=b:-1:1
            if mod(myFreeNodesVector(:,i),2)==0 % ie even (y axis)
                [k,l] = find(myDispNodeCoor(:,1)==myFreeNodesVector(:,i)/2);
                myDispNodeCoor(k,3) = myDispNodeCoor(k,3) + DF(:,i);
                %Extra line for printing
                myNodeDisp(myFreeNodesVector(1,i)/2, 2) = DF(1,i);
            else
                [k,l] = find(myDispNodeCoor(:,1)==(myFreeNodesVector(:,i)+1)/2);
                myDispNodeCoor(k,2) = myDispNodeCoor(k,2) + DF(:,i);
                %Extra line for printing
                myNodeDisp((myFreeNodesVector(1,i)+1)/2, 1) = DF(1,i);
            end
        end
        
        
        %Plot structure before load application
        xCoor=NodeCoor(:,2);
        yCoor=NodeCoor(:,3);
        %Plot myDispNodeCoor after load application
        xDispCoor=myDispNodeCoor(:,2);
        yDispCoor=myDispNodeCoor(:,3);
        %plot both Nodes together
        plot(xCoor,yCoor,'o','color','b')
        plot(xDispCoor,yDispCoor,'o','color','r')
        
        
        %Now plot members before load application
        %  for i2 = 1:m
        %     line([myCombinedMatrix(i2, 1), myCombinedMatrix(i2, 3)], [myCombinedMatrix(i2, 2), myCombinedMatrix(i2, 4)])
        %  end
        %Get myCombinedMatrix1 for Displaced nodes
        for i=1:m
            myCombinedMatrix1(i,1)= myDispNodeCoor(ElemConnectivityBar(i,2),2);
            myCombinedMatrix1(i,2)= myDispNodeCoor(ElemConnectivityBar(i,2),3);
            myCombinedMatrix1(i,3)= myDispNodeCoor(ElemConnectivityBar(i,3),2);
            myCombinedMatrix1(i,4)= myDispNodeCoor(ElemConnectivityBar(i,3),3);
        end
        %Now plot members after load application
        for i2 = 1:m
            line([myCombinedMatrix1(i2, 1), myCombinedMatrix1(i2, 3)], [myCombinedMatrix1(i2, 2), myCombinedMatrix1(i2, 4)],'color','r')
        end
        
        
        % We have myCombinedMatrix for element coor before load application and
        % myCombinedMatrix1 for element coor after load application
        
        % get the elongated length for each element and calculate the stress and strain in one
        % loop
        for i=1:m
            x1 = myCombinedMatrix1(i,1);
            y1 = myCombinedMatrix1(i,2);
            x2 = myCombinedMatrix1(i,3);
            y2 = myCombinedMatrix1(i,4);
            myElongatedLength(i,1) = x2 * myCosValue(1,i) + y2 * mySinValue(1,i) - ( x1 * myCosValue(1,i) + y1 * mySinValue(1,i));
            myLengthDifference (i,1) = myElongatedLength(i,1) -  myOriginalLength(i,1);
            myStrain (i,1) = myLengthDifference (i,1) / myOriginalLength(i,1)
            myStress (i,1) = EnNU(1,1) * myStrain (i,1)
        end
        
        %Some postprocessing to combine values from RE with their respective nodes
        myEssentialNodesVector = 1:a*2;
        for j=a*2 : -1: 1
            if max(j==myFreeNodesVector);
                myEssentialNodesVector(:,j)=[];
            end
        end
        %transpose RE so it matches myEssentialNodesVector
        RE = RE';
        %Convert myGlobal force into a matrix
        myGlobalForce = vec2mat(myGlobalForce,2);
        %get the Essential nodes and put their reaction force in the correct
        %index in the myGlobalForce index (Compare myEssentialNodesVector to RE)
        [j,k] = size(RE);
        for i = 1:k
            if mod(myEssentialNodesVector(1,i),2)==0 %ie even , Y Coor
                myGlobalForce (myEssentialNodesVector(1,i)/2,2) = RE (1,i);
            else
                myGlobalForce ((myEssentialNodesVector(1,i)+1)/2,1) = RE (1,i) ;
            end
        end
        
        
        %Open file to print to
        fid = fopen('MyOutputFile.txt','w');
        
        %Print Stress & Strain
        fprintf(fid,'Element\tStress\t\t\tStrain\n');
        
        for i = 1:size(myStress,1)
            fprintf(fid,'%d\t\t',i);
            fprintf(fid,'%f\t\t\t',myStress(i,:));
            fprintf(fid,'%f',myStrain(i,:));
            fprintf(fid,'\n');
        end
        
        fprintf(fid,'\n');
        
        %Print Stress & Strain
        fprintf(fid,'Node\tU_1\t\t\t\t\tU_2\t\t\t\t\tF_1\t\t\t\t\tF_2\n');
        
        %Print Displacements
        for i = 1:size(myNodeDisp,1)
            fprintf(fid,'%d\t\t',i);
            fprintf(fid,'%f\t\t\t',myNodeDisp(i,:));
            fprintf(fid,'%f\t\t\t',myGlobalForce(i,:));
            fprintf(fid,'\n');
        end
        
        fclose(fid);


end
