function [Network_Distance, Bijection, B] = NetworkCompare(Matrix1_Initial,Matrix2_Initial)%,View_Graph,Ids1,Ids2)
    
%     close all hidden
%     clear all
%     clc
% 
%    Matrix2_Initial = [     0     0     1     1     0
%      1     0     1     0     0
%      0     1     0     1     0
%      0     0     1     1     0
%      0     1     0     0     0];
%     
%     Matrix1_Initial = [     0     0     1     1     0
%      0     0     0     0     0
%      0     0     0     1     0
%      0     0     1     1     0
%      0     0     0     0     0];



    % Normalize both matrices so any scale differences don't play a factor
    Matrix1_Initial = Matrix1_Initial/max(max(Matrix1_Initial));                            
    Matrix2_Initial = Matrix2_Initial/max(max(Matrix2_Initial));



    % Set the smaller matrix to be the first one, so will map the smaller onto
    % the larger
    Switch = 0;
    Matrix1 = [];
    Matrix2 = [];
    if length(Matrix1_Initial) <= length(Matrix2_Initial)
        Matrix1 = Matrix1_Initial;
        Matrix2 = Matrix2_Initial;
        
%         Ids1 = Ids1_Initial;
%         Ids2 = Ids2_Initial;        
    else
        Matrix1 = Matrix2_Initial;
        Matrix2 = Matrix1_Initial;
        Switch = 1;
    
%         Ids1 = Ids2_Initial;
%         Ids2 = Ids1_Initial;
    end

    
%     % Pad the smaller matrix to be the same size as the larger one
%     if length(Matrix1) <= length(Matrix2)
%         
%         Matrix1_Extended = zeros(length(Matrix2));
%         for i = 1:length(Matrix2)
%             for j = 1:length(Matrix2)
%                 if i <= length(Matrix1) & j <= length(Matrix1)
%                     Matrix1_Extended(i,j) = Matrix1(i,j);
%                 end
%             end
%         end
%         Matrix1 = Matrix1_Extended;
%         
%     else
%         
%         Matrix2_Extended = zeros(length(Matrix1));
%         for i = 1:length(Matrix1)
%             for j = 1:length(Matrix1)
%                 if i <= length(Matrix2) & j <= length(Matrix2)
%                     Matrix2_Extended(i,j) = Matrix2(i,j);
%                 end
%             end
%         end
%         Matrix2 = Matrix2_Extended;
%         
%     end
            
%     This is for the ecological matrices ONLY!!!
%     [r,c] = size(Matrix1);
%     Matrix1 = Matrix1(2:r, 2:c);
% 
%     [r,c] = size(Matrix2);
%     Matrix2 = Matrix2(2:r, 2:c);

    
    % Create IDs
    Ids = [];
    for i = 1:max(length(Matrix1),length(Matrix2));
        Ids = [Ids; i];
    end
    Ids = cellstr(num2str(Ids));

    
    % Run the flow analysis for the number of edges, to ensure every node and
    % edge that can be reached froma given starting node will be (but if
    % the network has no connections then run it for one time-step because 
    % will get an error if try to do none)
    Duration = max([1,sum(sum(Matrix1)),sum(sum(Matrix2))]);
    
    Duration1 = graphallshortestpaths(sparse(Matrix1>0)); % Use Matrix1 > 0 because graphallshortestpaths assumes weights are durations of flow, not quantities 
    Duration1(isinf(Duration1))=0;
    Duration1 = 2 * max(max(Duration1));
    if Duration1 == 1 % Make sure that there are at least two time steps or no curve can be fit to Gini score vs Duration later on
        Duration1 = 2;
    end
    
    
    Duration2 = graphallshortestpaths(sparse(Matrix2>0)); % Use Matrix2 > 0 because graphallshortestpaths assumes weights are durations of flow, not quantities
    Duration2(isinf(Duration2))=0;
    Duration2 = 2 * max(max(Duration2));
    if Duration2 == 1 % Make sure that there are at least two time steps or no curve can be fit to Gini score vs Duration later on
        Duration2 = 2;
    end

   
D1 = Duration1;

Gini1 = zeros(D1,length(Matrix1)); 
for d = 1:D1    
    Duration1 = d;
    
    % Flow analysis for each node of the first matrix
    Matrix1_Flows = [];
    Matrix1_Flows_unsorted = [];
    parfor z = 1:length(Matrix1)

        % Will update all nodes only after doing the analysis on all of them,
        % so use a two-matrix system: Energy_Old & Energy_New
        Energy_Old1 = [];
        Energy_Added = [];
        for x = 1:length(Matrix1)
            Energy_Old1 = [Energy_Old1; x 0];
            Energy_Added = [Energy_Added; x 0];
        end

        % Initial conditions (IC) of 100 units of energy at starting node z
        IC = 100;
        Origin = z;
        Energy_Old1(Origin,2) = IC;

        [r c] = size(Energy_Old1);

        ids = Energy_Old1(:,1);

        Energy_New = Energy_Old1;
        Energy = [];
        Time = [];
        
% view(biograph(sparse(Matrix1)));

        % Let energy flow through the network starting at node z for 'Duration'
        % time steps
        for t = 1:Duration1
            for i = 1:r
% i
                Edges = [];

                % Only consider nodes with nonzero energy
                if Energy_Old1(i,2) ~= 0

                    % Find all edges for a node with energy in it
                    for j = 1:length(Matrix1)
                        if Matrix1(Energy_Old1(i),j) ~= 0
                            Edges = [Edges j];
                        end 
                    end 

                    % Let energy equally flow along the outgoing edges
                    if length(Edges) == 0
                        Flow = 0;
                    else
                        % Cannot have more than 1 unit of energy flow along an
                        % edge per time step
                        Flow = Matrix1(i,Edges)*(Energy_Old1(i,2)/sum(Matrix1(i,Edges))); %min(1,Energy_Old1(i,2)/length(Edges));
% Flow
                    end

                    % Add the flowed energy to their new nodes
                    for k = 1:length(Edges)
                        Energy_New(Edges(k),2) = Energy_New(Edges(k),2) + Flow(k);
                    end
% Energy_New;
                    % Keep track of just what passed through each node at
                    % each time step (what was added only, not what was
                    % lost)
                    for k = 1:length(Edges)
                        Energy_Added(Edges(k),2) = Energy_Added(Edges(k),2) + Flow(k);
                    end                    
% Energy_Added;

                    % Subtract off the lost energy from the node under
                    % consideration except the original node which only gains
                    % energy
%                     if i == Origin
%                     else
                        Energy_New(i,2) = Energy_New(i,2) - sum(Flow); %abs(Energy_New(i,2) - Flow*length(Edges));
%                     end
% Energy_New;

                end
% Edges                
% Flow
% pause
            
            end

            % Update the old energy levels to be the new ones
            Energy_Old1 = Energy_New; 

            % Keep track of     
            Time = [Time t];
            Energy = [Energy Energy_Added(:,2)];

% t
% Energy
% pause

        end
% error
        
        % Plot energy per node over time
        % plot(Time(1:Duration),Energy(:,[1:Duration]))

        % Assign final energies to their ID number
        if Duration1 == 1 % will get an error because sum(Energy')' becomes a single number instead of a column vector
            Energy = [ids Energy];
        else
            Energy = [ids sum(Energy')'];
        end
        Matrix1_Flows_unsorted = [Matrix1_Flows_unsorted Energy(:,2)];
% Origin
% Duration1
% Energy 
% sum(Energy)
% error

        % Sort the final node energies
        Sorted = sortrows(Energy,-2);
       

        % Assign the actual ID name to the final node energies
        ids_sorted = [];
        for i = 1:length(Sorted)
            ids_sorted = [ids_sorted; Ids(Sorted(i,1)) Energy(Sorted(i,1),2)];
        end

        
        % Add to growing matrix of final flow distributions for each starting
        % node
        Matrix1_Flows = [Matrix1_Flows cell2mat(ids_sorted(:,2))];

    end

    Matrix1_Flows;
    Matrix1_Flows_unsorted;
    
    for g = 1:length(Matrix1_Flows)
        Gini1(d,g) = Gini(Matrix1_Flows(:,g));
    end
    
    X1 = 0:1/(size(Gini1,1)-1):1;
    X1 = X1';
    
%     X1
%     Gini1
%     plot(X1,Gini1)
% %     pause
%     clc
    
end    


% Gini1
% X1









D2 = Duration2;    
Gini2 = zeros(D2,length(Matrix2)); 
for d = 1:D2    
    Duration2 = d;
    
    % Repeat the process, but for the second matrix (see comments above - the
    % process is identical)
    Matrix2_Flows = [];
    Matrix2_Flows_unsorted = [];
    parfor z = 1:length(Matrix2)

        % Flow analysis

        Energy_Old2 = [];
        Energy_Added = [];
        for x = 1:length(Matrix2)
            Energy_Old2 = [Energy_Old2; x 0];
            Energy_Added = [Energy_Added; x 0];            
        end


        IC = 100;
        Origin = z;
        Energy_Old2(Origin,2) = IC;

        [r c] = size(Energy_Old2);

        ids = Energy_Old2(:,1);

        Energy_New = Energy_Old2;
        Energy = [];
        Time = [];


        for t = 1:Duration2
            for i = 1:r

                Edges = [];

                if Energy_Old2(i,2) ~= 0
                    for j = 1:length(Matrix2)
                        if Matrix2(Energy_Old2(i),j) ~= 0
                            Edges = [Edges j];
                        end 
                    end


                    if length(Edges) == 0
                        Flow = 0;
                    else
                        Flow = Matrix2(i,Edges)*(Energy_Old2(i,2)/sum(Matrix2(i,Edges))); %min(1,Energy_Old2(i,2)/length(Edges));
                    end

                    for k = 1:length(Edges)
                        Energy_New(Edges(k),2) = Energy_New(Edges(k),2) + Flow(k); 
                    end
                    
                    % Keep track of just what passed through each node at
                    % each time step (what was added only, not what was
                    % lost)
                    for k = 1:length(Edges)
                        Energy_Added(Edges(k),2) = Energy_Added(Edges(k),2) + Flow(k);
                    end 

%                     if i == Origin
%                     else
                        Energy_New(i,2) = Energy_New(i,2) - sum(Flow); %abs(Energy_New(i,2) - Flow*length(Edges));
%                     end


                end
            end

            Energy_New; 
            Energy_Old2 = Energy_New;


            Time = [Time t];
            Energy = [Energy Energy_Added(:,2)];

        end


        
        Energy_New;
        Energy;
        sum(Energy_New(:,2));
        Time;

        % plot(Time(1:Duration),Energy(:,[1:Duration]))

        if Duration2 == 1
            Energy = [ids Energy];
        else
            Energy = [ids sum(Energy')'];
        end
        Matrix2_Flows_unsorted = [Matrix2_Flows_unsorted Energy(:,2)];
        % Energy_ES = Energy(34:length(Energy),:);

        Sorted = sortrows(Energy,-2);

        ids_sorted = [];
        for i = 1:length(Sorted)
            ids_sorted = [ids_sorted; Ids(Sorted(i,1)) Energy(Sorted(i,1),2)];

        end

        ids_sorted;

        One_Under_Consideration = Ids(Origin);

        Matrix2_Flows = [Matrix2_Flows cell2mat(ids_sorted(:,2))];
    end

    Matrix2_Flows;
    Matrix2_Flows_unsorted;
    
        for g = 1:length(Matrix2_Flows)
            Gini2(d,g) = Gini(Matrix2_Flows(:,g));
        end
        
    X2 = 0:1/(size(Gini2,1)-1):1;
    X2 = X2';
    
end





% disp('end of flow')






%     % Deal with a weird rounding error where zeros (which I'm sure should be 
%     % zeros) are sometimes computed as a very small number (i.e. 1e-14)
%     Matrix1_Flows = round(10^10*Matrix1_Flows)/10^10;
%     Matrix2_Flows = round(10^10*Matrix2_Flows)/10^10;



    Size1 = length(Matrix1);
    Size2 = length(Matrix2);

    % Compute pairwise distances between nodes
    Similarity = zeros(Size1,Size2)+10;
    for i = 1:Size1
        for j = 1:Size2
            
            % Deal with the situation of having no energy in any of the
            % nodes (how would this happen?)
%             if sum(Matrix1_Flows(:,i)==0) == length(Matrix1_Flows(:,i)) 
%                 Gini_Matrix1_Flows = 1;
%             else
%                 Gini_Matrix1_Flows = Gini(Matrix1_Flows(:,i));
%             end
%             
%             if sum(Matrix2_Flows(:,j)==0) == length(Matrix2_Flows(:,j)) 
%                 Gini_Matrix2_Flows = 1;
%             else
%                 Gini_Matrix2_Flows = Gini(Matrix2_Flows(:,j));
%             end
%             
%             Similarity(i,j) = abs(Gini_Matrix1_Flows - Gini_Matrix2_Flows);
% %             Similarity(i,j) = pdist2((Matrix1_Flows(:,i)/sum(Matrix1_Flows(:,i)))',(Matrix2_Flows(:,j)/sum(Matrix2_Flows(:,j)))'); %/VectorFar(Matrix1_Flows(:,i));

            Similarity(i,j) = AreaWhole(X1,Gini1(:,i),X2,Gini2(:,j),1000);
%             pause
        end
    end
    
% disp('starting hungarian algorithm')    
    
    % To deal with the parts of the Similarity matrix which must be square
    % but which don't exist because of unequal size matrices, just put the
    % value of 10 there (0 <= Gini <= 1)
    Similarity(isnan(Similarity)) = 10;
    
    % Hungarian Method
    [assignment,cost,dMat] = munkres(Similarity);
    % Similarity = dMat;
    
    % Reshape the similarity and dmat matrices
    Similarity = Similarity(1:Size1,1:Size2);
    dMat = dMat(1:Size1,1:Size2);
    
    
    % Create the bijection by sorting the pairwise node distances in descending
    % order
% % %     Bijection = [Energy_Old1(:,1)];      I CHANGED THIS JUST FOR parfor TO WORK ON JANUS
    B = [];
    for i = 1:Size1
% % %         B = [B; zeros(Size2,1)+i Energy_Old2(:,1) Similarity(i,:)' dMat(i,:)'];   I CHANGED THIS JUST FOR parfor TO WORK ON JANUS
        B = [B; zeros(Size2,1)+i (1:Size2)' Similarity(i,:)' dMat(i,:)'];
    end



    % Find the optimal solution and give them a value of -1 so will be at
    % the top of the list when sorted (sometimes there are more 0 distances
    % than just the optimal ones at the end of the Hungarian Method)
    H = mod(find(munkres(Similarity) == 1),length(Similarity));
    C = zeros(length(H),4);
    for i = 1:length(H)
        if H(i) == 0
            C(i,:) = [length(H) i Similarity(length(H),i) dMat(length(H),i)];
        else
            C(i,:) = [H(i) i Similarity(H(i),i) dMat(H(i),i)];
        end
    end
    C;

    % The part that gives the value of -1
    [r_B,c_B] = size(B);
    [r_C,c_C] = size(C);
    for i = 1:r_B
        for j = 1:r_C
            if sum(B(i,1:2) == C(j,1:2)) == 2
                B(i,4) = -1;
            end
        end
    end


    B = sortrows(B,4);

    Bijection = [B(1,:)]; %0];

    rows = 1;
    i = 2;
    Connectivity1_all = graphallshortestpaths(sparse(Matrix1));
    
    while rows < min(Size1,Size2) & i < length(B)
        Bijection_Trial = [Bijection; B(i,:)]; %0];
        [r_BT,c_BT] = size(Bijection_Trial);

% %         Connectivity1 = Connectivity1_all(Bijection_Trial(:,1),Bijection_Trial(:,1));
% %         
% %         % Build up the second matrix with the mappings in Bijection_Trial
% %         % so connectivity can include the entire network
% %         BT_Connectivity = [[Bijection_Trial(:,1); setdiff(cell2mat(Ids),Bijection_Trial(:,1))], [Bijection_Trial(:,2); setdiff(cell2mat(Ids),Bijection_Trial(:,2))]];
% %         BT_Connectivity = sortrows(BT_Connectivity,2);
% %         Matrix2_Connectivity = zeros(Size2);
% %         for j = 1:Size2
% %             for k = 1:Size2
% %                 Matrix2_Connectivity(j,k) = Matrix2(BT_Connectivity(j,1),BT_Connectivity(k,1));
% %             end
% %         end
% %         Connectivity2_all = graphallshortestpaths(sparse(Matrix2_Connectivity));
% %         Connectivity2 = Connectivity2_all(Bijection_Trial(:,2),Bijection_Trial(:,2));
% %         
% % %         Connectivity1 = graphallshortestpaths(sparse(Matrix1(Bijection_Trial(:,1),Bijection_Trial(:,1))));
% %         Connectivity1(isinf(Connectivity1))=length(Bijection_Trial)+1;
% % %         Connectivity2 = graphallshortestpaths(sparse(Matrix2(Bijection_Trial(:,2),Bijection_Trial(:,2))));
% %         Connectivity2(isinf(Connectivity2))=length(Bijection_Trial)+1;
% % 
% %         Connectivity_Threshold = 1000000000000000000; %10.5*r_BT^2;

        if sum(B(i,1) == Bijection(:,1)) == 0 & sum(B(i,2) == Bijection(:,2)) == 0 %%& sum(sum(abs(Connectivity1-Connectivity2))) <= Connectivity_Threshold
    %         sum(sum(abs(graphallshortestpaths(sparse(Matrix1(Bijection_Trial(:,1),Bijection_Trial(:,1))))(isnan(graphallshortestpaths(sparse(Matrix1(Bijection_Trial(:,1),Bijection_Trial(:,1))))))=99999-graphallshortestpaths(sparse(Matrix2(Bijection_Trial(:,2),Bijection_Trial(:,2))))(isnan(graphallshortestpaths(sparse(Matrix2(Bijection_Trial(:,2),Bijection_Trial(:,2))))))=99999))) <= Connectivity_Threshold


    %         isequal(Matrix1(Bijection_Trial(:,1),Bijection_Trial(:,1)),Matrix2(Bijection_Trial(:,2),Bijection_Trial(:,2))) == 1
            Bijection = [Bijection; B(i,:)]; %sum(sum(abs(Connectivity1-Connectivity2)))];
            rows = rows + 1;

    %         isequal(Matrix1(Bijection_Trial(:,1),Bijection_Trial(:,1)),Matrix2(Bijection_Trial(:,2),Bijection_Trial(:,2)))
    %         error
        end
        i = i + 1;
    end

    % Add the names back in, so the bijection output tells the user what the
    % mapping is
    C = Bijection;
    [r_C,c_C] = size(C);
    Bijection = [];
    for i = 1:r_C
        Bijection = [Bijection; Ids(C(i,1)) Ids(C(i,2)) C(i,3)]; %C(i,5)];
    end 

    % Arrange the output 
    if Switch == 0
        Bijection = sortrows(Bijection,1);
    else
        B = [B(:,2) B(:,1) B(:,3) B(:,4)];
        B = sortrows(B,1);
        Bijection = [Bijection(:,2) Bijection(:,1) Bijection(:,3)];
        Bijection = sortrows(Bijection,1);
    end
        
        
    % Output the total distance score
%     disp(['Score = ' num2str(sum(cell2mat(Bijection(:,3)))/r_C)])
    Network_Distance = sum(cell2mat(Bijection(:,3)))/r_C;
    

    % Output percent of nodes in the mapping (of the smaller network)
%     disp(['Percent mapping: ' num2str(100*(r_C/Size1)) '% (' num2str(r_C) '/' num2str(Size1) ')'])

    % Determine if an isomorphism was found, including the possibility that
    % there could be one but not for certain when nodes had multiple possible
    % perfect matches
%     if sum(cell2mat(Bijection(:,3))) == 0 & r_C == Size1
%         if sum(B(:,3) == 0) == Size1
%             disp('Isomorphism!')    
%         else
%             disp('Possible Isomorphism...')
%         end
%     else 
%         disp('No isomorphism detected.')
%     end
    
    
    
    
    View_Graph = 0;
    Ids1 = Ids;
    Ids2 = Ids;

    if View_Graph == 1
        Graph = biograph(Matrix1_Initial,Ids1(1:length(Matrix1_Initial)));%,'LayoutType','equilibrium')
        get(Graph.nodes,'ID');
        view(Graph);
        
        Graph = biograph(Matrix2_Initial,Ids2(1:length(Matrix2_Initial)));%,'LayoutType','equilibrium');
        get(Graph.nodes,'ID');
        view(Graph);
    end


%     Have the node labels include the mappings
%     Ids1 = Ids;
%     Ids2 = Ids;
    for i = 1:r_C
        x = Ids1(str2double(Bijection(i,1)));
        y = Ids2(str2double(Bijection(i,2)));

        Ids1(str2double(Bijection(i,1))) = {[x{:}, '      <--> ', y{:}]};    

        Ids2(str2double(Bijection(i,2))) = {[y{:}, '      <--> ', x{:}]};  

    end
    
    
%     View the 2 matrices, with colors denoting the distance from the bijection
    if View_Graph == 1
%         if Switch == 1
%            Bijection = [Bijection(:,2) Bijection(:,1) Bijection(:,3)]; 
%         end
        
        Graph = biograph(Matrix1_Initial,Ids1(1:length(Matrix1_Initial)),'LayoutType','radial');
        get(Graph.nodes,'ID');
        for i = 1:length(Bijection)
            Bijection = sortrows(Bijection,1);
            % Can't divide by zero, which occurs if the distances are all
            % zero
            if max(cell2mat(Bijection(:,3))) == 0
                set(Graph.nodes(str2double(Bijection(i,1))),'Color',[1 1 1]);    
            else
                set(Graph.nodes(str2double(Bijection(i,1))),'Color',[1-cell2mat(Bijection(i,3)) 1-cell2mat(Bijection(i,3)) 1-cell2mat(Bijection(i,3))]);
            end
        end
        view(Graph);
        
        Graph = biograph(Matrix2_Initial,Ids2(1:length(Matrix2_Initial)),'LayoutType','equilibrium');
        get(Graph.nodes,'ID');
        for i = 1:length(Bijection)
            Bijection = sortrows(Bijection,2);
            % Can't divide by zero, which occurs if the distances are all
            % zero
            if max(cell2mat(Bijection(:,3))) == 0
                set(Graph.nodes(str2double(Bijection(i,2))),'Color',[1 1 1]);
            else
                set(Graph.nodes(str2double(Bijection(i,2))),'Color',[1-cell2mat(Bijection(i,3)) 1-cell2mat(Bijection(i,3)) 1-cell2mat(Bijection(i,3))]);        
            end
        end
        view(Graph);
    end

    

    % Output total computation time
%     disp(['Compute time: ' num2str(cputime - Compute_Time) ' seconds'])




    %-------------------
    % Hungarian method


    % [assignment,cost,dMat] = munkres(Similarity)





    % H = mod(find(munkres(Similarity) == 1),length(Similarity))
    % 
    % C = zeros(length(H),3);
    % for i = 1:length(H)
    %     if H(i) == 0
    %         C(i,:) = [i length(H) Similarity(i,length(H))];
    %     else
    %         C(i,:) = [i H(i) Similarity(i,H(i))];
    %     end
    % end
    % 
    % C
    % 
    % [r_C,c_C] = size(C);
    % Bijection = [];
    % for i = 1:r_C
    %     Bijection = [Bijection; Ids(C(i,1)) Ids(C(i,2)) C(i,3)];
    % end 
    % 
    % Bijection   
    % 
    % 
    % 
    % 
    % 
    % 
    % 
    % Matrix1_New = zeros(length(Matrix1));
    % for i = 1:length(Matrix1)
    %     for j = 1:length(Matrix1)
    %         Matrix1_New(i,j) = Matrix1(find(C(:,2)==i),find(C(:,2)==j));
    %     end
    % end
    %    
    % Connectivity = abs(graphallshortestpaths(sparse(Matrix1))-graphallshortestpaths(sparse(Matrix1_New)))











    % Matrix1_New = zeros(length(Matrix1));
    % for i = 1:length(Matrix1)
    %     for j = 1:length(Matrix1)
    %         Matrix1_New(i,j) = Matrix1(find(C(:,2)==i),find(C(:,2)==j));
    %     end
    % end
    %     
    % Bijection
    % 
    % Connectivity = abs(graphallshortestpaths(sparse(Matrix1))-graphallshortestpaths(sparse(Matrix1_New)))

    % % Run it again!
    %     
    % Similarity = Similarity + Connectivity;
    % 
    % H = mod(find(munkres(Similarity) == 1),length(Similarity))
    % 
    % C = zeros(length(H),3);
    % for i = 1:length(H)
    %     if H(i) == 0
    %         C(i,:) = [i length(H) Similarity(i,length(H))];
    %     else
    %         C(i,:) = [i H(i) Similarity(i,H(i))];
    %     end
    % end
    % 
    % C
    % 
    % [r_C,c_C] = size(C);
    % Bijection = [];
    % for i = 1:r_C
    %     Bijection = [Bijection; Ids(C(i,1)) Ids(C(i,2)) C(i,3)];
    % end 
    % 
    % Bijection   
    
% end

% sum(mean(Matrix1_Flows')/Duration1)
% sum(mean(Matrix2_Flows')/Duration2)
% plot(mean(Matrix1_Flows')/Duration1,'red')
% hold on
% plot(mean(Matrix2_Flows')/Duration2)
% hold off

% Network_Distance_Gini = abs(Gini(mean(Matrix1_Flows'))-Gini(mean(Matrix2_Flows')));

% Network_Distance
% Bijection 
% % B

% n(a,:) = [N(a) Network_Distance];
% n(1:a,:)
end