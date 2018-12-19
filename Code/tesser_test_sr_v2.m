%TesserScan community memberships by node, last two numbers = bonudary nodes:
comm1_prim = [1, 2, 19, 20, 21]; 
comm1_bound = [3, 18];
connected_bound1 = [3, 4]; 
comm2_prim = [5, 6, 7, 8, 9];
comm2_bound = [4, 10];
connected_bound2 = [10, 11]; 
comm3_prim = [12, 13, 14, 15, 16]; 
comm3_bound = [11, 17];
connected_bound3 = [17, 18]; 

%% I = identity matrix
I = eye(21, 21);

%% A = adjacency matrix, where 1s for every connected node in 21-node temporal community structure 
A = zeros(21, 21);

% connecting primary and boundary nodes of community
A(comm1_prim, comm1_prim) = 1;
A(comm1_prim, comm1_bound) = 1;
A(comm1_bound, comm1_prim) = 1;
A(comm2_prim, comm2_prim) = 1;
A(comm2_prim, comm2_bound) = 1;
A(comm2_bound, comm2_prim) = 1;
A(comm3_prim, comm3_prim) = 1;
A(comm3_prim, comm3_bound) = 1;
A(comm3_bound, comm3_prim) = 1;

% connecting boundary nodes across communities
A(connected_bound1, connected_bound1) = 1;
A(connected_bound2, connected_bound2) = 1;
A(connected_bound3, connected_bound3) = 1;

%make sure that all the diagonals (i.e. items to iteslf) are 0
A(logical(eye(size(A)))) = 0;

%% T = transition probability matrix, where probability to move to every connected node in 21-node temporal community structure 
T = A/6; 

%% Calculating and plotting SR with A matrix 
sr_A1 = inv(I - (.4 * A));
sr_A2 = inv(I - (.8 * A));

figure()
subplot(3,1,1)
imagesc(A); axis square; colorbar
title('A = adjacency matrix')
subplot(3,1,2)
imagesc(sr_A1, [0 .2]); axis square; colorbar
title('SR (gamma=0.4)')
subplot(3,1,3);
imagesc(sr_A2, [0 .2]); axis square; colorbar
title('SR (gamma=0.8)')

%% Calculating and plotting SR with T matrix
sr_T1 = inv(I - (.4 * T));
sr_T2 = inv(I - (.8 * T));

figure()
subplot(3,1,1)
imagesc(T); axis square; colorbar
title('T = transition matrix')
subplot(3,1,2)
imagesc(sr_T1, [0 .2]); axis square; colorbar
title('SR (gamma=0.4)')
subplot(3,1,3);
imagesc(sr_T2, [0 .2]); axis square; colorbar
title('SR (gamma=0.8)')
