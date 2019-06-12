%TesserScan community memberships by node, last two numbers = bonudary nodes:
comm1 = [1, 2, 19, 20, 21, 3, 18]; 
comm2 = [5, 6, 7, 8, 9, 4, 10];
comm3 = [12, 13, 14, 15, 16, 11, 17]; 

% everything is adjacent to itself
A = eye(21, 21);

% internal
A(2:6,2:6) = 1;
A(9:13,9:13) = 1;
A(16:20,16:20) = 1;

% boundary
A(1,2:6) = 1;
A(2:6,1) = 1;
A(1,end) = 1;
A(end,1) = 1;

A(7,2:8) = 1;
A(2:8,7) = 1;

A(8,9:13) = 1;
A(9:13,8) = 1;

A(14,9:15) = 1;
A(9:15,14) = 1;

A(15,16:20) = 1;
A(16:20,15) = 1;

A(21,16:20) = 1;
A(16:20,21) = 1;

%imagesc(A); axis square

T = A/7; % transition probability

sr1 = inv(eye(21) - .4 * T);
sr2 = inv(eye(21) - .8 * T);

close all
subplot(1,3,1)
imagesc(A); axis square; colorbar
title('adjacency')
subplot(1,3,2)
imagesc(sr1, [0 .2]); axis square; colorbar
title('SR (gamma=0.4)')
subplot(1,3,3);
imagesc(sr2, [0 .2]); axis square; colorbar
title('SR (gamma=0.8)')
