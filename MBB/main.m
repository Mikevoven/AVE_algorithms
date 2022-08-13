load data.mat
%% 
A = GAVETestData.data3.size_3600.A;
B = GAVETestData.data3.size_3600.B;
b = GAVETestData.data3.size_3600.b;
%% MBB
tic
[~,i,err] = MBB(A,B,b);
toc
disp(i)
disp(err)
%% SOR-like
%[miniter, opt_w] = SOR_like_w(A,B,b)

w = 0.8;
tic
[~,i,err] = SOR_like(A,B,b,w);
toc
disp(i)
disp(err)

%% GN
tic
[~,i,err] = GN(A,B,b);
toc
disp(i)
disp(err)

%% SignAccord
tic
[x,S,err,iter,flag] = SignAccord(A, B, b);
toc
disp(i)
disp(err)