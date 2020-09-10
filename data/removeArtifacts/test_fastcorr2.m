
N = 10000;
M = 1000;

tic;
fprintf('[*] corr2...\n')
z1 = zeros(1,N);
for i = 1:N
    z1(i) = corr2(rand(1,M),rand(1,M));
end
t1 = toc;
fprintf('Elapsed time is %f seconds.\n',t1)

tic;
fprintf('[*] fastcorr2...\n')
z2 = zeros(1,N);
for i = 1:N
    z2(i) = fastcorr2(rand(1,M)',rand(1,M)');
end
t2 = toc;
fprintf('Elapsed time is %f seconds.\n',t2)

fprintf('Time savings factor is %f.\n',t1/t2)

e = zeros(1,N);
for i = 1:N
    a = rand(1,M);
    b = rand(1,M);
    e(i) = sqrt((fastcorr2(a',b')-corr2(a,b))^2);
end
fprintf('[*] average error: %d\n',mean(e))
