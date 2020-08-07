function [ Y ] = classify_v2( X )

x1_t1 = 10; % uV
x1_t2 = 2000; % uV
x2_t = 400; % uV
% ratio [0,1]
%x3_t = 0.6;
%x4_t = 0.4;

if (X(1) < x1_t1)
    Y = 1;
elseif (X(1) > x1_t2)
    Y = 2;
elseif (X(2) > x2_t)
    Y = 3;
% elseif (X(3) > x3_t)
%     Y = 4;
% elseif (X(4) > x4_t)
%     Y = 5;
else
    Y = 0;
end

end

