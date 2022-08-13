function F = MBBCheck(A,B)
%MBBCHECK 用于检查 GAVE 是否满足 MBB 方法的使用条件
%   MBBCHECK 的原理是比较 A,B 对称部分的特征值，具体可见论文
%   MBBCHECK 的方法是充分条件，即便该检查不通过，MBB 方法也有可能收敛

%   Version 1.0
%   Copyright 2022 S. Yang

 if eigs(A+A',1,'smallestabs') - eigs(B+B',1) > 0
     F = true;
 else
     F = false;
end

