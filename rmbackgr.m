function out = rmbackgr(input)
%
% ±³¾°È¥³ý

    [m,n] = size(input);
    out = input - mean(input' )' * ones( 1, n );
return
