function design = conditions2design(conditions)
% Convert a vector of condition numbers to a binary design matrix
%
% design = conditions2design(conditions)
nonblanks = find(conditions);
sz = [length(conditions), max(conditions)];
linearInd = sub2ind(sz, (nonblanks), conditions(nonblanks));

design = zeros(sz);
design(linearInd) = 1;




