function tag = getTag(str, delimiter)
%
% function tag = getTag(str, delimiter)
%
% the tag is the string that comes after the last delimiter and before the
% dot in the file extension
%
% Rachel Denison
% 2014

if nargin < 2
    delimiter = '_';
end

delimIdx = strfind(str, delimiter);
dotIdx = strfind(str, '.');
tag = str(delimIdx(end)+1:dotIdx(end)-1);