function outName = createUniqueFile(basePath, ext)
% createUniqueFile  Create a new empty file ensuring no overwrite.
% outName = createUniqueFile(basePath, ext)
% basePath: full path without extension or folder+basename
% ext: extension with leading dot, e.g. '.txt'

if nargin<2, ext = '.txt'; end
% start with the given name
candidate = [basePath ext];
i = 1;
while isfile(candidate)
    candidate = sprintf('%s_%d%s', basePath, i, ext);
    i = i + 1;
end
% create the file (empty)
fid = fopen(candidate, 'w');
if fid == -1
    error('Could not create file: %s', candidate);
end
fclose(fid);

outName = candidate;
end
