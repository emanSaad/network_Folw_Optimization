%% read the network from file
% this snippet of code from jan:stackoverflow
out = [];
len = 20;
part = zeros(2,len);
ipart = 0;
fid = fopen('graph20N.txt', 'r');
if fid < 0, error('Cannot open file'); end
s1 = fgets(fid);
  if ischar(s1)
    data = sscanf(s1, '%g %g', 2);
    numberOfNodes=data(1,1);
    numberOfLinks=data(2,1);
  end
while 1  % Infinite loop
  s = fgets(fid);
  if ischar(s)
    data = sscanf(s, '%g %g', 2);
    
    if length(data) == 2
      ipart = ipart + 1;
      part(:, ipart) = data;
      if ipart == len
        out = cat(2, out, part);
        ipart = 0;
      end
    end
  else  % End of file:
    break;
  end
end

out = cat(2, out, part(:, 1:ipart));
network=out';

%%
