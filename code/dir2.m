function listing = dir2(name)

listing = dir(name);

to_keep = ones(1, length(listing));
inds = [];

for i = 1:length(listing)
    current_file = listing(i).name;
    if current_file(1) == '.'
        inds(end + 1) = i;
    end
end
listing(inds) = [];
end