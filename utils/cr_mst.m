function cable = cr_mst(layout_coords, wf)

n = size(layout_coords,1);
sub = wf.substation;

coords = [layout_coords; sub];

D = squareform(pdist(coords));
T = minspantree(graph(D));

[s,t] = findedge(T);
len = D(sub2ind(size(D),s,t));

s(s==n+1) = 0;
t(t==n+1) = 0;

cable.connections = [s t];
cable.lengths = len;
cable.types = ones(numel(len),1);

if isfield(wf,'innercable_price') && ~isempty(wf.innercable_price)
    price = wf.innercable_price(end);
    cable.total_cost = sum(len) * price;
else
    cable.total_cost = NaN;
end

cable.MST = true;

end