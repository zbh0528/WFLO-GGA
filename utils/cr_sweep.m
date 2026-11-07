function cable = cr_sweep(layout_coords, wf)

sub = wf.substation;
n = size(layout_coords,1);

if isfield(wf,'cable_capacity') && ~isempty(wf.cable_capacity)
    cap = wf.cable_capacity(end);
else
    cap = 5;
end

if n <= cap
    [conn,len] = build_group_mst((1:n)', layout_coords, sub);
    cable.connections = conn;
    cable.lengths = len;
    cable = assign_cable_types(cable, n, wf.cable_capacity);
    if isfield(wf,'innercable_price') && ~isempty(wf.innercable_price)
        cable.total_cost = sum(cable.lengths .* wf.innercable_price(cable.types)');
    end
    return
end

v = layout_coords - sub;
ang = mod(atan2(v(:,2),v(:,1)),2*pi);
[~,ord] = sort(ang);
order = (1:n)';
order = order(ord);

groups = {};
cur = [];

for k = 1:n
    cur(end+1,1) = order(k);
    if numel(cur) == cap
        groups{end+1,1} = cur;
        cur = [];
    end
end

if ~isempty(cur)
    groups{end+1,1} = cur;
end

conn_all = [];
len_all = [];

for g = 1:numel(groups)
    gidx = groups{g};
    [conn,len] = build_group_mst(gidx, layout_coords, sub);
    conn_all = [conn_all; conn];
    len_all  = [len_all; len];
end

cable.connections = conn_all;
cable.lengths = len_all;

cable = assign_cable_types(cable, n, wf.cable_capacity);

if isfield(wf,'innercable_price') && ~isempty(wf.innercable_price)
    cable.total_cost = sum(cable.lengths .* wf.innercable_price(cable.types)');
end

end


function [connections,lengths] = build_group_mst(idx, layout_coords, sub)

coords = layout_coords(idx,:);
m = size(coords,1);

if m == 1
    connections = [0 idx];
    lengths = vecnorm(coords - sub,2,2);
    return
end

D = squareform(pdist(coords));
T = minspantree(graph(D));

[s,t] = findedge(T);
conn1 = [idx(s) idx(t)];
len1  = D(sub2ind(size(D),s,t));

dsub = vecnorm(coords - sub,2,2);
[~,j] = min(dsub);

connections = [conn1; 0 idx(j)];
lengths = [len1(:); dsub(j)];

end


function cable = assign_cable_types(cable, n_turb, cap)

conn = cable.connections;
conn(conn==0) = n_turb+1;

G = graph(conn(:,1),conn(:,2));
[edges,~] = dfsearch(G,n_turb+1,'edgetonew');

T = digraph(edges(:,1),edges(:,2));

n = n_turb+1;
subtree = ones(n,1);

order = flip(dfsearch(G,n_turb+1));

for i = 1:length(order)
    u = order(i);
    ch = successors(T,u);
    for c = ch'
        subtree(u) = subtree(u) + subtree(c);
    end
end

types = zeros(size(conn,1),1);

for i = 1:size(conn,1)

    u = conn(i,1);
    v = conn(i,2);

    if ismember(v,successors(T,u))
        child = v;
    else
        child = u;
    end

    load = subtree(child);

    if load <= cap(1)
        types(i) = 1;
    elseif load <= cap(2)
        types(i) = 2;
    elseif load <= cap(3)
        types(i) = 3;
    else
        types(i) = -1;
    end
end

cable.types = types;

end