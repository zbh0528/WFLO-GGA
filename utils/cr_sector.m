function cable = cr_sector(layout_coords, wf)

sub = wf.substation;

n = size(layout_coords,1);
cap = wf.cable_capacity(end);

idx = (1:n)';

if n <= cap
    [conn,len] = build_group_mst(idx, layout_coords, sub);
    cable.connections = conn;
    cable.lengths = len;
    cable = assign_cable_types(cable, n, wf.cable_capacity);
    cable.total_cost = sum(cable.lengths .* wf.innercable_price(cable.types)');
    return
end

vec = layout_coords - sub;
ang = mod(atan2(vec(:,2),vec(:,1)),2*pi);
[~,ord] = sort(ang);
idx_sorted = idx(ord);

K = ceil(n/cap);

best_cost = Inf;

for s = 1:cap

    order = [idx_sorted(s:end); idx_sorted(1:s-1)];
    split = round(linspace(0,n,K+1));

    conn_all = [];
    len_all = [];

    for g = 1:K
        i1 = split(g)+1;
        i2 = split(g+1);
        group = order(i1:i2);

        if isempty(group)
            continue
        end

        [conn,len] = build_group_mst(group, layout_coords, sub);
        conn_all = [conn_all; conn];
        len_all = [len_all; len];
    end

    cable_tmp.connections = conn_all;
    cable_tmp.lengths = len_all;

    cable_tmp = assign_cable_types(cable_tmp, n, wf.cable_capacity);

    cost = sum(cable_tmp.lengths .* wf.innercable_price(cable_tmp.types)');

    if cost < best_cost
        best_cost = cost;
        best_conn = cable_tmp.connections;
        best_len  = cable_tmp.lengths;
        best_type = cable_tmp.types;
    end
end

cable.connections = best_conn;
cable.lengths = best_len;
cable.types = best_type;
cable.total_cost = best_cost;

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
len1 = D(sub2ind(size(D),s,t));

dist = vecnorm(coords - sub,2,2);
[~,k] = min(dist);

connections = [conn1; 0 idx(k)];
lengths = [len1(:); dist(k)];

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
        types(i)=1;
    elseif load <= cap(2)
        types(i)=2;
    elseif load <= cap(3)
        types(i)=3;
    else
        types(i)=-1;
    end
end

cable.types = types;

end