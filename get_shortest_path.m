function SP = get_shortest_path(u, v, pred)
% FUNCTION: get the nodes in the shortest path
% INPUT:
% u is the head node;
% v is the tail node;
% pred is the father node;
% OUTPUT:
% SP is the shortest path;

SP = v;
while pred(v)~=0
    SP = [pred(v), SP]; 
    v = pred(v);
end
    
end