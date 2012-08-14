confs := function(m,n)
  local g;
  g:=DihedralGroup(IsPermGroup, 2*(m+n));
  return List(Orbits(g,Combinations([1..m+n],n),OnSets),x->x[1]);
end;

prtconfs := function(m,n,f)
  local v, x, c;
  c:=confs(m,n);
  PrintTo(f, m, " ", n, " ", Length(c), "\n");
  for v in c do
     for x in v do
       AppendTo(f,x, " ");
     od;
     AppendTo(f,"\n");
  od;
end;
