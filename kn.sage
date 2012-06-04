def kncr(n):
  g=graphs.EmptyGraph()
  chords=[]
  for i in range(n-1):
     for j in range(i+2,n):
       if abs(i-j)>1 and abs(i-j)<n-1:
          chords.append((i,j))
  g.add_vertices(chords)
  edges=[]
  for x in chords:
    for y in chords:
      if x!=y:
        if x[0]<y[0]<x[1]<y[1]: 
           edges.append((x,y))
  g.add_edges(edges) 
  return g

#Use it e.g. as follows:

#sage: load 'kn.sage'
#sage: g=kncr(5)
#sage: g.max_cut()
#4.0

