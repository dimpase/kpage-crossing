# given an undirected graph and an ordering of vertices, set up the corresponding
# subgraph of the complement of the line graph, to do its max_cut, or max-k-cut
#

def axgr(h,p):

  def opa(x):
    x0 = p.index(x[0])
    x1 = p.index(x[1])
    if x0 > x1:
       x0, x1 = x1, x0
    return x0,x1

  g=graphs.EmptyGraph()
  e=h.edges()
  g.add_vertices(e)
  for x in e: 
    x0, x1 = opa(x) 
    for y in e: 
       y0, y1 = opa(y) 
       if x0<y0<x1<y1 or y0<x0<y1<x1:
           g.add_edge(x,y)
  return g 
       
def gray_code(n):
  if n<2:
    return [[0],[1]]
  x = gray_code(n-1)
  return [[0]+y for y in x]+[[1]+y for y in reversed(x)]      

def axrg_hypercube(n):
  return axgr(graphs.CubeGraph(n), 
              [reduce(lambda x,y: str(x)+str(y),t) for t in gray_code(n)])

def hyprebound(n):
   return (5/32)*4^n-2^(n-2)*floor((n^2+1)/2)
  
# this is the old code that does it for K_n
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

