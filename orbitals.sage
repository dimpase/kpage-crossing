# given a list G of permutations on W=[0..n-1], find the orbitals of the 
# group generated by them, i.e. the orbits on WxW. This only needs O(GW^2) operations
def orbitals(G, result="comp"):
  n = len(G[0])
  dO = {(i,j): [(i,j)] for j in xrange(n) for i in xrange(n)}
  O = [[(i,j) for j in xrange(n)] for i in xrange(n)]
  update = True 
  while update:
   update = False
   dOk = copy(dO.keys())
   for i,j in dOk:
      for g in G:
        s, t = g[i], g[j]
        # this does not preserve the natural pairing
        mi, ma = (min(O[i][j], O[s][t]), max(O[i][j], O[s][t]))
        if mi<ma:
           update = True
           dO[mi]+=dO[ma]
           for p,q in dO[ma]: 
              O[p][q] = mi
           dO.pop(ma)
        
  pa = {(i,j): O[j][i] for i, j in dO.keys()} # pairing of orbitals
 
  if result == "comp":
     return dO.keys(),pa
  else: 
     if result == "raw":
        return O
     else: 
        return O, dO, pa 

def orbmats(G): # for testing purposes
  oo=orbitals(G, result="raw")
  rep=list(set(flatten(oo,max_level=1)))
  n=len(G[0])
  A=[zero_matrix(n,n,sparse=True) for i in xrange(len(rep))]
  for k in xrange(len(rep)):
      for i in xrange(n):
          for j in xrange(n):
              A[rep.index(oo[i][j])][i,j]=1
  return A

def printorbitals(d,f):
  keys = sorted(d.keys())
  f.write(str(len(keys))+'\n')
  for i in keys:
    f.write(" "+str(len(d[i])))
  f.write('\n')
  for i in keys:
    for p,q in sorted(d[i]):
      f.write(str(p)+' '+str(q)+'\n')

