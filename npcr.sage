##########################################################################
# gadgets
##########################################################################
def gadgets(n,k):
# blue vertices 0,1,...,n-1 are placed clockwise on a circle
# a red vertex is placed somewhere between i and i+1 (mod n), 
# which is marked by i. 
# Arcs from the red vertex are are coloured in k colours 0..k-1, so that 
# the arcs in the k-th sheet get colour k. This way we get k^n different gadgets
# for each placement of the red vertex. In total we'll have nk^n gadgets.
  return [(x, y) for x in [0..n-1] for y in tuples([0..k-1],n)]

##########################################################################
# automorphisms
# 1. flips: any simultaneous changing of colours (so there are k! of these)
# 2. symmetries of the n-gon, dihedral (so there are 2n of these)
##########################################################################
def gadgets_flip(g,p): # p must be from S_k (type 1 automorphisms)
   return [g.index((x, [p(y+1)-1 for y in s])) for x,s in g]

def gadgets_ia(g,p): # p must be from D_n (type 2 automorphisms)
# there is a subtle bug in the following: 
#   return [g.index((p(x+1)-1,permutation_action(p,s))) for x,s in g]
# namely, the action on the 1st  coordinate of the tuple is wrong!
#
# in fact, this action must be the action on the edges (i,i+1) of the n-gon,
# for indices taken mod n; the edges are labelled by i, i=0..n-1.
# e.g. for n=3 we have the involition (0,2) 
# taking the edge (0,1) to the edge (1,2), and not to (2,0).
    n = len(g[0][1])
    
    def ea(j):
      i1 = p(j+1)-1
      if j<n-1:
         i2 = p(j+2)-1
      else:
         i2 = p(1)-1
      if set([i1,i2])==set([0,n-1]):
         return n-1
      else:  
         return min(i1,i2) 
    
    return [g.index((ea(x),permutation_action(p,s))) for x,s in g]
           

def ggens(g,k):
    H = SymmetricGroup(k)
    G = DihedralGroup(len(g[0][1]))
    gg=[G.gens()[1],G.gens()[1]*G.gens()[0]]
    # print G.gens()
    return [gadgets_flip(g,p) for p in H.gens()]+[gadgets_ia(g,p) for p in gg]

##########################################################################
# Crossing number of two gadgets.
##########################################################################
# helper functions
def sca(s,k): # scaling
   c=[[] for i in xrange(k)]
   for i in xrange(len(s)):
       c[s[i]]+=[3*i]
   return c

def crk(a,b,k):
   i,s0 = a
   j,q0 = b  
   i = 3*i + 1 
   j = 3*j + 1 
   z = zip(sca(s0,k), sca(q0,k))
#
# two gadgets (i,s) and (j,q) represent a 2-page drawing of K_{2,n}
# when i != j, we just count crossings in each hemisphere;
   if i != j:
      return sum([crcount(i,s,j,q) for s,q in z])
# otherwise, we try both possible orderings of i and j and take the minimum
   else:
      return min(sum([crcount(i,   s, i+1, q) for s, q in z]), 
                 sum([crcount(i+1,   s, i, q) for s, q in z])) 

def crcount(i,s,j,q):
#   print i, " ", s, " ", j, " ", q, "\n"
   c = 0
   for x in s:
     for y in q:
       if i<j<x<y or j<x<y<i or x<y<i<j or y<i<j<x or \
          i<y<x<j or y<x<j<i or x<j<i<y or j<i<y<x: 
         c += 1
   return c

def crkmat(n,k):
   g = gadgets(n,k)
   return g,[[crk(a, b, k) for b in g] for a in g]
  
def permmat(p):
   n = len(p)
   def delt(i,j):
      if i==j: 
         return 1
      else:
         return 0
   return matrix(ZZ, n, n, lambda i, j: delt(p[j],i)) #, sparse=True)

def testM(n,k): # test that we have symmetries we should have
   g,M=crkmat(n,k)
   m=matrix(M)
   return [m*permmat(p)==permmat(p)*m for p in ggens(g,k)]

############# output ###########
load orbitals.sage

# find the coefficients of expression of a sum of orbitals          
def mexpress(M,O,test=False):
  r = dict()
  if test==True:
    n = len(M[0])
    for i in xrange(n):
      for j in xrange(n):
        o = O[i][j]
        if o==(i,j):
          if M[i][j] != 0:
             r[(i,j)] = M[i][j]
        else:
          if M[i][j]!=M[o[0]][o[1]]:
             return False,i,j
    return r
  else:
    for i,j in O:
      if M[i][j] != 0:
        r[(i,j)] = M[i][j]
    return r
    
    
def crkprt(n,k,fn):
    f = file(fn,"w")
    g,M=crkmat(n,k)
    oo=orbitals(ggens(g,k), result="c") #"raw")
#    f.write(str(M))
    f.write(str(len(M))+"\n")
    printorbitals(oo[1],f)
    keys=sorted(oo[1].keys())
    e = mexpress(M,keys)
#    print sorted(e.keys())
    f.write(str(len(e.keys()))+'\n')
    for i in sorted(e.keys()):
       f.write(str(keys.index(i))+" "+str(e[i])+'\n')
    f.close()
#    f.write(str(oo))
#    f.write("\n")
#    f.write(str(g))
#    f.write("\n")
#    f.close()

   
