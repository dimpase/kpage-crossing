##########################################################################
# gadgets
##########################################################################
def gadgets(n,k):
# blue vertices 0,1,...,n-1 are placed clockwise on a circle
# a red vertex is placed somewhere between i and i+1 (mod n), 
# which is marked by i. 
# Arcs from the red vertex are are coloured in k colours 0..k-1, so that 
# the arcs in the k-th sheet get colour k. This way we get k^n different gadgets
# each placement of the red vertex. In total we'll have nk^n gadgets.
  return [(x, y) for x in [0..n-1] for y in tuples([0..k-1],n)]

##########################################################################
# automorphisms
# 1. flips: any simultaneous changing of colours (so there are k! of these)
# 2. symmetries of the n-gon, dihedral (so there are 2n of these)
##########################################################################
def gadgets_flip(g,p): # p must be from S_k (type 1 automorphisms)
   return [g.index((x, [p(y+1)-1 for y in s])) for x,s in g]

def gadgets_ia(g,p): # p must be from D_n (type 2 automorphisms)
   return [g.index((p(x+1)-1,permutation_action(p,s))) for x,s in g]

def ggens(g,k):
    H = SymmetricGroup(k)
    G = DihedralGroup(len(g[0][1]))
    gg=[G.gens()[1],G.gens()[1]*G.gens()[0]]
    # print G.gens()
    return [gadgets_flip(g,p) for p in H.gens()]+[gadgets_ia(g,p) for p in gg]

##########################################################################
# Crossing number of two gadgets.
##########################################################################
def crk(a,b):
   i,s = a
   j,q = b  
   if i>j: 
     return crk(b,a)
   if i==j:
     return min(crk0(i,s,j,q), crk0(i,q,j,s)) 
   return crk0(i,s,j,q)

# we assume here that j is to the left from i (when i=j)
# there are 2 types of possible crossings of ir and jt
#  a) i<t<r<j, where '<' has the usual meaning
#  b) j<r<t<i, where '<' means the order in the segment (j,j+1,...,n-1,0,1,...,i)
def crk0(i,q,j,s): 
   n = len(s)
   c = 0 # of crossings
   for t in [i+1..j-2]: # type a)
     for r in [t+1..j-1]:
        if s[r] == q[t]: c+=1
   for r in [j+1..i+n-2]: # type b)
     for t in [r+1..i+n-1]:
        if s[r % n] == q[t % n]: c+=1
   return c   

def crkmat(n,k):
   g = gadgets(n,k)
   return g,[[crk(a, b) for b in g] for a in g]
  
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
def crkprt(n,k,fn):
    f = file(fn,"w")
    g,M=crkmat(n,k)
    oo=orbitals(ggens(g,k), result="raw")
    f.write(str(M))
    f.write("\n")
    f.write(str(oo))
    f.write("\n")
    f.write(str(g))
    f.write("\n")
    f.close()

   
