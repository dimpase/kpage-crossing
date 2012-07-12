#from sage.rings.real_double import RDF
#from sage.matrix.constructor import matrix
#from sage.modules.free_module_element import vector
#from sage.all import * 

##########################################################################
# gadgets
##########################################################################
def gadgets(n):
# blue vertices 0,1,...,n-1 are placed clockwise on a circle
# a red vertex is placed somewhere between i and i+1 (mod n), 
# which is marked by i. Arcs in the northern hemisphere 
# from the red vertex are specified by a subset s of [0..n-1]
# the arcs in the southern hemisphere are implicit: [0..n-1]\s
#   
#   return reduce(lambda a,b: a+b, 
#                 [[(i,s) for s in reduce(lambda z, x: z + [y + [x] for y in z], 
#                                         range(n), [[]])] 
#                  for i in range(n)])
   return reduce(lambda a,b: a+b, 
                 [[(i,s)  for i in range(n)] 
                        for s in reduce(lambda z, x: z + [y + [x] for y in z], 
                                         range(n), [[]])])

##########################################################################
# automorphisms
# 1. flip: (x,s)->(x,[0..n-1]\s)
##########################################################################
def gadgets_flip(n,g):
   return [g.index((x,sorted(list(set(range(n)).difference(s))))) for x,s in g]

##########################################################################
# 2. indiced action of p in D_{2n}
#     p is given is a list of images of [0..n-1]
##########################################################################
def gadgets_indact(n,g,p):
    return [g.index((p[x],sorted(list(set([p[y] for y in s]))))) for x,s in g]
#    return [g.index((p[(x+1)%n],list(set([p[y] for y in s])))) for x,s in g]

##########################################################################
# rearrange gadgets to make the order block-cyclic
##########################################################################

def c_ord(n,g):
   gg0 =  PermutationGroup([[y+1 for y in # '+1' to make GAP happy...
                             gadgets_indact(n,g,[(x+1)%n for x in range(n)])],
                           [y+1 for y in gadgets_flip(n,g)]])
   print gg0.gen(1)*gg0.gen(0)==gg0.gen(0)*gg0.gen(1)
   print gg0.gen(0),gg0.gen(1)
   gg = PermutationGroup([gg0.gen(1)*gg0.gen(0)]) # would give the biggest
                                                  # block-cyclic structure
   nord = []
   olens = []
   for x in gg.orbits(): 
      nord += x
      olens.append(len(x))
   return olens,[g[x-1] for x in nord]

##########################################################################
# the group generators
##########################################################################
def gadgets_group(n,g): # Does it work?!
    Dn = DihedralGroup(n).gens()
    x0 = [x-1 for x in Dn[0].list()]
    x1 = [x-1 for x in Dn[1].list()]
    return PermutationGroup(
           [Permutation([x+1 for x in gadgets_flip(n,g)]),
            Permutation([x+1 for x in gadgets_indact(n,g,x0)]),
            Permutation([x+1 for x in gadgets_indact(n,g,x1)])])

##########################################################################
# Crossing number of from two gadgets.
##########################################################################
def cr2(n,a,b):
   i,s = a
   j,q = b  
# to make place to accommodate both a and b, we scale n using sca()
   s = sca(s)
   q = sca(q)
   cs = sorted(list(set(sca(range(n))).difference(s)))
   cq = sorted(list(set(sca(range(n))).difference(q)))
   i = 3*i + 1 
   j = 3*j + 1 
#
# two gadgets (i,s) and (j,q) represent a 2-page drawing of K_{2,n}
# when i != j, we just count crossings in each hemisphere;
   if i != j:
      return crcount(i,s,j,q)+crcount(i,cs,j,cq)
# otherwise, we try both possible orderings of i and j and take the minimum
   else:
      return min(crcount(i,   s, i+1, q)+crcount(i,   cs, i+1, cq),
                 crcount(i+1, s, i,   q)+crcount(i+1, cs, i,   cq))


# helper functions
def sca(s): return [3*x for x in s]  # scaling
 
def crcount(i,s,j,q):
#   print i, " ", s, " ", j, " ", q, "\n"
   c = 0
   for x in s:
     for y in q:
#       print i, j, x, y, "\n"
       if i<j<x<y or j<x<y<i or x<y<i<j or y<i<j<x or \
          i<y<x<j or y<x<j<i or x<j<i<y or j<i<y<x: 
         c += 1
   return c

def permmat(p):
   n = len(p)
   def delt(i,j):
      if i==j: 
         return 1
      else:
         return 0
   return matrix(QQ, n, n, lambda i, j: delt(p[j],i), sparse=True)

def cr2mat(n, raw = False):
   olens,g = c_ord(n,gadgets(n))
   reps = []
   pos = 0
   for x in olens:
      reps.append(g[pos])
      pos = pos + x
   if raw:
      reps = g
   return olens,g,[[cr2(n, a, b) for b in g] for a in reps]
  
def cr2prt(n,fn):
    f = file(fn,"w")
    olens,a,b=cr2mat(n)
    f.write("gadgets=\n")
    f.write(str(a))
    f.write("\n matrix=\n")
    f.write(str(b))
    f.write("\n blocksizes=\n")
    f.write(str(olens))
    f.close()
