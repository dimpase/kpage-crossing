from sage.rings.real_double import RDF
from sage.matrix.constructor import matrix
from sage.modules.free_module_element import vector
from sage.all import * 

##########################################################################
# computing the Goemans-Williamson sdp upper bound for 
# MAXCUT on an (edge-weighted) graph specified by its adjacency matrix A
#
# we maximize Tr(AX) subject to X psd, and X_{ii}=1 for all i=1..d
#   
##########################################################################
def maxcutrel(A, solver=None): 
   from cvxopt.base import matrix as m
   from cvxopt.base import spmatrix as spm
   from cvxopt import solvers
   if solver=='dsdp': 
      from cvxopt import dsdp 
   solvers.options['show_progress']=True
#   solvers.options['show_progress']=False

# data type and dimensions
   try: 
      A = A.weighted_adjacency_matrix()
      return maxcutrel(A, solver=solver)
   except (AttributeError,TypeError):
      try:
          A = A.adjacency_matrix()
          return maxcutrel(A, solver=solver)
      except (AttributeError,TypeError):
          pass
      pass
   try: 
      d = A.nrows()

# SDP constraints 
      c = m([-1.]*d)
      G = [spm(1., [i*(d+1) for i in range(d)], range(d), (d*d,d))]
   
      sol = solvers.sdp(c, Gs=G, hs=[m(1.*A.numpy())], solver=solver)
      return 0.25*(sum(sum(A))+sol['primal objective']),sol
   #return sol['primal objective'], sol['status']
   except AttributeError:
      print "Wrong matrix format"
      raise
   else:
      print "Unknown error"
      raise

##########################################################################
# generating NEOS BiqMac MAXCUT solver input
#  http://www.neos-server.org/neos/solvers/co:BiqMac/SPARSE.html
#
#  An example input file is at:
#  http://www.neos-server.org/neos/solvers/co:BiqMac/Mac_SPARSE.txt
# 
##########################################################################

def genbiqmacinput(g,fname="/Users/dima/Desktop/blah"):
   f = open(fname, "w")
# file header
   f.write(str(g.num_verts())+" "+str(g.num_edges())+"\n")

   v = dict()
   x = g.vertex_iterator()
   i = 0
   try: 
     while (1):
        i += 1
        v[x.next()] = i
   except StopIteration:
     pass
   
   x = g.edge_iterator()
   try: 
     while (1):
        e = x.next()
        f.write(str(v[e[0]])+" "+str(v[e[1]])+" 1\n")
   except StopIteration:
     pass
   f.close()

