load cutgraph.sage
k45=graphs.CompleteBipartiteGraph(4,5)
a=axgr(k45,[0..8])
print a.clique_maximum()
