import networkx as nx
import Queue as qx
import math
from scipy.stats.mstats import mquantiles

edge_weight='weight'  # This is the marker for importance of each edge

def spc(G):
	# this will compute non-normalized SPC value in log-scale
	nodes=nx.topological_sort(G)
	l=len(nodes)
	src=source(G)
	snk=sink(G)
	for w in src:
		G.node[w]['Nm']=1
	for w in snk:
		G.node[w]['Np']=1
	for i in range(l):
		j=-i-1
		n1=nodes[i]
		n2=nodes[j]
		s=0
		for w in G.predecessors(n1):
			s=s+G.node[w]['Nm']
		G.node[n1]['Nm']=max(s,1)
		s=0
		for w in G.neighbors(n2):
			s=s+G.node[w]['Np']
		G.node[n2]['Np']=max(s,1)
	for eg in G.edges_iter():
		n1=eg[0]
		n2=eg[1]
		G[n1][n2]['weight']=math.log(G.node[n1]['Nm'])+math.log(G.node[n2]['Np'])
		G.node[n1]['sig']=math.log(G.node[n1]['Nm'])+math.log(G.node[n1]['Np'])
		G.node[n2]['sig']=math.log(G.node[n2]['Nm'])+math.log(G.node[n2]['Np'])
	return G
		
def source(G):
	s=[]
	for w in G.nodes_iter():
		if G.in_degree(w)==0:
			s.append(w)
	return s
		

def sink(G):
	s=[]
	for w in G.nodes_iter():
		if G.out_degree(w)==0:
			s.append(w)
	return s
	
def add_source(G):
	w=source(G)
	if len(w) > 1:
		G.add_node('source')
		for n in w:
			G.add_edge('source',n)
	else:
		print 'The graph already contains single source: '+w[0]
	return G

def remove_source(G):
	w=source(G)
	if len(w) > 1:
		print 'Warning: graph has multiple sources'
		for n in w:
			G.add_node(n)
	else:
		G.remove_node(w[0])
	return G
		
def add_sink(G):
	w=sink(G)
	if len(w) > 1:
		G.add_node('sink')
		for n in w:
			G.add_edge(n,'sink')
	else:
		print 'The graph already contains single sink: '+w[0]
	return G
	
def remove_sink(G):
	w=sink(G)
	if len(w) > 1:
		print 'Warning: graph has multiple sinks'
		for n in w:
			G.add_node(n)
	else:
		G.remove_node(w[0])
	return G
	

def del_indices(ls,ind):
	w=[i for j, i in enumerate(ls) if j not in ind]
	return w
	

def order(v,w):
	a=zip(v,w)
	a.sort()
	l=zip(*a)
	v=list(l[0])
	w=list(l[1])
	return [v,w]
	
def path_weight(G,path,mark):
	s=0.0
	for i in range(len(path)-1):
		w=G[path[i]][path[i+1]][mark]
		s=s+w
	return s

def main_path(*args, **kwds):
	# flags is the list containing desired methods to combine
	# fwd to do forward main path from all nodes in fwd list
	# bwd to do backward main path from all nodes in bwd list

	G=args[0]
	H=nx.DiGraph()
	try:
		fwd=kwds['fwd']
	except KeyError:
		fwd=False
	try:
		bwd=kwds['bwd']
	except KeyError:
		bwd=False
	try:
		gl=kwds['gl']
	except KeyError:
		gl=False

	if not fwd and not bwd and not gl:
		print 'Global main paths is returned as a default'
	if not gl:
		if fwd:
			qf=qx.Queue()
			if type(fwd)!=list:
				fwd=[fwd]
			for node in fwd:
				qf.put(node)

			while not qf.empty():
				s=qf.get(False)
				wei=[]
				ls=[]
				for w in G.successors_iter(s):
					ls.append(w)
					wei.append(-G[s][w][edge_weight])
				if len(ls)==0:continue
				wei,ls=order(wei,ls)
				win=ls[0]
				H.add_node(s,G.node[s])
				H.add_node(win,G.node[win])
				H.add_edge(s,win,G[s][win])
				qf.put(win)
				spc=-wei[0]
				for i,d in enumerate(wei):
					if i==0 :
						i=1
						continue
					d=-d
					if d==spc :
						win=ls[i]
						H.add_node(s,G.node[s])
						H.add_node(win,G.node[win])
						H.add_edge(s,win,G[s][win])
						qf.put(win)
					else:
						break
		if bwd:
			qb=qx.Queue()
			if type(bwd)!=list :
				bwd=[bwd]
			for node in bwd:
				qb.put(node)
			while not qb.empty():
				s=qb.get(False)
				ls=[]
				wei=[]
				for w in G.predecessors_iter(s):
					ls.append(w)
					wei.append(-G[w][s][edge_weight])
				if len(ls)==0:continue
				wei,ls=order(wei,ls)
				win=ls[0]
				H.add_node(s,G.node[s])
				H.add_node(win,G.node[win])
				H.add_edge(win,s,G[win][s])
				qb.put(win)
				spc=-wei[0]
				for i,d in enumerate(wei):
					if i==0 :
						i=1
						continue
					d=-d
					if d==spc :
						win=ls[i]
						H.add_node(s,G.node[s])
						H.add_node(win,G.node[win])
						H.add_edge(win,s,G[win][s])
						qb.put(win)
					else:
						break
		if not fwd and not bwd :
			print 'er'
			tarf=[]
			tarb=[]
			for w in G.__iter__():
				if G.out_degree(w) == 0:
					tarf.append(w)
				if G.in_degree(w) == 0:
					tarb.append(w)	
			for edge in G.edges_iter():
				G[edge[0]][edge[1]]['nwei']=-G[edge[0]][edge[1]][edge_weight]
			met=0
			cur=0
			for n1 in tarb:
				for n2 in tarf:
					pp=nx.shortest_path(G,n1,n2,'nwei')
					met=path_weight(G,pp,'nwei')
					if met < cur:
						cur=met
						path=pp
			for i in range(len(path)-1):
				n1=path[i]
				n2=path[i+1]
				H.add_node(n1,G.node[n1])
				H.add_node(n2,G.node[n2])
				H.add_edge(n1,n2,G[n1][n2])
	else:
		qg=qx.Queue()
		tarf=[]
		tarb=[]
		for w in G.__iter__():
			if G.out_degree(w) == 0:
				tarf.append(w)
			if G.in_degree(w) == 0:
				tarb.append(w)
		if type(gl)==bool:
				path=nx.shortest_path(G,node,tar,weight='nwei')
		else:
			if type(fwd)!=bool or type(bwd)!=bool:
				print ''' ERROR: Your input is inconsistent. 
While you have set gl key to be a node or node list,
fwd and bwd keys can only be boolean.'''
				return H
			if not fwd and not bwd:
				fwd=True
				print 'Default behaviour in global mode is forward search'
			for edge in G.edges_iter():
				G[edge[0]][edge[1]]['nwei']=-G[edge[0]][edge[1]][edge_weight]
			if type(gl)!=list:
				gl=[gl]
			for w in gl:
				qg.put(w)
			while not qg.empty():
				node=qg.get(False)
				if fwd:
					for tar in tarf:
						path=nx.shortest_path(G,node,tar,weight='nwei')
						for i in range(len(path)-1):
							n1=path[i]
							n2=path[i+1]
							H.add_node(n1,G.node[n1])
							H.add_node(n2,G.node[n2])
							H.add_edge(n1,n2,G[n1][n2])
				if bwd:
					for tar in tarb:
						path=nx.shortest_path(G,tar,node,weight='nwei')
						for i in range(len(path)-1):
							n1=path[i]
							n2=path[i+1]
							H.add_node(n1,G.node[n1])
							H.add_node(n2,G.node[n2])
							H.add_edge(n1,n2,G[n1][n2])
	return H

def key_route(G):
	temp=0
	cur=0
	for w in G.edges_iter():
		can=G[w[0]][w[1]][edge_weight]
		if can>cur:
			cur=can
			win=w
	n1=win[0]
	n2=win[1]
	H=main_path(G,fwd=n2,bwd=n1)
	H.add_edge(n1,n2,G[n1][n2])
	return H
	
def phyl(G):
	H=nx.DiGraph()
	so=source(G)
	si=sink(G)
	print 'Walking forward/backward...'
	Hfb=main_path(G,fwd=so,bwd=si)
	print 'Walking global...'
	Hg=main_path(G,gl=so,fwd=True)
	print 'Key-Route...'
	Hk=key_route(G)
	print 'Combining...'
	S=set(Hfb.edges()+Hg.edges()+Hk.edges())
	print 'MP # of edges: '+str(len(S))
	for ed in S:
		n1=ed[0]
		n2=ed[1]
		H.add_node(n1,G.node[n1])
		H.add_node(n2,G.node[n2])
		H.add_edge(n1,n2,G[n1][n2])
	return H
	
def prune(B,H=None,per=.2,deg=True,cap='in'):
	# prune Graph based on degree or properties in cap
	# node based
	# H will always be preserved
	G=nx.DiGraph(B)
	if deg:
		if cap=='in':
			seq=G.in_degree().values()
			cut=mquantiles(seq,1-per)
			mk=1
		elif cap=='out':
			seq=G.out_degree().values()
			cut=mquantiles(seq,1-per)
			mk=2
		else:
			print 'Error in setting cap string'
			return None
	else:
		seq=[]
		S=zip(*G.nodes(True))[1]
		for w in S:
			try:
				seq.append(w[cap])
				mk=3
			except KeyError:
				print 'Some nodes lack cap string'
				return None
		cut=mquantiles(seq,1-per)
	to_del=[]
	for n in G.nodes_iter():
		G.node[n]['size']=2*G.in_degree(n)
		if H is not None:
			if n in H.nodes():
				continue
		if mk==1:
			t=G.in_degree(n)
			if t < cut:
				to_del.append(n)
		elif mk==2:
			t=G.out_degree(n)
			if t < cut:
				to_del.append(n)
		else:
			t=G.node[n][cap]
			if t < cut:
				to_del.append(n)
	G.remove_nodes_from(to_del)
	return G
	
	
	
			
			

	
	
