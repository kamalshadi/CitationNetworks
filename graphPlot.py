import networkx as nx
import pylab as pl

#position given scenario "x coordinate" and "y coordinate"
# label color for nodes : "Interior color"
def order(v,w):
	a=zip(v,w)
	a.sort()
	l=zip(*a)
	v=list(l[0])
	w=list(l[1])
	return [v,w]
	
def myDraw(G,posF=False,H=None,labels=None):
	#It will draw a graph considering following parameters
	# weight for edge
	# label/color/size for nodes
	# Note all these parameters are optional
	color="color"
	weight="weight1"
	alpha='alpha'
	aD=0.1
	dS=100
	dC='blue'
	Ew=1
	pos={}
	aw=False
	jum=True
	if posF:
		if type(posF)==dict:
			pos=posF
			jum=False
	for w in G.nodes(True):
		if posF and jum:
			try:
				pos[w[0]]=(G.node[w[0]]['x coordinate'],G.node[w[0]]['y coordinate'])
			except KeyError:
				G.remove_node(w[0])
				continue
		if color not in w[1].keys():
			G.node[w[0]][color]=dC
		if 'size' not in w[1].keys():
			G.node[w[0]]['size']=dS
		if alpha not in w[1].keys():
			G.node[w[0]][alpha]=aD
	if G.order()==0:
		return
	if not posF:
		pos=nx.graphviz_layout(G,prog='dot')
	if H is None:
		N=nx.get_node_attributes(G,color)
		nC=N.values()
		z=N.keys()
		z,nC=order(z,nC)
		N=nx.get_node_attributes(G,'size')
		nS=N.values()
		z=N.keys()
		z,nS=order(z,nS)
		nx.draw_networkx_nodes(G, pos,nodelist=z,node_size=nS,node_color=nC)
		if labels:
			t=G.nodes(True)
			a1=[w[0] for w in t ]
			a2=[w[1][labels] if labels in w[1].keys() else w[0] for w in t]
			nx.draw_networkx_labels(G,pos, dict(zip(a1,a2)),font_size=24)
		elist=[]
		ew=[]
		aw=[]
		for w in G.edges_iter():
			elist.append(w)
			if weight not in G[w[0]][w[1]]:
				ew.append(Ew)
			else:
				ew.append(G[w[0]][w[1]][weight])
		nx.draw_networkx_edges(G, pos, edgelist=elist, width=ew,arrows=aw)
	else:
		for nw in G.nodes():
			nC=G.node[nw][color]
			nS=G.node[nw]['size']
			try:
				unused=H[nw]
				nA=1
			except KeyError:
				nA=aD
			nx.draw_networkx_nodes(G, pos,nodelist=[nw],node_size=[nS],node_color=[nC],alpha=nA)
			if labels:
				a1=[nw]
				a2=[G.node[nw][labels] if labels in G.node[nw].keys() else a1[0]]
				nx.draw_networkx_labels(G,pos, dict(zip(a1,a2)),alpha=0)
		elist=[]
		ew=[]
		aw=[]
		for w in G.edges_iter():
			elist=[w]
			if weight not in G[w[0]][w[1]]:
				ew.append(Ew)
			else:
				ew.append(G[w[0]][w[1]][weight])
			try:
				unused=H[w[0]][w[1]]
				anow=1
			except KeyError:
				anow=aD
			nx.draw_networkx_edges(G, pos, edgelist=elist, width=ew,alpha=anow,arrows=aw)
	ax = pl.gca()
	ax.yaxis.set_visible(False)
	ax.xaxis.set_visible(False)
	
def myLayout(G,l):
	dic={}
	for nod in G.nodes_iter():
		si=str(G.node[nod][l])
		try:
			dic[si].append(nod)
		except KeyError:
			dic[si]=[nod]
	levels=sorted(dic.keys())
	y=0
	xs=10.0
	pos={}
	for e in levels:
		ls=dic[e]
		m=len(ls)
		dx=xs/m
		seq=G.degree(ls).values()
		un,p=order(seq,ls)
		for i,n in enumerate(p):
			if i==0:
				x=0
			elif i%2==0:
				x=(-i/2)*dx
			else:
				x=((i/2)+1)*dx
			pos[n]=(x,y)
		y=y+1
	return pos
		
		
		
