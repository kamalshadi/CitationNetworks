# this function return the phenolegy graph 
import networkx as nx
import Queue as qx
from main_path import *
import pylab as pl
from itertools import combinations

level0='epoch'

def del_indices(ls,ind):
	w=[i for j, i in enumerate(ls) if j not in ind]
	return w
	

def phenolegy(G,similarity=False):
	H=nx.DiGraph()
	dic={}
	w=source(G)
	if similarity:
		if len(w)>1:
			print 'Currently similarity could be set only for single source matrices'
	l1=set(nx.neighbors(G,w[0])+G.predecessors(w[0])+w)
	for node in G.nodes_iter():
		key=G.node[node][level0]
		sig=G.node[node]['sig']
		if similarity:
			l2=set(nx.neighbors(G,node)+G.predecessors(node)+[node])
			sim=jaccard(l1,l2)
		try:
			if similarity:
				dic[key].append((sig,sim,node))
			else:
				dic[key].append((sig,node))
		except KeyError:
			if similarity:
				dic[key]=[(sig,sim,node)]
			else:
				dic[key]=[(sig,node)]
	l0=min(dic.keys())
	le=max(dic.keys())
	K=le-l0+1
	print 'Detected epoches: '+str(K) 
	ls=[]
	win=[]
	for l in range(l0,le+1):
		if similarity:
			a=list(zip(*sorted(dic[l],reverse=True))[-1])
			b=list(zip(*sorted(dic[l],reverse=True,key=lambda tup:tup[1]))[-1])
			tempx=two_sort(a,b)
			ls.append(tempx)
			win.append([tempx[0]])
		else:
			tempx=list(zip(*sorted(dic[l],reverse=True))[-1])
			ls.append(tempx)
			win.append([tempx[0]])
	# dic shows the chosen numbers of
	state=dict(zip(*[range(l0,le+1),[1]*K]))
	fg=True
	while fg:
		upward=True
		downward=True
		for l in range(l0,le):
			nb=win[l]
			nu=win[l+1]
			for wb in nb:
				H.add_node(wb,G.node[wb])
				for wu in nu:
					try:
						H.add_edge(wu,wb,G[wu][wb])
					except KeyError:
						if H.in_degree(wb)==0 and nu[-1]==wu:
							try:
								j=state[l+1]
								state[l+1]=state[l+1]+1
								win[l+1].append(ls[l+1][j])
								upward=False
							except IndexError:
								pass
		print 'Upward once'
		for l in reversed(range(l0+1,le+1)):
			nb=win[l-1]
			nu=win[l]
			for wu in nu:
				H.add_node(wu,G.node[wu])
				for wb in nb:
					try:
						H.add_edge(wu,wb,G[wu][wb])
					except KeyError:
						if H.out_degree(wu)==0 and nb[-1]==wb:
							try:
								j=state[l-1]
								state[l-1]=state[l-1]+1
								win[l-1].append(ls[l-1][j])
								downward=False
							except IndexError:
								pass
		if upward and downward:
		#~ if upward:
			break
	return H

def subgraph(G,node):
	print G.node[node]
	bunch=[]
	for w in G.nodes_iter():
		if w==node:
			continue
		if nx.has_path(G,node,w):
			bunch.append(w)
	S1=G.subgraph(bunch+[node])
	S=nx.DiGraph(S1)
	S.node[node]['color']='red'
	return S
	
def jaccard(l1,l2):
	s1=set(l1)
	s2=set(l2)
	n1=float(len(s1.intersection(s2)))
	n2=float(len(s1.union(s2)))
	try:
		e=n1/n2
		return e
	except ZeroDivisionError:
		return 0
	
	
def flatten(ls):
	return [item for sublist in ls for item in sublist]
		
def two_sort(l1,l2):
	qb=qx.Queue()
	for w in l2:
		qb.put(w)
	l=len(l1)
	l3=[]
	a=[]
	b=[]
	q=qx.Queue()
	for i in range(l):
		if type(l1[i])!=list:
			a.append(l1[i])
		else:
			a=a+flatten(l1[i])
		while len(b) < len(a):
			bn=qb.get(False)
			if type(bn)!=list:
				b.append(bn)
			else:
				b=b+flatten(bn)
		for w in a:
			if w in b:
				l3.append(w)
				q.put(w,False)
		to_dela=[]
		to_delb=[]
		while not q.empty():
			s=q.get(False)
			for i in range(len(a)):
				if a[i]==s:
					to_dela.append(i)
			for i in range(len(b)):
				if b[i]==s:
					to_delb.append(i)
		if len(to_dela)>0:
			a=del_indices(a,to_dela)
		if len(to_delb)>0:
			b=del_indices(b,to_delb)
	return l3			
	
								
def similarity(H):
	G=nx.DiGraph(H)
	for ed in G.edges_iter():
		n1=ed[0]
		n2=ed[1]
		s11=G.neighbors(n1)
		s12=G.predecessors(n1)
		s21=G.neighbors(n2)
		s22=G.predecessors(n2)
		G[n1][n2]['sim']=max([jaccard(s11,s21),jaccard(s11,s22),jaccard(s12,s21),jaccard(s12,s22)])
	return G
	
def level(H,key='y coordinate'):
	G=nx.DiGraph(H)
	h=[]
	for nod in G.nodes_iter():
		try:
			h.append(G.node[nod][key])
		except KeyError:
			pass
	L=sorted(list(set(h)))
	a=L[0]-1
	b=L[0]+1
	L=[a]+L+[b]
	for nod in G.nodes_iter():
		try:
			w=G.node[nod][key]
			G.node[nod]['level']=L.index(w)
		except KeyError:
			if G.in_degree(nod)==0:
				G.node[nod]['level']=len(L)-1
			elif G.out_degree(nod)==0:
				G.node[nod]['level']=0
			else:
				print 'Error in level program'
				return None
	return G

def make_ready(fName,suffix='_C'):
	G=nx.read_graphml('scimetric/'+fName+'.xml')
	G1=add_source(G)
	G2=add_sink(G1)
	G3=level(G2)
	print G.size()
	G4=spc(G3)
	G5=similarity(G4)
	nx.write_graphml(G5,'scimetric/'+fName+suffix+'.xml')
	
def tag_node(f1,f2,tag='date',label='Label'):
	G=nx.read_graphml('scimetric/'+f1+'.xml')
	dic={}
	with open('scimetric/'+f2) as f:
		for line in f:
			iD,date=line.split()
			dic[iD]=date
	for nod in G.nodes_iter():
		w=G.node[nod][label]
		G.node[nod][tag]=dic[w]
	nx.write_graphml(G,'scimetric/'+f1+'.xml')

def date2num(a):
	y,m,d=a.split('-')
	return int(y)*12+int(m)
	
def level_by_date(G,n):
	#sources and sinks will be assigned to level 0-n
	# n is months epoch
	date1=[]
	for w in G.nodes_iter():
		date1.append(G.node[w]['date'])
	date=sorted(list(set(date1)))
	dic={}
	l=0
	cur=date[0]
	dic[cur]=l
	z=[]
	for ymd in date:
		if (date2num(ymd)-date2num(cur)) < n:
			dic[ymd]=l
		else:
			l=l+1
			dic[ymd]=l
			cur=ymd
	for w in G.nodes_iter():
		t=G.node[w]['date']
		G.node[w]['epoch']=dic[t]
	return G
		
def phylogeny(G,imp='sig',similarity=False):
	H=nx.DiGraph()
	dic={}
	w=source(G)
	if similarity:
		if len(w)>1:
			print 'Currently similarity could be set only for single source matrices'
	l1=set(nx.neighbors(G,w[0])+G.predecessors(w[0])+w)
	for node in G.nodes_iter():
		key=G.node[node][level0]
		sig=G.node[node][imp]
		if similarity:
			l2=set(nx.neighbors(G,node)+G.predecessors(node)+[node])
			sim=jaccard(l1,l2)
		try:
			if similarity:
				dic[key].append((sig,sim,node))
			else:
				dic[key].append((sig,node))
		except KeyError:
			if similarity:
				dic[key]=[(sig,sim,node)]
			else:
				dic[key]=[(sig,node)]
	l0=min(dic.keys())
	le=max(dic.keys())
	K=le-l0+1
	print 'Detected epoches: '+str(K) 
	ls=[]
	win=[]
	for l in range(l0,le+1):
		if similarity:
			a=list(zip(*sorted(dic[l],reverse=True))[-1])
			b=list(zip(*sorted(dic[l],reverse=True,key=lambda tup:tup[1]))[-1])
			tempx=two_sort(a,b)
			ls.append(tempx)
			win.append([tempx[0]])
		else:
			tempx=list(zip(*sorted(dic[l],reverse=True))[-1])
			ls.append(tempx)
			win.append([tempx[0]])
	# dic shows the chosen numbers of
	state=dict(zip(*[range(l0,le+1),[1]*K]))
	fg=True
	while fg:
		upward=True
		downward=True
		for l in range(l0,le):
			nb=win[l]
			nu=win[l+1]
			for wb in nb:
				#~ H.add_node(wb,G.node[wb])
				for wu in nu:
					try:
						H.add_edge(wu,wb,G[wu][wb])
						H.add_node(wu,G.node[wu])
						H.add_node(wb,G.node[wb])
					except KeyError:
						if H.in_degree(wb)==0 and nu[-1]==wu:
							try:
								tt=-1
								her=None
								for xx in G.predecessors(wb):
									if G.node[xx][level0]==l+1:
										if G.node[xx][imp]>tt:
											her=xx
											tt=G.node[xx][imp]
								if her:
									upward=False
									win[l+1].append(her)
							except IndexError:
								pass
		for l in reversed(range(l0+1,le+1)):
			nb=win[l-1]
			nu=win[l]
			for wu in nu:
				#~ H.add_node(wu,G.node[wu])
				for wb in nb:
					try:
						H.add_edge(wu,wb,G[wu][wb])
						H.add_node(wu,G.node[wu])
						H.add_node(wb,G.node[wb])
					except KeyError:
						if H.out_degree(wu)==0 and nb[-1]==wb:
							try:
								tt=-1
								her=None
								for xx in G.neighbors(wu):
									if G.node[xx][level0]==l-1:
										if G.node[xx][imp]>tt:
											her=xx
											tt=G.node[xx][imp]
								if her:
									downward=False
									win[l-1].append(her)
							except IndexError:
								pass
		if upward and downward:
			break
	nod=H.nodes()
	for com in combinations(nod,2):
		n1=com[0]
		n2=com[1]
		try:
			di=G[n1][n2]
			if abs(G.node[n1][level]-G.node[n2][level])<2:
				H.add_edge(n1,n2,di)
		except KeyError:
			try:
				di=G[n2][n1]
				if abs(G.node[n1][level]-G.node[n2][level])<2:
					H.add_edge(n1,n2,di)
			except KeyError:
				pass
	return H
