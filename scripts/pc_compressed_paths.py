"""Exact path-score multiplicities for disjoint-ring/tree molecular graphs.

The unrestricted upstream score with igeom=0 is retained, including the order
of floating-point additions along each path. Identical (depth,score) states
are combined using arbitrary-precision integer multiplicities. No paths are
discarded and no charges are assigned by this module.
"""
from collections import Counter


def blocks_and_incidence(mol):
    if any(a.GetDegree()>=6 for a in mol.GetAtoms()):
        raise ValueError('Upstream endpoint-recording convention requires degree below6')
    blocks=[];covered=set();ring_atoms=set()
    for ring in mol.GetRingInfo().AtomRings():
        vertices=set(ring)
        if len(vertices)!=6 or ring_atoms & vertices:
            raise ValueError('Only disjoint simple six-membered rings are reviewed')
        ring_atoms.update(vertices)
        edges={tuple(sorted((b.GetBeginAtomIdx(),b.GetEndAtomIdx()))) for b in mol.GetBonds()
               if b.GetBeginAtomIdx() in vertices and b.GetEndAtomIdx() in vertices}
        if len(edges)!=6 or any(sum(v in e for e in edges)!=2 for v in vertices):
            raise ValueError('Ring has chords or is not a simple cycle')
        covered.update(edges);blocks.append(edges)
    for bond in mol.GetBonds():
        edge=tuple(sorted((bond.GetBeginAtomIdx(),bond.GetEndAtomIdx())))
        if edge not in covered:blocks.append({edge})
    incidence=[[] for _ in mol.GetAtoms()];adjacency=[]
    # The bipartite block/vertex incidence graph must be a connected tree.
    tree={('v',i):set() for i in range(mol.GetNumAtoms())}
    for i,edges in enumerate(blocks):
        local={v:[] for edge in edges for v in edge}
        for a,b in sorted(edges):local[a].append(b);local[b].append(a)
        adjacency.append(local);tree[('b',i)]=set()
        for v in local:
            incidence[v].append(i);tree[('v',v)].add(('b',i));tree[('b',i)].add(('v',v))
    seen={('v',0)};queue=[('v',0)]
    for node in queue:
        for neighbor in tree[node]:
            if neighbor not in seen:seen.add(neighbor);queue.append(neighbor)
    if len(seen)!=len(tree) or sum(map(len,tree.values()))!=2*(len(tree)-1):
        raise ValueError('Block incidence is not a connected tree')
    return adjacency,incidence


def score_histograms(mol):
    adjacency,incidence=blocks_and_incidence(mol)
    elements=[a.GetAtomicNum() for a in mol.GetAtoms()]
    local_paths={}
    for index,graph in enumerate(adjacency):
        for start in graph:
            paths=[]
            def visit(path):
                for neighbor in graph[path[-1]]:
                    if neighbor not in path:
                        extended=path+(neighbor,);paths.append(extended[1:]);visit(extended)
            visit((start,));local_paths[index,start]=paths
    result=[]
    for root in range(mol.GetNumAtoms()):
        histogram=Counter()
        def walk(vertex,blocked,states):
            for (_,score),multiplicity in states.items():histogram[score]+=multiplicity
            for block in incidence[vertex]:
                if block==blocked:continue
                exits={}
                for path in local_paths[block,vertex]:
                    target=exits.setdefault(path[-1],Counter())
                    for (depth,score),multiplicity in states.items():
                        for offset,atom in enumerate(path,1):
                            # Same binary64 operations and left-to-right sum as
                            # pinned scorepath with coefficient exactly1.0.
                            score+=((depth+offset)*0.11+elements[atom]*0.08)
                        target[depth+len(path),score]+=multiplicity
                for endpoint,values in exits.items():walk(endpoint,block,values)
        walk(root,None,Counter({(1,0.11+elements[root]*0.08):1}))
        result.append(histogram)
    return result


def representatives(histograms):
    seen={};ids=[]
    for i,histogram in enumerate(histograms,1):
        signature=tuple(sorted(histogram.items()))
        if signature not in seen:seen[signature]=i
        ids.append(seen[signature])
    return ids
