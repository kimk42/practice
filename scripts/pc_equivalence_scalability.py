"""Prove a lower bound on unrestricted PC atom-path enumeration without running it.

Disjoint simple rings connected as a tree give two vertex-simple routes through
each ring on an end-to-end path. This is a resource gate, not a charge method.
"""
from datetime import datetime, timezone
import json
from pathlib import Path

from rdkit import Chem

from pc_length_parameters import ROOT, preflight, sha
from run_pc_sqm_checkpointed import write_json


def path_lower_bound(mol):
    rings=[set(ring) for ring in mol.GetRingInfo().AtomRings()]
    membership={}
    for i,ring in enumerate(rings):
        if len(ring)!=6 or any(a in membership for a in ring):
            raise ValueError('Expected disjoint six-membered rings')
        for atom in ring:
            neighbors=[n.GetIdx() for n in mol.GetAtomWithIdx(atom).GetNeighbors() if n.GetIdx() in ring]
            if len(neighbors)!=2:
                raise ValueError('Ring is not a simple cycle')
            membership[atom]=('ring',i)
    nodes={a.GetIdx():membership.get(a.GetIdx(),('atom',a.GetIdx())) for a in mol.GetAtoms()}
    graph={node:set() for node in nodes.values()}
    for bond in mol.GetBonds():
        a=nodes[bond.GetBeginAtomIdx()];b=nodes[bond.GetEndAtomIdx()]
        if a==b:continue
        if b in graph[a]:raise ValueError('Parallel inter-block edges invalidate the tree proof')
        graph[a].add(b);graph[b].add(a)
    ends=[a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum()==8 and a.GetDegree()==2
          and sorted(n.GetAtomicNum() for n in a.GetNeighbors())==[1,6]
          and any(n.GetIsAromatic() for n in a.GetNeighbors())]
    if len(ends)!=2:raise ValueError('Two phenolic OH ends required')
    first,last=[nodes[i] for i in ends];parents={first:None};queue=[first]
    for node in queue:
        for neighbor in graph[node]:
            if neighbor not in parents:parents[neighbor]=node;queue.append(neighbor)
    if len(parents)!=len(graph) or sum(map(len,graph.values()))!=2*(len(graph)-1):
        raise ValueError('Contracted graph is not a connected tree')
    path=[];node=last
    while node is not None:path.append(node);node=parents[node]
    count=sum(node[0]=='ring' for node in path)
    bound=2**count
    return {'phenolic_O_end_indices_1based':[i+1 for i in ends],
            'disjoint_simple_rings':len(rings),'rings_on_end_to_end_path':count,
            'end_to_end_simple_paths_lower_bound':bound,
            'score_array_bytes_lower_bound_for_one_atom':8*bound,
            'exceeds_signed_int32_path_counter':bound>2**31-1,
            'proof':'contracted graph is a tree;each simple disjoint ring on the end-to-end route admits two distinct internal paths;their independent choices multiply'}


def require_supported_full_path_count(dp):
    mol,_=preflight(dp);bound=path_lower_bound(mol)
    if bound['exceeds_signed_int32_path_counter']:
        raise ValueError(f'PC{dp} unrestricted equivalence held: at least {bound["end_to_end_simple_paths_lower_bound"]} paths from one end atom exceed the pinned upstream signed-int counter;review an alternative before postprocessing')
    return bound


def main():
    report={'checked_utc':datetime.now(timezone.utc).isoformat(),'cases':{},
            'method_scope':'graph lower bound only;no charge assignment,no path truncation,no QM interruption',
            'installed_build_caveat':'pinned upstream identifies int counters and retained double arrays;exact installed binary source is not established',
            'decision':'Hold automatic unrestricted PC16 equivalence;keep existing QM running. PC4/PC8 existing validated candidates unchanged.'}
    source=Path('/workspace/md-tools/antechamber-source-review/equatom.c')
    text=source.read_text()
    if 'int pathnumindex = 0;' not in text or 'double *pathscore[MAX_EQU_ATOM];' not in text:
        raise ValueError('Pinned upstream counter/storage evidence changed')
    report['pinned_equatom_source_sha256']=sha(source)
    for dp in (4,8,16):
        mol,info=preflight(dp);result=path_lower_bound(mol)
        result['molecule_sha256']=info['sha256']
        if dp in (4,8):
            reference=json.loads((ROOT/f'runs/equatom-qsort-validation-v1/pc{dp}/qsort-path-1.json').read_text())
            actual=reference['path_counts'][result['phenolic_O_end_indices_1based'][0]-1]
            if actual<result['end_to_end_simple_paths_lower_bound']:
                raise ValueError('Graph lower bound contradicts complete enumeration')
            result['actual_full_path_count_for_first_end_O']=actual
        report['cases'][str(dp)]=result
    write_json(ROOT/'models/surface-revision-v2/pc-path-scalability-audit.json',report)
    print(json.dumps(report,indent=2))


if __name__=='__main__':main()
