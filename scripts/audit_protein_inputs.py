"""Audit protein scope and prepare a labeled mature FXII coordinate candidate.
No missing loops, glycans or experimental conformations are invented here.
"""
from collections import defaultdict
import hashlib
import json
from pathlib import Path
import numpy as np
import gemmi

ROOT=Path(__file__).resolve().parents[1]
PROTEINS=ROOT/'research/proteins'
OUT=ROOT/'models/proteins';OUT.mkdir(parents=True,exist_ok=True)
AA={'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLN':'Q','GLU':'E','GLY':'G','HIS':'H','ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S','THR':'T','TRP':'W','TYR':'Y','VAL':'V'}

def intervals(numbers):
    result=[]
    for x in sorted(set(numbers)):
        if result and x==result[-1][1]+1:result[-1][1]=x
        else:result.append([x,x])
    return result

def fxii():
    raw=PROTEINS/'fxii/raw';meta=json.loads((raw/'alphafold-api.json').read_text())[0]
    features=json.loads((raw/'uniprot-features.json').read_text())['features']
    signal=next(f for f in features if f['type']=='SIGNAL')
    start=int(signal['end'])+1
    source=raw/'AF-P00748-F1-model_v6.pdb'
    structure=gemmi.read_structure(str(source));chain=structure[0]['A']
    observed=''.join(AA[r.name] for r in chain if r.name in AA)
    if observed!=meta['sequence']:raise ValueError('FXII sequence mismatch')
    if hashlib.md5(observed.encode()).hexdigest()!=meta['sequenceChecksum']:raise ValueError('FXII official sequence checksum mismatch')
    residues={r.seqid.num:r for r in chain if r.name in AA}
    mature=[n for n in residues if n>=start]
    confidence={n:next(a.b_iso for a in residues[n] if a.name=='CA') for n in mature}
    disulfides=[]
    for f in features:
        if f['type']!='DISULFID':continue
        a,b=int(f['begin']),int(f['end'])
        first=next((x for x in residues[a] if x.name=='SG'),None);second=next((x for x in residues[b] if x.name=='SG'),None)
        distance=first.pos.dist(second.pos) if first and second else None
        disulfides.append({'residues':[a,b],'SG_distance_A':distance,'geometry_close':distance is not None and 1.7<=distance<=2.5,'annotation':'UniProt feature; chemistry must be retained during topology preparation'})
    coordinates=np.array([[a.pos.x,a.pos.y,a.pos.z] for n in mature for a in residues[n]])
    output=OUT/'fxii-mature-unglycosylated-af-candidate.pdb'
    lines=source.read_text().splitlines()
    selected=[x for x in lines if x.startswith('ATOM  ') and int(x[22:26])>=start]
    output.write_text('REMARK   1 FULL MATURE FXII CANDIDATE FROM ALPHAFOLD; NOT EXPERIMENTAL\nREMARK   1 SIGNAL PEPTIDE REMOVED BY UNIPROT ANNOTATION; GLYCANS ABSENT\n'+'\n'.join(selected)+'\nTER\nEND\n')
    report={'source':'AF-P00748-F1-model_v6, AlphaFoldDB via official API','method':'predicted structure; not experimentally validated conformation',
            'source_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),'official_sequence_checksum_verified':True,
            'precursor_residues':len(observed),'signal_peptide':[1,start-1],'mature_range':[start,max(mature)],'mature_residue_count':len(mature),
            'zymogen_state':'single chain retained; no activation-site cleavage introduced',
            'heavy_atom_count':len(selected),'coordinate_axis_extent_nm':(np.ptp(coordinates,axis=0)/10).tolist(),
            'mean_mature_pLDDT':float(np.mean(list(confidence.values()))),
            'mature_pLDDT_below50_ranges':intervals([n for n,v in confidence.items() if v<50]),
            'mature_pLDDT_below70_fraction':sum(v<70 for v in confidence.values())/len(mature),
            'disulfide_annotations':disulfides,
            'glycosylation_annotations':[{'position':f['begin'],'description':f['description']} for f in features if f['type']=='CARBOHYD'],
            'glycans_present_in_candidate':False,'protein_force_field_assigned':False,'ready_for_production_MD':False,
            'remaining_work':['glycoform choice/modeling','low-confidence/PAE/domain-orientation assessment','disulfide topology and geometry checks','protonation/termini assignment','protein-only equilibration'],
            'candidate_file':str(output.relative_to(ROOT))}
    pae_source=raw/Path(meta['paeDocUrl']).name
    pae=np.array(json.loads(pae_source.read_text())[0]['predicted_aligned_error'])
    if pae.shape!=(len(observed),len(observed)):raise ValueError('PAE dimensions mismatch')
    matrix=pae[start-1:,start-1:];off=matrix[~np.eye(len(matrix),dtype=bool)]
    pae_report={'source_file':str(pae_source.relative_to(ROOT)),
                'definition':'AlphaFold predicted aligned error in angstroms, not measured structural error',
                'mature_offdiagonal_median_A':float(np.median(off)),
                'mature_offdiagonal_fraction_above10A':float((off>10).mean()),
                'interpretation':'Assess relative domain orientation and flexible regions before adopting a single rigid full-length starting conformation; this scalar summary does not validate domain poses.'}
    (OUT/'fxii-pae-audit.json').write_text(json.dumps(pae_report,indent=2)+'\n')
    (OUT/'fxii-audit.json').write_text(json.dumps(report,indent=2)+'\n');return report

def fibrinogen():
    raw=PROTEINS/'fibrinogen/raw';block=gemmi.cif.read_file(str(raw/'3GHG-assembly1.cif')).sole_block();structure=gemmi.make_structure_from_block(block)
    poly=block.get_mmcif_category('_entity_poly.')
    stats=[];observed_total=0;expected_total=0;xyz=[]
    for chain in structure[0]:
        residues=[r for r in chain if r.name in AA and any(a.name=='CA' for a in r)]
        if chain.name not in ['A','B','C','D','E','F'] or not residues:continue
        entity=next(i for i,s in enumerate(poly['pdbx_strand_id']) if chain.name in s.split(','))
        declared=len(poly['pdbx_seq_one_letter_code_can'][entity].replace('\n',''))
        labels=[r.label_seq for r in residues if r.label_seq is not None]
        missing=sorted(set(range(1,declared+1))-set(labels))
        expected_total+=declared;observed_total+=len(residues)
        xyz.extend([[a.pos.x,a.pos.y,a.pos.z] for r in residues for a in r])
        stats.append({'chain':chain.name,'entity_id':poly['entity_id'][entity],'declared_construct_residues':declared,'observed_residues':len(residues),'unmodeled_label_seq_ranges':intervals(missing)})
    assert len(stats)==6
    report={'source':'3GHG biological assembly1, RCSB','source_sha256':hashlib.sha256((raw/'3GHG-assembly1.cif').read_bytes()).hexdigest(),
            'scope':'experimental construct with extensive unmodeled regions; not full native fibrinogen',
            'protein_chains':stats,'declared_construct_residues_total':expected_total,'observed_protein_residues_total':observed_total,
            'unmodeled_construct_residues_total':expected_total-observed_total,
            'observed_core_axis_extent_nm':(np.ptp(np.array(xyz),axis=0)/10).tolist(),
            'other_components':'bound peptides, glycans/ligands and ions remain in raw assembly; not silently removed for a simulation',
            'ready_for_full_length_MD':False,'remaining_work':['select full-length/disordered-region modeling strategy or explicitly scoped core construct','retain or justify glycan/ligand choices','handle internal gaps and disulfides without inventing missing loops']}
    (OUT/'fibrinogen-audit.json').write_text(json.dumps(report,indent=2)+'\n');return report

if __name__=='__main__':
    a=fxii();b=fibrinogen()
    print(json.dumps({'FXII':{'mature_residues':a['mature_residue_count'],'mean_pLDDT':a['mean_mature_pLDDT'],'below70_fraction':a['mature_pLDDT_below70_fraction'],'disulfides':len(a['disulfide_annotations']),'disulfides_geometrically_close':sum(x['geometry_close'] for x in a['disulfide_annotations'])},
                      'fibrinogen':{'observed_residues':b['observed_protein_residues_total'],'unmodeled_construct_residues':b['unmodeled_construct_residues_total'],'extent_nm':b['observed_core_axis_extent_nm']}},indent=2))
