"""Audit explicit fibrinogen sequence candidates, preserving fibrinopeptides.
This supplies no missing coordinates and does not select a verified sample isoform.
"""
import json,re,hashlib
from pathlib import Path
import gemmi
ROOT=Path(__file__).resolve().parents[1]
def main():
    raw=ROOT/'research/proteins/fibrinogen/raw';out=ROOT/'models/proteins'
    block=gemmi.cif.read_file(str(raw/'3GHG-assembly1.cif')).sole_block();poly=block.get_mmcif_category('_entity_poly.')
    rows=[];fasta=''
    for accession,base,entity,name in [('P02671-2','P02671','1','Aalpha'),('P02675','P02675','2','Bbeta'),('P02679-2','P02679','3','gammaA')]:
        data=json.loads((raw/(base+'-uniprot.json')).read_text());signal=next(f for f in data['features'] if f['type']=='SIGNAL');start=int(signal['end'])+1
        file=raw/(accession+'.fasta')
        if accession!=base and not file.exists():raise ValueError('Required isoform FASTA absent; canonical sequence is not a substitute: '+accession)
        if file.exists() and accession not in file.read_text().splitlines()[0]:raise ValueError('FASTA accession mismatch')
        sequence=''.join(file.read_text().splitlines()[1:]) if file.exists() else data['sequence']['sequence']
        mature=sequence[start-1:];index=poly['entity_id'].index(entity);construct=''.join(poly['pdbx_seq_one_letter_code_can'][index].split())
        scoring=gemmi.AlignmentScoring();scoring.match=2;scoring.mismatch=-3;scoring.gapo=-8;scoring.gape=-1
        result=gemmi.align_string_sequences(list(mature),list(construct),[],scoring)
        rows.append({'candidate_chain':name,'accession':accession,'precursor_residues':len(sequence),'signal_removed_only':[1,start-1],'mature_candidate_residues':len(mature),'3GHG_entity':entity,'construct_residues':len(construct),'alignment_cigar':result.cigar_str(),'alignment_matches':result.match_count,'candidate_sequence_sha256':hashlib.sha256(mature.encode()).hexdigest(),'caution':'isoform candidate, not sample-verified; sequence alignment adds no missing coordinates; fibrinopeptides retained because target is fibrinogen, not fibrin'})
        fasta+=f'>{name}|{accession}|signal_removed_only|provisional_isoform_candidate\n'+ '\n'.join(mature[i:i+80] for i in range(0,len(mature),80))+'\n'
    report={'chains':rows,'six_chain_candidate_residue_count':2*sum(r['mature_candidate_residues'] for r in rows),'observed_3GHG_residues':1913,'gap_note':'955 unmodeled residues is relative to the3GHG declared construct, not necessarily this full-length isoform candidate','structure_model_built':False,'scope':'provisional sequence references only, no reconstructed full-length structure'}
    (out/'fibrinogen-isoform-candidates.fasta').write_text(fasta);(out/'fibrinogen-isoform-audit.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report))
if __name__=='__main__':main()
