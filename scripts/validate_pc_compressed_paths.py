"""Compare compressed full-path scores with explicit upstream enumeration."""
from datetime import datetime, timezone
import json
from pathlib import Path
import resource
import shutil
import subprocess

from rdkit import Chem

from build_equatom_qsort import DRIVER
from pc_compressed_paths import score_histograms,representatives
from pc_length_parameters import ROOT,preflight,sha
from run_pc_sqm_checkpointed import write_json

OUT=ROOT/'runs/pc-compressed-path-validation-v1'
BUILD=Path('/workspace/md-tools/equatom-qsort-v1')


def limits():
    resource.setrlimit(resource.RLIMIT_AS,(3*1024**3,3*1024**3))
    resource.setrlimit(resource.RLIMIT_CPU,(300,300))


def main():
    OUT.mkdir();build=OUT/'helper';build.mkdir()
    for file in (BUILD/'qsort').iterdir():
        if file.suffix in ('.c','.h'):shutil.copy2(file,build/file.name)
    original='fprintf(out,"]}\\n");fclose(out);return 0;'
    extra=r'''fprintf(out,"],\"histograms\":[");
    for(int i=0;i<n;i++) {
        fprintf(out,"%s[",i?",":"");int first=1;
        for(int j=0;j<pathnum[i];) {
            int end=j+1;while(end<pathnum[i] && pathscore[i][end]==pathscore[i][j])end++;
            fprintf(out,"%s[\"%a\",%d]",first?"":",",pathscore[i][j],end-j);first=0;j=end;
        }
        fprintf(out,"]");
    }
    fprintf(out,"]}\n");fclose(out);return 0;'''
    if original not in DRIVER:raise ValueError('Helper driver marker changed')
    (build/'driver.c').write_text(DRIVER.replace(original,extra))
    subprocess.run(['gcc','-O2','-fno-fast-math','-ffp-contract=off','driver.c','eprintf.c','-lm','-o','histogram-helper'],cwd=build,check=True)
    results={}
    for dp in (4,8,16):
        mol,info=preflight(dp);compressed=score_histograms(mol);ids=representatives(compressed)
        result={'molecule':info,'unrestricted_path_counts':[sum(h.values()) for h in compressed],
                'distinct_score_counts':[len(h) for h in compressed],
                'representatives_1based':ids,'equivalence_classes':len(set(ids))}
        ranks=list(Chem.CanonicalRankAtoms(mol,breakTies=False,includeChirality=False))
        result['canonical_graph_partition_identical']=ids==[ranks.index(r)+1 for r in ranks]
        if dp in (4,8):
            graph=ROOT/f'runs/equatom-qsort-validation-v1/pc{dp}/graph.txt'
            dest=OUT/f'pc{dp}-enumerated.json'
            with (OUT/f'pc{dp}-enumeration.log').open('x') as log:
                subprocess.run([str(build/'histogram-helper'),str(graph),str(dest),'-1'],
                               stdout=log,stderr=subprocess.STDOUT,check=True,timeout=330,preexec_fn=limits)
            enumerated=json.loads(dest.read_text())
            for i,(expected,actual) in enumerate(zip(enumerated['histograms'],compressed)):
                full={float.fromhex(score):count for score,count in expected}
                if full!=actual:raise ValueError(f'PC{dp} atom{i+1} score/multiplicity mismatch')
            if enumerated['path_counts']!=result['unrestricted_path_counts'] or enumerated['representatives_1based']!=ids:
                raise ValueError('Explicit path count or equivalence differs')
            result['every_exact_binary64_score_and_integer_multiplicity_matches_enumeration']=True
            result['enumeration_sha256']=sha(dest)
        else:
            result['explicit_enumeration_performed']=False
        write_json(OUT/f'pc{dp}-compressed.json',[[[score.hex(),count] for score,count in sorted(h.items())] for h in compressed])
        result['compressed_histograms_sha256']=sha(OUT/f'pc{dp}-compressed.json')
        results[str(dp)]=result
        print(json.dumps({'dp':dp,'classes':len(set(ids)),'max_full_paths':max(result['unrestricted_path_counts']),
                          'max_distinct_scores':max(result['distinct_score_counts']),
                          'explicit_comparison_passed':dp in (4,8)}),flush=True)
    report={'checked_utc':datetime.now(timezone.utc).isoformat(),'cases':results,
            'source_algorithm':'pinned full-path igeom=0 score;no path limit,no float rounding',
            'helper_sha256':sha(build/'histogram-helper'),'compressed_algorithm_sha256':sha(ROOT/'scripts/pc_compressed_paths.py'),
            'physical_approval':False,'charges_assigned':False,
            'limitation':'explicit full enumeration comparison covers PC4 and PC8;PC16 relies on the validated tree-of-disjoint-rings decomposition,not a direct billions-of-paths execution'}
    write_json(OUT/'result.json',report)
    write_json(ROOT/'models/surface-revision-v2/pc-compressed-path-validation.json',report)


if __name__=='__main__':main()
