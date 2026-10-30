"""Build matched upstream atom-equivalence helpers differing only in sorting.

The installed Antechamber binary is never modified. Full path enumeration and
exact score comparisons remain as in the pinned upstream source.
"""
import hashlib
import json
from pathlib import Path
import shutil
import subprocess

SOURCE=Path('/workspace/md-tools/antechamber-source-review')
OUT=Path('/workspace/md-tools/equatom-qsort-v1')

DRIVER=r'''
#include <assert.h>
#include <ctype.h>
#include <errno.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "define.h"
#include "atom.h"
#include "eprintf.h"
#include "equatom.c"

int main(int argc,char **argv) {
    if(argc != 4) return 2;
    FILE *in=efopen(argv[1],"r");int n;
    if(fscanf(in,"%d",&n)!=1 || n<1 || n>MAX_EQU_ATOM) return 3;
    ATOM *atoms=ecalloc(n,sizeof(ATOM));int *equ=emalloc(n*sizeof(int));
    for(int i=0;i<n;i++) {
        if(fscanf(in,"%d %9s",&atoms[i].atomicnum,atoms[i].name)!=2) return 4;
        for(int j=0;j<6;j++) {
            if(fscanf(in,"%d",&atoms[i].con[j])!=1) return 5;
            if(atoms[i].con[j]>=n || atoms[i].con[j]<-1 || atoms[i].con[j]==i) return 6;
        }
    }
    fclose(in);
    identify_equatom(n,atoms,equ,atoi(argv[3]),0,NULL,0);
    FILE *out=efopen(argv[2],"w");
    fprintf(out,"{\"representatives_1based\":[");
    for(int i=0;i<n;i++) fprintf(out,"%s%d",i?",":"",(equ[i]<0?i:equ[i])+1);
    fprintf(out,"],\"path_counts\":[");
    for(int i=0;i<n;i++) fprintf(out,"%s%d",i?",":"",pathnum[i]);
    fprintf(out,"]}\n");fclose(out);return 0;
}
'''

NEW_SORT=r'''static int compare_path_score(const void *a, const void *b)
{
    const double x=*(const double *)a, y=*(const double *)b;
    return (x > y) - (x < y);
}

void sort(double array[], int elemnum)
{
    for(int i=0;i<elemnum;i++)
        if(!isfinite(array[i])) eprintf("Nonfinite path score; no equivalence output.");
    qsort(array, (size_t)elemnum, sizeof(double), compare_path_score);
}

'''


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    OUT.mkdir()
    sources={}
    for name in ['define.h','atom.h','eprintf.h','eprintf.c','rotate.c','equatom.c']:
        p=SOURCE/name;metadata=json.loads((SOURCE/(name+'.provenance.json')).read_text())
        if sha(p)!=metadata['sha256']:
            raise ValueError('Pinned source checksum mismatch')
        sources[name]=metadata
    original=(SOURCE/'equatom.c').read_text()
    begin=original.index('void sort(double array[], int elemnum)')
    end=original.index('void equatom(',begin)
    modified=original[:begin]+NEW_SORT+original[end:]
    binaries={}
    for label,text in [('reference',original),('qsort',modified)]:
        folder=OUT/label;folder.mkdir()
        for name in sources:
            shutil.copy2(SOURCE/name,folder/name)
        (folder/'equatom.c').write_text(text)
        (folder/'driver.c').write_text(DRIVER)
        command=['gcc','-O2','-fno-fast-math','-ffp-contract=off','driver.c','eprintf.c','-lm','-o','equatom-helper']
        result=subprocess.run(command,cwd=folder,capture_output=True,text=True)
        (folder/'build.log').write_text(result.stdout+result.stderr)
        result.check_returncode()
        binaries[label]={'path':str(folder/'equatom-helper'),'sha256':sha(folder/'equatom-helper'),
                         'equatom_source_sha256':sha(folder/'equatom.c')}
    report={'pinned_sources':sources,'binaries':binaries,
            'change':'Only sort implementation differs between helpers;full path enumeration,score arithmetic,equality and representative selection unchanged.',
            'installed_antechamber_replaced':False,'validation_pending':True}
    (OUT/'build.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(binaries))


if __name__=='__main__':
    main()
