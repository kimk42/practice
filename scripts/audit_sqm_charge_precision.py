"""Read-only tracing of rounded SQM charge output;never manufacture recovered digits."""
import json,re,hashlib
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
def main():
    base=ROOT/json.loads((ROOT/'runs/buffer-charge-active.json').read_text())['directory'];reports=[]
    species=json.loads((ROOT/'models/buffer-parameters/polyatomic-species-audit.json').read_text())['candidates']
    for item in species:
        if item['name']=='azide-':continue
        folder=base/item['name'];path=folder/'sqm.out';lines=path.read_text().splitlines();blocks=[]
        for index,line in enumerate(lines):
            if 'Element' not in line or 'Mulliken Charge' not in line:continue
            values=[];digits=[];total=None
            for row in lines[index+1:]:
                match=re.match(r'^\s*\d+\s+[A-Za-z]+\s+(-?\d+\.\d+)\s*$',row)
                if match:values.append(float(match.group(1)));digits.append(len(match.group(1).split('.')[1]))
                match_total=re.search(r'Total Mulliken Charge\s*=\s*(-?\d+\.\d+)',row)
                if match_total:total=float(match_total.group(1));break
            if total is not None and len(values)==item['atoms']:blocks.append((values,digits,total))
        if not blocks:raise ValueError('Expected complete charge table missing:'+item['name'])
        values,digits,total=blocks[-1];raw_sum=sum(values);bound=sum(.5*10**(-d) for d in digits)
        pre=folder/'ANTECHAMBER_AM1BCC_PRE.AC';pre_values=[float(x.split()[-2]) for x in pre.read_text().splitlines() if x.startswith('ATOM')]
        mol2=(folder/'candidate.mol2').read_text().split('@<TRIPOS>ATOM')[1].split('@<TRIPOS>')[0];bcc=sum(float(x.split()[8]) for x in mol2.splitlines() if x.strip())
        row={'species':item['name'],'expected_integer_charge':item['formal_charge'],'QM_reported_total':total,'printed_atomic_decimal_places':sorted(set(digits)),'sum_printed_atomic_charges':raw_sum,'sum_BCC_input_after_equivalence_handling':sum(pre_values),'sum_BCC_output':bcc,'maximum_sum_rounding_uncertainty_e':bound,'printed_total_residual_e':raw_sum-total,'pre_BCC_sum_matches_printed_sum':abs(sum(pre_values)-raw_sum)<1e-5,'BCC_sum_matches_pre_BCC_sum':abs(bcc-sum(pre_values))<1e-5,'compatible_with_decimal_rounding':abs(raw_sum-total)<=bound+1e-12,'source_sqm_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'charges_changed':False}
        reports.append(row)
    result={'state':'printed-charge precision chain inspected without new calculations','species':reports,'interpretation':'phosphate net-charge residual is consistent with3-decimal SQM atomic-charge export propagated into charge-conserving BCC;not an independently verified high-precision electronic charge set','action':'preserve originals;recover higher-precision supported output or use a documented independently validated charge model before accepting phosphate topology','manual_renormalization_applied':False,'MD_topology_exported':False,'full_buffer_ready':False}
    (ROOT/'models/buffer-parameters/charge-precision-audit.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps([{'species':x['species'],'printed_sum':x['sum_printed_atomic_charges'],'QM_total':x['QM_reported_total'],'BCC_sum':x['sum_BCC_output'],'rounding_compatible':x['compatible_with_decimal_rounding']} for x in reports]))
if __name__=='__main__':main()
