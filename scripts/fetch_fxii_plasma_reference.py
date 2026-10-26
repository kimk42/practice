"""Record a sourced physiological FXII reference, separate from user-measured assay data."""
import json,hashlib,urllib.request,xml.etree.ElementTree as ET
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
def main():
    pmc='PMC9742269';url=f'https://www.ebi.ac.uk/europepmc/webservices/rest/{pmc}/fullTextXML';path=ROOT/'research/proteins/fxii/raw'/(pmc+'.xml')
    if not path.exists():
        with urllib.request.urlopen(url,timeout=30) as response:path.write_bytes(response.read())
    data=path.read_bytes();root=ET.fromstring(data)
    evidence=[''.join(p.itertext()) for p in root.findall('.//p') if 'physiological concentration (375 nM)' in ''.join(p.itertext()) and 'human plasma-derived FXII' in ''.join(p.itertext())]
    if not evidence:raise ValueError('Expected physiological-reference evidence absent')
    result={'source_title':''.join(root.find('.//article-title').itertext()),'doi':'10.3389/fcvm.2022.1008410','pmcid':pmc,'source_url':url,'source_sha256':hashlib.sha256(data).hexdigest(),'short_evidence':'human plasma-derived FXII at physiological concentration (375 nM)','context':'human FXII reconstitution of F12-deficient mouse whole blood;literature physiological reference,not clinical reference-interval determination or user assay measurement','physiological_reference_nM':375,'user_target_fraction':0.01,'provisional_1percent_target_nM':3.75,'mass_conversion_if_80kDa_ug_mL':0.3,'80kDa_status':'conditional glycoprotein-mass convention,not a mass measured by this source/user or the unglycosylated MD model','user_assay_measured_concentration':False}
    out=ROOT/'research/proteins/modeling-evidence/fxii-plasma-reference.json';out.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result))
if __name__=='__main__':main()
