export function parsePDB(text) {
  const atoms = [];
  for (const line of text.split(/\r?\n/)) {
    if (line.startsWith('ENDMDL')) break;
    if (!line.startsWith('ATOM  ') || line.slice(12,16).trim() !== 'CA') continue;
    if (![' ', 'A'].includes(line[16])) continue;
    const xyz = [30,38,46].map(i => Number(line.slice(i,i+8).trim() || NaN));
    if (!xyz.every(Number.isFinite)) continue;
    atoms.push({x:xyz[0],y:xyz[1],z:xyz[2],chain:line[21],residue:line.slice(17,20).trim(),id:line.slice(22,27).trim()});
  }
  if (!atoms.length) throw new Error('첫 모델에서 유효한 Cα 원자를 찾지 못했습니다. 단백질 PDB 파일을 선택하세요.');
  const center = ['x','y','z'].map(k => atoms.reduce((s,a)=>s+a[k],0)/atoms.length);
  return atoms.map(a=>({...a,x:a.x-center[0],y:a.y-center[1],z:a.z-center[2]}));
}
export function rotate(a, angle) {
  return {...a,x:a.x*Math.cos(angle)-a.z*Math.sin(angle),z:a.x*Math.sin(angle)+a.z*Math.cos(angle)};
}
export function contacts(atoms, angle, gap, cutoff) {
  const turned=atoms.map(a=>rotate(a,angle));
  const floor=Math.min(...turned.map(a=>a.z))-gap;
  return {floor,atoms:turned.map(a=>({...a,distance:a.z-floor,contact:a.z-floor<=cutoff}))};
}
