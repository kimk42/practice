import test from 'node:test';
import assert from 'node:assert/strict';
import {parsePDB,contacts} from '../model.mjs';
const row=(serial,name,x,y,z)=>`ATOM  ${String(serial).padStart(5)}  ${name.padEnd(3)} ALA A   1    ${x.toFixed(3).padStart(8)}${y.toFixed(3).padStart(8)}${z.toFixed(3).padStart(8)}  1.00 20.00           C`;
test('PDB selects C alpha atoms from first model and centers coordinates',()=>{
 const data=[row(1,'CA',10,20,30),row(2,'N',0,0,0),row(3,'CA',20,20,30),'ENDMDL',row(4,'CA',99,99,99)].join('\n');
 const atoms=parsePDB(data);assert.equal(atoms.length,2);assert.equal(atoms[0].x,-5);assert.equal(atoms[1].x,5);assert.equal(atoms[0].chain,'A');
});
test('invalid input fails explicitly',()=>assert.throws(()=>parsePDB('invalid')));
test('plane placement maintains gap and cutoff identifies expected residue',()=>{
 const atoms=[{x:0,y:0,z:0},{x:0,y:0,z:10}];
 const r=contacts(atoms,0,5,8);assert.equal(r.floor,-5);assert.deepEqual(r.atoms.map(a=>a.contact),[true,false]);
 const turned=contacts(atoms,Math.PI,5,8);assert.deepEqual(turned.atoms.map(a=>a.contact),[false,true]);assert.equal(Math.min(...turned.atoms.map(a=>a.distance)),5);
});
