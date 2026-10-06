"""Read supplied workbook data only; preserve inputs and stream through blank-format rows."""
from pathlib import Path
import csv, zipfile, xml.etree.ElementTree as E, os
base=Path(os.environ.get('RNASEQ_PROJECT',Path(__file__).resolve().parents[1]))
ns={'m':'http://schemas.openxmlformats.org/spreadsheetml/2006/main'}
out=base/'results';out.mkdir(exist_ok=True)
for p in (base/'inputs').glob('*.xlsx'):
    kind='BMDC' if 'BMDC' in p.name else 'PBMC'
    records=[]
    with zipfile.ZipFile(p) as z:
        ss=[]
        if 'xl/sharedStrings.xml' in z.namelist():
            ss=[''.join(e.itertext()) for e in E.fromstring(z.read('xl/sharedStrings.xml'))]
        book=E.fromstring(z.read('xl/workbook.xml'))
        rels=E.fromstring(z.read('xl/_rels/workbook.xml.rels'))
        targets={r.attrib['Id']:r.attrib['Target'] for r in rels}
        for sh in book.find('m:sheets',ns):
            sheet=sh.attrib['name']
            target=targets[sh.attrib['{http://schemas.openxmlformats.org/officeDocument/2006/relationships}id']]
            target=target.lstrip('/') if target.startswith('/') else 'xl/'+target
            headers={}
            with z.open(target) as f:
                for _,row in E.iterparse(f,events=['end']):
                    if row.tag!='{'+ns['m']+'}row': continue
                    cells={}
                    for c in row:
                        col=''.join(x for x in c.attrib.get('r','') if x.isalpha())
                        if len(col)>1: continue # Relevant data are columns A:V.
                        v=c.find('m:v',ns); t=c.find('m:is',ns)
                        val=v.text if v is not None else ''.join(t.itertext()) if t is not None else ''
                        if c.attrib.get('t')=='s': val=ss[int(val)]
                        if val: cells[col]=val
                    if cells.get('A')=='Panel': headers=cells
                    elif cells.get('J','').startswith(('GO:','hsa','mmu')):
                        rec={headers.get(col,col):val for col,val in cells.items()}
                        rec.update(Sheet=sheet,Row=row.attrib['r']);records.append(rec)
                    row.clear()
    keys=['Sheet','Row','Panel','Comparison','Direction','Database','Selected','PlotOrder','ID','Description','DisplayLabel','Theme','GeneName','geneID','Count','pvalue','p.adjust','GeneRatio','BgRatio']
    with (out/(kind+'_excel_rows.tsv')).open('w') as f:
        w=csv.DictWriter(f,fieldnames=keys,delimiter='\t',extrasaction='ignore');w.writeheader();w.writerows(records)
    chosen=[r for r in records if r['Sheet']==('Sheet4' if kind=='BMDC' else 'Sheet5')]
    ids=list(dict.fromkeys(r['ID'] for r in chosen))
    if len(ids)!=10: raise ValueError(f'{kind}: expected 10 unique selected pathway IDs, found {len(ids)}')
    with (out/(kind+'_selected_pathways.tsv')).open('w') as f:
        w=csv.writer(f,delimiter='\t');w.writerow(['ID','Description','Direction','SourceSheet'])
        for id in ids:
            rr=next(r for r in chosen if r['ID']==id)
            w.writerow([id,rr['Description'],'Up',rr['Sheet']])
    print(kind, 'selected',len(ids),'unique pathways')
