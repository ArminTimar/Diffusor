"""Archive source parameters and literal formulas for independent age checks."""
from pathlib import Path
import json
import hashlib
import re
import openpyxl

ROOT=Path('papers/Supplementaries and data/gordeychik2018')
OUT=Path(__file__).resolve().parent
def read(name,sheet):
    p=ROOT/name
    w=openpyxl.load_workbook(p,data_only=True)
    f=openpyxl.load_workbook(p,data_only=False)
    return p,list(w[sheet].values),list(f[sheet].values)

p,rr,_=read('41598_2018_30133_MOESM4_ESM.xlsx','Table SM2-B')
orientations={}
orientation_names={}
orientation_apr={}
normalise=lambda s:re.sub(r'[^a-z0-9]','',s.lower())
for i,r in enumerate(rr):
    if isinstance(r[0],str) and r[0].startswith('SHIV'):
        orientations[normalise(r[0])]=[rr[i+j][9] for j in range(3)]
        orientation_names[normalise(r[0])]=r[0]
        # the 'Profile' row below a/b/c holds the geometric factor of this traverse
        orientation_apr[normalise(r[0])]=rr[i+3][10]
p,rr,ff=read('41598_2018_30133_MOESM5_ESM.xlsx','Table SM3-C')
conditions={r[1]:dict(Fo=r[2],low=dict(fo2_Pa=r[3],P_Pa=r[4],T_C=r[5],D=r[6]),
                     high=dict(fo2_Pa=r[8],P_Pa=r[9],T_C=r[10],D=r[11]))
            for r in rr[3:8]}
records=[]
for index,group in [(6,'outer_core'),(7,'advanced_core'),(8,'core_overgrowth')]:
    name=f'41598_2018_30133_MOESM{index}_ESM.xlsx'
    sheet=f'Table SM{index-2}-A'
    p,rr,ff=read(name,sheet)
    for i,r in enumerate(rr[3:],4):
        if not isinstance(r[0],(int,float)) or not isinstance(r[1],str):continue
        offset=0 if index<8 else 5
        source=conditions['Transition zone' if index<8 else 'Overgrowth']
        # pair an orientation only when its geometric factor is the one the age row uses
        name=normalise(r[1])
        same_apr=name in orientation_apr and abs(orientation_apr[name]-r[7+offset])<1e-9
        other=[orientation_names[k] for k,v in orientation_apr.items()
               if abs(v-r[7+offset])<1e-9 and k[:-1]==name[:-1]]
        item=dict(group=group,sample=r[1],part=r[2],row=i,
                  source_file=p.name,source_sheet=sheet,source_sha256=hashlib.sha256(p.read_bytes()).hexdigest(),
                  cosine_squared=orientations.get(name) if same_apr else None,
                  orientation_source_name=orientation_names.get(name) if same_apr else None,
                  orientation_note=None if same_apr or name not in orientation_apr else
                      (f"Table SM2-B gives {orientation_names[name]} the geometric factor {orientation_apr[name]:.4f}; "
                       f"this age row uses {r[7+offset]:.4f}"
                       + (f", the factor SM2-B lists for {', '.join(other)}" if other else "")
                       + ". No orientation is paired with it."),
                  Apr=r[7+offset],conditions=source,
                  published_Fo_days=[r[10+offset],r[11+offset]],
                  published_Ni_days=[r[14+offset],r[15+offset]],
                  D_Fo=[r[8+offset],r[9+offset]],D_Ni=[r[12+offset],r[13+offset]],
                  source_Fo_formulas=[ff[i-1][10+offset],ff[i-1][11+offset]],
                  source_Ni_formulas=[ff[i-1][14+offset],ff[i-1][15+offset]])
        if index==7:item['Dt_mm2']=r[6]
        else:
            item['width_Fo_mm']=r[5]
            item['width_Ni_mm']=r[6] if index==6 else r[9]
        records.append(item)
(OUT/'gordeychik_parameters.json').write_text(json.dumps(records,indent=2,ensure_ascii=False)+'\n',encoding='utf8')
print(len(records),'published model rows')
