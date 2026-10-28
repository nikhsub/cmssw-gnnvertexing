#!/usr/bin/env python3
"""Write prompt-D matched-track and pair-fit sanity diagnostics from a NanoAOD ROOT file."""
import argparse
import csv
import math
from collections import Counter, defaultdict

import awkward as ak
import numpy as np
import uproot


def branch(tree, table, name):
    target = table + "_" + name
    hits = [k for k in tree.keys() if k == target or k.endswith("_" + target)]
    if len(hits) != 1:
        raise ValueError("Expected one branch ending %s, found %s" % (target, hits))
    return ak.to_list(tree[hits[0]].array(library="ak"))


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument("root");p.add_argument("--prefix",default="hadr_decay_sanity")
    a=p.parse_args()
    with uproot.open(a.root) as f:
        t=f["Events"]
        h={n:branch(t,"QTrainHadron",n) for n in ("hadronIdx","pdgId","charmOrigin","prodX","prodY","prodZ","prodVertexValid","decayX","decayY","decayZ","decayVertexValid")}
        for n in ("trueLxy", "trueDecayLengthValid"):
            hits=[k for k in t.keys() if k=="QTrainHadron_"+n or k.endswith("_QTrainHadron_"+n)]
            if len(hits)==1: h[n]=ak.to_list(t[hits[0]].array(library="ak"))
        tr={n:branch(t,"QTruthMatchedTrack",n) for n in ("trackIdx","hadronIdx","matchStatus","daughterRelation")}
        pa={n:branch(t,"QTruthMatchedPairFit",n) for n in ("trackIdx1","trackIdx2","hadronIdx","hadronCategory","nMatchedTracks","daughterRelation1","daughterRelation2","commonImmediateDecayVertex","passedOrdinaryPairSelection","vtxFitValid","vtxX","vtxY","vtxZ","vtxLxyValid","vtxLxy","vtxCovarianceValid")}
    rows=[]; totals=Counter(); multiplicity=Counter(); directCascade=Counter(); residuals=defaultdict(list); byMultiplicity=defaultdict(Counter); residualsByN=defaultdict(lambda:defaultdict(list))
    for ie in range(len(h["hadronIdx"])):
        hmap={idx:i for i,idx in enumerate(h["hadronIdx"][ie])}
        for ih,idx in enumerate(h["hadronIdx"][ie]):
            if h["charmOrigin"][ie][ih] != 1: continue
            ntrk=sum(x==idx for x in tr["hadronIdx"][ie]); multiplicity[ntrk]+=1; totals["promptD"]+=1
            if ntrk>=2: totals["promptD_N2plus"]+=1
        for ip,idx in enumerate(pa["hadronIdx"][ie]):
            ih=hmap.get(idx)
            if ih is None or h["charmOrigin"][ie][ih] != 1: continue
            totals["promptD_pairs"]+=1
            if pa["passedOrdinaryPairSelection"][ie][ip]: totals["pairs_in_trackPair"]+=1
            d1,d2=pa["daughterRelation1"][ie][ip],pa["daughterRelation2"][ie][ip]
            topology="direct_common_vertex" if pa["commonImmediateDecayVertex"][ie][ip] else "direct_ancestry_vertex_unresolved" if d1==1 and d2==1 else "contains_cascade" if d1 in (1,2) and d2 in (1,2) else "unresolved"
            directCascade[topology]+=1
            if topology=="direct_common_vertex":
                totals["direct_direct_pairs"]+=1
                if pa["passedOrdinaryPairSelection"][ie][ip]: totals["direct_direct_in_trackPair"]+=1
            nmatched=pa["nMatchedTracks"][ie][ip]
            byMultiplicity[nmatched]["pairs"]+=1
            byMultiplicity[nmatched]["covariance_valid"]+=bool(pa["vtxCovarianceValid"][ie][ip])
            byMultiplicity[nmatched]["fit_success"]+=bool(pa["vtxFitValid"][ie][ip])
            if not pa["vtxFitValid"][ie][ip]: totals["promptD_fit_failed"]+=1; continue
            totals["promptD_fit_success"]+=1
            px,py,pz=h["prodX"][ie][ih],h["prodY"][ie][ih],h["prodZ"][ie][ih]
            dx,dy,dz=h["decayX"][ie][ih],h["decayY"][ie][ih],h["decayZ"][ie][ih]
            fx,fy,fz=pa["vtxX"][ie][ip],pa["vtxY"][ie][ip],pa["vtxZ"][ie][ip]
            valid=bool(h["prodVertexValid"][ie][ih] and h["decayVertexValid"][ie][ih])
            row={"eventEntry":ie,"hadronIdx":idx,"hadronPdgId":h["pdgId"][ie][ih],"nMatchedTracks":nmatched,
                 "trackIdx1":pa["trackIdx1"][ie][ip],"trackIdx2":pa["trackIdx2"][ie][ip],"daughterRelation1":d1,"daughterRelation2":d2,
                 "commonImmediateDecayVertex":pa["commonImmediateDecayVertex"][ie][ip],"pairTopology":topology,"passedOrdinaryPairSelection":pa["passedOrdinaryPairSelection"][ie][ip],"fitValid":1,
                 "fitX":fx,"fitY":fy,"fitZ":fz,"trueProdX":px,"trueProdY":py,"trueProdZ":pz,"trueDecayX":dx,"trueDecayY":dy,"trueDecayZ":dz,
                 "trueDecayGeometryValid":int(valid),"residualX_cm":fx-dx if valid else float("nan"),"residualY_cm":fy-dy if valid else float("nan"),
                 "residualZ_cm":fz-dz if valid else float("nan"),"fitLxyFromTrueProduction_cm":math.hypot(fx-px,fy-py) if valid else float("nan"),
                 "trueLxy_cm":h["trueLxy"][ie][ih] if "trueLxy" in h and h["trueDecayLengthValid"][ie][ih] else float("nan"),
                 "recoPvLxyValid":pa["vtxLxyValid"][ie][ip],"recoPvLxy_cm":pa["vtxLxy"][ie][ip],"covarianceValid":pa["vtxCovarianceValid"][ie][ip]}
            rows.append(row)
            if valid:
                for k in ("residualX_cm","residualY_cm","residualZ_cm"):
                    residuals[k].append(row[k]);residualsByN[nmatched][k].append(row[k])
    csv_path=a.prefix+"_promptD_pairs.csv"
    if rows:
        with open(csv_path,"w",newline="") as f:
            w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    print("Prompt-D hadron matched-track multiplicity:")
    for n in sorted(multiplicity): print("  N=%s: %d"%(n,multiplicity[n]))
    print("Prompt-D pair categories: %s"%dict(directCascade))
    print("Prompt-D pair fit summary: %s"%dict(totals))
    print("Pair selection survival: %.6f"%(totals["pairs_in_trackPair"]/totals["promptD_pairs"] if totals["promptD_pairs"] else 0.))
    print("Direct/direct pair selection survival: %.6f"%(totals["direct_direct_in_trackPair"]/totals["direct_direct_pairs"] if totals["direct_direct_pairs"] else 0.))
    print("Supplementary fit success: %.6f"%(totals["promptD_fit_success"]/totals["promptD_pairs"] if totals["promptD_pairs"] else 0.))
    print("Prompt-D pair covariance validity by matched-track multiplicity:")
    for n in sorted(byMultiplicity):
        c=byMultiplicity[n];print("  N=%d: fit-success %d/%d, covariance-valid %d/%d"%(n,c["fit_success"],c["pairs"],c["covariance_valid"],c["pairs"]))
    summary_path=a.prefix+"_by_multiplicity.csv"
    with open(summary_path,"w",newline="") as f:
        fields=["nMatchedTracks","nPairs","fitSuccess","covarianceValid","meanResidualX_cm","stdResidualX_cm","meanResidualY_cm","stdResidualY_cm","meanResidualZ_cm","stdResidualZ_cm"]
        w=csv.DictWriter(f,fieldnames=fields);w.writeheader()
        for n,c in sorted(byMultiplicity.items()):
            row={"nMatchedTracks":n,"nPairs":c["pairs"],"fitSuccess":c["fit_success"],"covarianceValid":c["covariance_valid"]}
            for k,label in (("residualX_cm","X"),("residualY_cm","Y"),("residualZ_cm","Z")):
                vals=np.asarray(residualsByN[n][k]);row["meanResidual%s_cm"%label]=float(np.mean(vals)) if len(vals) else float("nan");row["stdResidual%s_cm"%label]=float(np.std(vals)) if len(vals) else float("nan")
            w.writerow(row)
    print("Wrote %s"%summary_path)
    print("Wrote %s (%d prompt-D fitted pairs)"%(csv_path,len(rows)))
    if rows:
        try:
            import matplotlib.pyplot as plt
            fig,axes=plt.subplots(1,3,figsize=(13,4))
            for ax,key,label in zip(axes,("residualX_cm","residualY_cm","residualZ_cm"),("x","y","z")):
                vals=np.asarray(residuals[key]); vals=vals[np.isfinite(vals)]
                if len(vals): ax.hist(vals,bins=60,histtype="step");ax.set_xlabel("fitted - true decay %s [cm]"%label);ax.set_ylabel("pairs")
            fig.tight_layout();fig.savefig(a.prefix+"_residuals.png",dpi=150);plt.close(fig)
            print("Wrote %s_residuals.png"%a.prefix)
        except ImportError:
            print("matplotlib unavailable; skipped residual plot")


if __name__=="__main__": main()
