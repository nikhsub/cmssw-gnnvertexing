#!/usr/bin/env python3
"""Summarize the compact Stage-1 QFrag, QHadron, and QGenCharged tables."""
import argparse
from collections import Counter

import ROOT

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("file")
args = parser.parse_args()
root_file = ROOT.TFile.Open(args.file)
if not root_file or root_file.IsZombie():
    raise SystemExit("Could not open ROOT file: " + args.file)
events = root_file.Get("Events")
if not events:
    raise SystemExit("ROOT file has no Events tree")

counts = Counter()
for event in events:
    counts["quarks"] += len(event.QFrag_pdgId)
    for pdg in event.QFrag_pdgId:
        flavour = abs(pdg)
        if flavour == 4:
            counts["charm"] += 1
        elif flavour == 5:
            counts["bottom"] += 1
        elif 1 <= flavour <= 3:
            counts["light"] += 1
    for pdg, origin, quality, quark in zip(event.QHadron_pdgId, event.QHadron_charmOrigin,
                                           event.QHadron_associationQuality, event.QHadron_quarkIdx):
        counts["B" if abs(pdg) in (511, 521, 531, 541, 5122, 5132, 5232, 5332, 5142, 5242, 5342, 5512, 5532, 5542, 5554) else "D"] += 1
        counts["hadronOrigin_%d" % origin] += 1
        counts["hadronQuality_%d" % quality] += 1
        counts["hadronMapped"] += quark >= 0
    for quality, membership in zip(event.QGenCharged_associationQuality, event.QGenCharged_membershipType):
        counts["chargedQuality_%d" % quality] += 1
        counts["membership_%d" % membership] += 1
    counts["charged"] += len(event.QGenCharged_pdgId)

print("QFrag: charm=%d bottom=%d light=%d total=%d" %
      (counts["charm"], counts["bottom"], counts["light"], counts["quarks"]))
print("QHadron: B=%d D=%d mapped=%d" % (counts["B"], counts["D"], counts["hadronMapped"]))
print("QHadron charm origin enum: 0=n/a 1=prompt 2=B-to-D 3=unknown; counts=" +
      str({key[len("hadronOrigin_"):]: value for key, value in counts.items() if key.startswith("hadronOrigin_")}))
print("QHadron association quality enum: 0=missing 1=unique 2=ambiguous 3=unknown; counts=" +
      str({key[len("hadronQuality_"):]: value for key, value in counts.items() if key.startswith("hadronQuality_")}))
print("QGenCharged: total=%d" % counts["charged"])
print("QGenCharged membership enum: 0=unassigned 1=heavy-hadron decay 2=additional candidate 3=ambiguous; counts=" +
      str({key[len("membership_"):]: value for key, value in counts.items() if key.startswith("membership_")}))
print("QGenCharged association quality enum: 0=missing 1=unique 2=ambiguous 3=unknown; counts=" +
      str({key[len("chargedQuality_"):]: value for key, value in counts.items() if key.startswith("chargedQuality_")}))
root_file.Close()
