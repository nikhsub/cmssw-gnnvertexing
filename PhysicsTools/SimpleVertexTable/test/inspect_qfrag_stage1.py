#!/usr/bin/env python3
import argparse
import ROOT

parser = argparse.ArgumentParser(description="Summarise Stage-1 QFrag NanoAOD FlatTables")
parser.add_argument("file")
args = parser.parse_args()

root = ROOT.TFile.Open(args.file)
events = root.Get("Events")
counters = {name: 0 for name in ("charm", "bottom", "light", "B", "D", "unique", "ambiguous", "missing", "prompt", "bToD", "unknown", "hfDecay", "fragmentation", "chargedAmbiguous", "unassigned")}
for event in events:
    for value in event.QFrag_flavour:
        if value == 4: counters["charm"] += 1
        elif value == 5: counters["bottom"] += 1
        elif 1 <= value <= 3: counters["light"] += 1
    for value in event.QHadron_flavour:
        if value == 5: counters["B"] += 1
        elif value == 4: counters["D"] += 1
    for value in event.QHadron_associationQuality:
        counters[{0: "missing", 1: "unique", 2: "ambiguous"}.get(value, "missing")] += 1
    for flavour, origin in zip(event.QHadron_flavour, event.QHadron_charmOrigin):
        if flavour == 4:
            counters[{1: "prompt", 2: "bToD", 3: "unknown"}.get(origin, "unknown")] += 1
    for value in event.QGenCharged_membershipType:
        counters[{0: "unassigned", 1: "hfDecay", 2: "fragmentation", 3: "chargedAmbiguous"}.get(value, "unassigned")] += 1

print("QFrag systems: charm=%d bottom=%d light=%d" % (counters["charm"], counters["bottom"], counters["light"]))
print("QHadron: B=%d D=%d" % (counters["B"], counters["D"]))
print("QHadron quark association: unique=%d ambiguous=%d missing=%d" % (counters["unique"], counters["ambiguous"], counters["missing"]))
denom = counters["D"]
print("D origin (D hadrons only): prompt=%d B-to-D=%d unknown=%d" % (counters["prompt"], counters["bToD"], counters["unknown"]))
if denom:
    print("D-origin fractions: prompt=%.4f B-to-D=%.4f unknown=%.4f" % (counters["prompt"]/denom, counters["bToD"]/denom, counters["unknown"]/denom))
print("QGenCharged: HF-decay=%d retained-fragmentation=%d ambiguous=%d unassigned=%d" % (counters["hfDecay"], counters["fragmentation"], counters["chargedAmbiguous"], counters["unassigned"]))
