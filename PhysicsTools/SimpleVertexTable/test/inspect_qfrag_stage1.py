#!/usr/bin/env python3
import argparse
import ROOT

parser = argparse.ArgumentParser(description="Summarise Stage-1 QFrag NanoAOD FlatTables")
parser.add_argument("file")
args = parser.parse_args()

root = ROOT.TFile.Open(args.file)
events = root.Get("Events")
counters = {name: 0 for name in ("charm", "bottom", "light", "B", "D", "unique", "ambiguous", "missing", "prompt", "bToD", "unknown", "hfDecay", "fragmentation", "chargedAmbiguous", "unassigned", "lineageUnique", "lineageAmbiguous", "lineageMissing", "lineageUnknown", "lineagePrompt", "lineageBtoD", "lineageOriginUnknown")}
old_new = {}
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
    for old, new in zip(event.QHadron_associationQuality, event.QHadron_lineageAssociationQuality):
        old_new[(old, new)] = old_new.get((old, new), 0) + 1
        counters[{0: "lineageMissing", 1: "lineageUnique", 2: "lineageAmbiguous", 3: "lineageUnknown"}.get(new, "lineageUnknown")] += 1
    for flavour, origin in zip(event.QHadron_flavour, event.QHadron_charmOrigin):
        if flavour == 4:
            counters[{1: "prompt", 2: "bToD", 3: "unknown"}.get(origin, "unknown")] += 1
    for flavour, origin in zip(event.QHadron_flavour, event.QHadron_lineageOriginCategory):
        if flavour == 4:
            counters[{1: "lineagePrompt", 2: "lineageBtoD", 3: "lineageOriginUnknown"}.get(origin, "lineageOriginUnknown")] += 1
    for value in event.QGenCharged_membershipType:
        counters[{0: "unassigned", 1: "hfDecay", 2: "fragmentation", 3: "chargedAmbiguous"}.get(value, "unassigned")] += 1

print("QFrag systems: charm=%d bottom=%d light=%d" % (counters["charm"], counters["bottom"], counters["light"]))
print("QHadron: B=%d D=%d" % (counters["B"], counters["D"]))
print("QHadron quark association: unique=%d ambiguous=%d missing=%d" % (counters["unique"], counters["ambiguous"], counters["missing"]))
print("QHadron alternative lineage: unique=%d ambiguous=%d missing=%d unknown=%d" % (counters["lineageUnique"], counters["lineageAmbiguous"], counters["lineageMissing"], counters["lineageUnknown"]))
print("Old/new association quality counts (old,new): %s" % sorted(old_new.items()))
denom = counters["D"]
print("D origin (D hadrons only): prompt=%d B-to-D=%d unknown=%d" % (counters["prompt"], counters["bToD"], counters["unknown"]))
print("Alternative D lineage origin: prompt=%d B-to-D=%d unknown=%d" % (counters["lineagePrompt"], counters["lineageBtoD"], counters["lineageOriginUnknown"]))
if denom:
    print("D-origin fractions: prompt=%.4f B-to-D=%.4f unknown=%.4f" % (counters["prompt"]/denom, counters["bToD"]/denom, counters["unknown"]/denom))
print("QGenCharged: HF-decay=%d retained-fragmentation=%d ambiguous=%d unassigned=%d" % (counters["hfDecay"], counters["fragmentation"], counters["chargedAmbiguous"], counters["unassigned"]))
