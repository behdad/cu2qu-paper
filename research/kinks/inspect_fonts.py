from fontTools.ttLib import TTFont

for path, num in [
    ("/usr/share/fonts/abattis-cantarell-vf-fonts/Cantarell-VF.otf", -1),
    ("/usr/share/fonts/google-noto-sans-cjk-vf-fonts/NotoSansCJK-VF.ttc", 0),
    ("/usr/share/fonts/google-noto-serif-cjk-vf-fonts/NotoSerifCJK-VF.ttc", 0),
]:
    f = TTFont(path, fontNumber=num) if num >= 0 else TTFont(path)
    print(path)
    print("  tables:", sorted(f.keys())[:40])
    print("  upem", f["head"].unitsPerEm, "glyphs", len(f.getGlyphOrder()))
    for a in f["fvar"].axes:
        print("  axis", a.axisTag, a.minValue, a.defaultValue, a.maxValue)
    if "avar" in f:
        print("  avar", {k: dict(v) for k, v in f["avar"].segments.items()})
    cff2 = f["CFF2"].cff
    td = cff2.topDictIndex[0]
    vs = td.VarStore.otVarStore
    print(
        "  regions:",
        [
            [(r.StartCoord, r.PeakCoord, r.EndCoord) for r in reg.VarRegionAxis]
            for reg in vs.VarRegionList.Region
        ],
    )
    print("  vardata:", [(vd.VarRegionCount, vd.VarRegionIndex) for vd in vs.VarData][:5])
