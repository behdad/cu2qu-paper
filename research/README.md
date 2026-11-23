# Earlier tools

Scripts and output from the work that led to the note: the HarfBuzz, fontTools
and kurbo changes, and research runs by AI agents on other converters, error
measures and weights. They are kept as they were. Paths in them point to
where they ran, and most will not run elsewhere without editing.

None of this is needed to check the note; `scripts/` is. And none of it was
reviewed the way the note was. Where a file here gives numbers for other
software (kurbo, lyon, Skia, FontForge, cubic2quad), some come from that
software and some from reimplementations written for the comparison; they
should not be quoted as results about it.

- `closed-form/`: first checks of the closed form, and the first look at
  other weights.
- `harfbuzz/`, `fonttools/`, `kurbo/`: tools for the three sets of changes.
- `survey/`, `error/`, `weights/`: the three research runs summarized in
  `notes/research-agents.md`.
- `kinks/`, `fit-test/`: the first measurement of kinks on masters, and the
  study of the fit test that does not return.
