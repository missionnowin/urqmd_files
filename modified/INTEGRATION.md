# Flucton Collision-Physics Patch — Integration Guide

## 1. The physics bug in the existing patch

Your partial patch introduced the flucton as ityp=71 and wired it into `sigmaLN`
row 16 — but it reused **iline numbers 13, 14, 15** to dispatch flucton–N
scattering inside `make22`/`crossx`.

Those iline labels are already used in UrQMD as the generic baryon–baryon
`AQM-elastic / AQM-inelastic / string` branches. So every flucton–N collision
was silently executing the *string* or *AQM* code path, producing a
flucton-plus-nucleon two-body final state. The flucton **never broke up**
into 3 nucleons, and because it decays back to N+N elastically at rest,
the only nucleons produced in the lab inherited roughly half of the
flucton's 4-momentum each — i.e. almost exactly the same kinematics as
ordinary scattering. That is why your backward-hemisphere impulse
distribution looked unchanged.

Backward ("cumulative") nucleons in p+A data (Kopeliovich's "Buddha's
light", Hen/JLab SRC measurements, and the recent arXiv 2509.22433
review) come from the internal nucleon momentum inside the flucton —
p\* ≈ √((M_F/2)² − m_N²) ≈ 0.35 GeV/c — being boosted to the lab by the
flucton's collective motion. You only see it if the flucton actually
**disintegrates** into three (or more) nucleons in the final state.

## 2. Fix at a glance

Introduce three new iline numbers **65 / 66 / 67** reserved for
flucton–N processes and wire dedicated cross-section and exit-channel
code to them.

| iline | process              | cross section         | final state                |
|-------|----------------------|-----------------------|----------------------------|
| 65    | F + N total          | siglookup(13, e)      | (selector only)            |
| 66    | F + N elastic        | siglookup(14, e)      | F + N (2-body)             |
| 67    | **F + N → N+N+N**    | siglookup(15, e)      | **3 nucleons (nbodydec)**  |

Channel **67 is where the cumulative protons are born.** `nbodydec`
performs uniform 3-body phase-space decay of the total (F+N)
4-momentum, populating `pnew(1..4, 1..3)` in the (F+N) rest frame;
`scatfinal` then applies `rotbos(th, phi2, betax, betay, betaz, ...)`
so the boost to the lab frame is already in place.

## 3. File-by-file changes (9 files in `/patch_new/`)

Each file is a drop-in replacement for its counterpart in `/modified/`
(or `/all/` / `/base/` if `/modified/` has no copy).

### 3.1 `comres.f`
- `maxbrf = 2`   (was 1) — now two flucton decay channels
- `nsigs = 16`    (was 12) — sigmaLN now has 16 reaction classes

### 3.2 `blockres.f`
- Flucton mass/width set to **2.0 GeV / 0.200 GeV** (was 0 width — which made
  `anndec` never fire for it). Finite width means `dectim` schedules the
  decay as a finite lifetime instead of instantly.
- `sigmaLN(*, *, 16) = 3, 65, 66, 67, 0...` — dispatches flucton–N through
  the **new** iline labels.
- `sigmainf(15, 3) = -3`   (3-body flag) plus row 15 points at the new
  breakup bookkeeping.
- New data blocks:
  - `bftype(1..4, 0..maxbrf)` — decay products for F → NN and F → NNπ
  - `branfluc(0..maxbrf)` = (1.0, 0.85, 0.15) — 85% NN, 15% NNπ
  - `lbf(0..maxbrf)` = (0, 0, 1) — L of decay products
- `flbr`, `fbran`, `mminit`, `brange`, `b3type` all extended with
  `if(ia.ge.minfluc.and.ia.le.maxfluc)` branches so the flucton is a
  first-class resonance from UrQMD's decay machinery point of view.

### 3.3 `scatter.f`
- New `collclass = 16` when exactly one of {i1, i2} is a flucton and the
  other is an ordinary nucleon.
- Flucton+flucton set to `collclass = 0` (no-op for now — you may revisit
  when you extend to heavy ion where FF collisions become relevant).

### 3.4 `make22.f` — the heart of the patch
Three distinct additions:

1. **Dispatcher `goto` at line 111 extended** with `65,66,67` so that
   `goto(...)io` actually branches to the new labels.
2. **Cross-section dispatcher `goto` at line 4393 extended** the same way
   inside `crossx`.
3. **Three new exit-channel labels (≈ line 3097) and three new crossx
   labels (≈ line 6204).**

Critical new block — the breakup (`io = 67`):

```fortran
 67   continue
c ===== Flucton + N BREAKUP : F + N -> N + N + N =====
      if(e.le.3.d0*massit(minnuc)+1.d-4)then
         call setizm(...); goto 2001   ! below threshold -> elastic
      endif
      nstring1=1; nstring2=2; nexit=3
      itypnew(1..3) = 3 nucleons
      i3new(1..3)   = (+1, +1, iz_target)  for F=pp cluster
      do i=1,3
         pnew(5,i) = massit(iabs(itypnew(i)))
         mstring(i) = pnew(5,i)
      enddo
      call nbodydec(e)   ! uniform 3-body phase space in F+N c.m.
      return
```

`nbodydec` (from `jdecay2.f`) populates `pnew(1..4, 1..nexit)` in the
rest frame of invariant mass `e`. Downstream, `scatfinal` calls
`rotbos(th, phi2, betax, betay, betaz, pnew(1), pnew(2), pnew(3),
pnew(4), nexit)` — the (betax, betay, betaz) are the CM→lab velocities
computed in `scatter.f`. So the 3-body final state is correctly
Lorentz-boosted back into the lab, and the backward-going piece of the
phase-space distribution survives as the cumulative tail.

### 3.5 `anndec.f`
- Added a **flucton decay dispatch** branch
  `else if(iabs(i1).ge.minfluc.and.iabs(i1).le.maxfluc)` which calls
  `anndex(..., maxbrf, minfluc, maxfluc, bftype, branfluc)`. Without
  this, the flucton would sit in the list forever (or crash) when its
  lifetime expires.
- Guard in the meson–baryon annihilation block so a flucton on the i2
  slot is skipped (flucton has no antibaryon meson-annihilation channel).

### 3.6 `cascinit.f`
- **Reverted to `/modified/cascinit.f` verbatim.** No changes. The 7
  fluctons per nucleus that you already place are kept as-is; this patch
  is purely about collision physics.

### 3.7–3.9 `dectim.f`, `ityp2pdg.f`, `itypdata.f`
- Already consistent with the above; included in `/patch_new/` unchanged
  so the package is self-contained.

## 4. Verification done

- All nine Fortran files pass `gfortran -fsyntax-only -ffixed-form
  -std=legacy -w` with no errors and no warnings against the full UrQMD
  source tree (staged from `/base/ + /all/` plus these overlays).
- `cascinit.f` is byte-identical to `/modified/cascinit.f`.
- Computed-goto extensions at both `make22` (line 115) and `crossx`
  (line 4393) reach label 67.
- `sigmaLN` row 16 confirmed to reference 65/66/67, not 13/14/15.

## 5. Expected physics and a test plan

### 5.1 Kinematic expectation
- Internal flucton momentum: `p* = sqrt((M_F/2)² − m_N²) ≈ 0.348 GeV/c`.
- For `p_lab = 5 GeV/c` proton on a stationary flucton:
  `sqrt(s) ≈ 5 GeV` → phase-space max `|p*| ≈ 2 GeV/c` per final nucleon
  in the F+N c.m. After Lorentz boost, O(0.3–1.0 GeV/c) backward
  protons appear in the lab — well above the Fermi cutoff (≈ 0.25 GeV/c)
  that limits ordinary cascade.

### 5.2 Suggested smoke test (matches your "large ions" plan)
1. Run `p + Pb` (or `p + Au`) at 5–10 GeV/c lab momentum, ~10 000 events.
2. Select **final-state protons with `p_z < 0`** in the target rest
   frame and histogram `|p_z|`.
3. Before the patch you see only the Fermi-momentum shoulder
   (< 0.3 GeV/c). After the patch you should see a **distinct tail out
   to ~1 GeV/c** that scales with the number of fluctons (currently 7
   per nucleus — you can later scale this with A for heavy systems).
4. Cross-check: the flucton–N breakup cross section `siglookup(15, e)`
   must be non-zero at your running `e`. Print its value once at
   initialization to confirm.
5. Multiplicity check: total baryon number must increase by exactly
   `+2` per flucton that breaks up (F → NNN vs. F which decays later
   anyway to NN). Monitor `nbar` before/after.

### 5.3 Tuning knobs
- `branfluc(1) = 0.85`, `branfluc(2) = 0.15` — pion-emission fraction.
- Flucton width `0.200` GeV — sets the mean flucton lifetime; smaller
  width keeps the flucton alive through more scatterings (more chances
  to break up) but may suppress late-time contributions.
- `numfluc` (you set 7) — linear scaling of the cumulative tail.
- The breakup cross section itself (`siglookup(15, e)`) is currently a
  table lookup; if the tail is too weak you can rescale it in
  `sigtab.f` / `tabinit.f`.

## 6. What to do next

1. Drop the 9 files from `/patch_new/` into your build tree (they
   overwrite the corresponding files in `/modified/`).
2. Rebuild UrQMD.
3. Run the smoke test above.
4. If you want, I can commit these to `missionnowin/urqmd_files` via
   `gh` on a new branch (e.g. `flucton-breakup`) with a clear commit
   message and push it — just say the word.

## 7. References

- V. Kopeliovich, "Buddha's light and the cumulative effect",
  http://quarks.inr.ac.ru/2014/proceedings/www/p5/Kopelio.pdf
- Flucton/SRC review, arXiv:2509.22433,
  https://arxiv.org/pdf/2509.22433
- UrQMD collaboration, standard 3-body phase space via `nbodydec` in
  `jdecay2.f`.
