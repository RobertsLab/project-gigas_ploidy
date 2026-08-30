# Known Issues in MethylKit Analysis

Issues identified in `bisulfite_analysis/WGBS/code/2_WGBS_Methylkit.R`. This file tracks bugs and methodological concerns to address before final publication.

---

## 🔴 Bugs (script will error or produce wrong output on re-run)

### BUG-1: Undefined variable `fileName` in plot filename generation
**Lines:** 108–112  
**Symptom:** Script errors with `object 'fileName' not found` when generating per-sample plots  
**Cause:** `fileNameCov5` and `fileNameCov10` were defined above (lines 101–106), but `fileName$nameBase` etc. are used instead.  
**Fix:**
```r
# Replace fileName$nameBase with fileNameCov5$nameBase (and fileNameCov10$nameBase for Cov10 block)
fileNameCov5$actualFileName1 <- paste(fileNameCov5$nameBase, "-Filtered", "-5xCoverage", "-Sample", fileNameCov5$sample.ID, ".jpeg", sep = "")
```

---

### BUG-2: Undefined variable `differentialMethylationStatsTreatment` — results depend on pre-saved .RData
**Lines:** 209, 220, 225, 354, 479  
**Symptom:** Script errors if run fresh; DML output silently comes from `load("20210528_methylKit_DMLs_cov10.RData")` (line 44, Windows path)  
**Cause:** The correct variable is `differentialMethylationStatsTrtCov10` / `diffMethStatsTrt50Cov10` throughout.  
**Fix:** Replace all occurrences:
- `differentialMethylationStatsTreatment` → `differentialMethylationStatsTrtCov10`
- `diffMethStatsTreatment50` → `diffMethStatsTrt50Cov10`

**Additional action needed:** Document or upload the `.RData` session file to gannet so results are reproducible.

---

### BUG-3: Coverage level mismatch — Cov10 DMLs subsetted from Cov5 matrix for PCA/heatmap
**Lines:** 354–361  
**Symptom:** PCA and heatmap may show incorrect loci — `which()` matching Cov10 DML positions against Cov5 CpG positions can return wrong indices or `integer(0)`  
**Fix:**
```r
# Change Cov5 references to Cov10:
DMLPositions[i] <- which(getData(diffMethStatsTrt50Cov10)$start[i] == getData(methylationInformationFilteredCov10)$start)
DMLMatrix <- methylationInformationFilteredCov10[DMLPositions,]
```
**Action:** Regenerate `plots/2021-05-26-DML-Only-PCA.pdf` and `plots/2021-05-27-DML-Only-Heatmap.pdf` after fix.

---

## 🟠 Wrong metadata/labels (results exist but are mislabeled)

### BUG-4: CpG counts per chromosome use *C. virginica* values, not *C. gigas*
**Lines:** 260–276  
**Evidence:** Comment at line 274 reads "Need to update these numbers for c.gigas"; the bash snippet in the comment explicitly references `C_virginica-3.0_CG-motif.bed`  
**Impact:** The DML-per-CpG normalization in `plots/2021-05-25-DML-and-Gene-Distribution.pdf` uses incorrect denominators — the relative enrichment per chromosome is wrong.  
**Fix:** Re-count CpGs per chromosome from the Roslin v1 genome used in this study.

---

### BUG-5: Heatmap and PCA legend labels copied from a pCO2/ocean acidification study
**Lines:** 459, 445–446  
**Evidence:**
```r
annotation_col = data.frame(pCO2 = factor(rep(c("Ambient","Treatment"), each = 5)))  # line 459
legend = c("Control", "Elevated")  # line 445
```
**Impact:** Heatmap annotations and PCA legend describe diploid as "Ambient/Control" and triploid as "Treatment/Elevated" — misleading for this ploidy comparison study.  
**Fix:** Replace with `Ploidy = factor(rep(c("Diploid","Triploid"), each = 5))` and update legend labels.

---

## 🟡 Methodological concerns (may affect result interpretation)

### CONCERN-1: No overdispersion correction in `calculateDiffMeth()`
**Lines:** 193–195  
**Issue:** Default logistic regression used; for WGBS data with biological replicates, betabinomial (`overdispersion = "MN"`) is more appropriate and reduces false positives.  
**Context:** The commented-out line 197 shows this was considered. The 60,721 DMLs at ≥20% diff may be inflated.  
**Recommendation:** Rerun with `overdispersion = "MN", test = "Chisq"` and compare DML counts.

---

### CONCERN-2: `destrand = FALSE` in `unite()` — CpGs not destranded
**Lines:** 150–160  
**Issue:** Standard practice for CpG methylation is to destrand (merge +/- strand CpGs at the same position). Using `destrand = FALSE` keeps them separate, treating the same CpG dinucleotide as two loci.  
**Recommendation:** Consider rerunning with `destrand = TRUE` to follow convention, and check whether DML count changes substantially.
