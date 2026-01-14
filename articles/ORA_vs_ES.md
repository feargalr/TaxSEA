# Analysis types

## ORA vs Enrichment (two ways to run TaxSEA)

TaxSEA supports two complementary analysis modes: **Enrichment** and
**ORA (Over-Representation Analysis)**.  
They answer related but distinct questions and are appropriate in
different situations.

## Enrichment (rank-based)

### What it does (in plain terms)

You provide TaxSEA with a **ranked list of taxa** (for example from most
increased to most decreased).  
TaxSEA then asks:

> Do the members of a given taxon set tend to appear unusually high or
> low in this ranking?

Importantly, **all taxa are used**, not just those passing a
significance threshold.

### When to use it (recommended in most cases)

Use **Enrichment** when:

- You have differential abundance results (e.g. log fold changes, test
  statistics, signed p-values)
- Taxa can be meaningfully ordered by strength and direction of effect
- You want to avoid hard cutoffs such as “significant vs not
  significant”

### Pros

- Uses all taxa, not only a selected subset  
- Less sensitive to arbitrary thresholds  
- More stable when signal is spread across many taxa  
- Can detect subtle, coordinated shifts across a taxon set

### Cons

- Requires a sensible ranking metric  
- Not appropriate if the data are strictly presence/absence with no
  meaningful ordering

## ORA (Over-Representation Analysis)

### What it does (in plain terms)

You provide TaxSEA with a **list of taxa of interest** (your “hits”),
along with a background universe.  
TaxSEA then asks:

> Are taxa from this set appearing in my hit list more often than
> expected by chance?

This is a classic enrichment-of-hits approach using contingency tables.

### When to use it

Use **ORA** when:

- Your data are naturally **binary** (presence vs absence, detected vs
  not detected)
- A ranked statistic is not meaningful or not available
- You have a strong, biologically justified threshold defining your taxa
  of interest

### Pros

- Simple and intuitive interpretation  
- Works well for presence/absence style data  
- Appropriate when only a hit list is available

### Cons

- Highly sensitive to how “hits” are defined  
- Discards information about effect size and ordering  
- Results can change substantially with small changes to thresholds

## Recommended practice

In most microbiome differential abundance analyses, **Enrichment is
recommended** whenever a reasonable ranking can be constructed (e.g. log
fold change, Wald statistic, signed −log10 p-value). It makes fuller use
of the data and avoids arbitrary cutoffs.

**ORA** should be used when ranking is not possible or not meaningful,
particularly in **presence/absence** scenarios or when the scientific
question explicitly concerns membership in a predefined hit list rather
than graded shifts.

Both methods are valid; the choice depends on the structure of the data
and the biological question being asked.
