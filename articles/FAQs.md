# Frequently Asked Questions

## What problem does TaxSEA actually solve?

Traditional differential abundance looks for individual microbes that
change. TaxSEA looks for groups of microbes with shared biology that
change together.

That shift matters because human microbiomes are highly individual —
chasing single species is noisy and rarely reproducible. TaxSEA helps
you detect the biology that’s consistent across people.

## Is TaxSEA a replacement for differential abundance?

No, it’s a complement, not a replacement.

In practice, DA often tells you where to look. TaxSEA helps you
understand what it means.

## Can I use TaxSEA even when there are no individual species level differences?

No, it’s a complement, not a replacement.

In practice, DA often tells you where to look. TaxSEA helps you
understand what it means.

## How is this different from functional profiling tools like HUMAnN3?

They answer different questions. - HUMAnN3 → what genes and pathways are
present? - TaxSEA → what known biological traits or signatures are
shifting at the taxon level?

Sometimes those overlap. Sometimes they don’t as we don’t always have
good gene annotations.

TaxSEA is especially useful when traits are known biologically but not
well captured by gene annotations.

## Does this only work for the human gut microbiome?

That’s where it works best right now, because most curated databases
focus there. However since we introduced the BacDive database it has
more utility in non-human microbiomes.

But the framework itself is general if you have: - taxa - ranks (logFC,
correlations, etc.) - meaningful trait groupings

you can use TaxSEA in other systems too.

## Do I need differential abundance results first?

Not necessarily.TaxSEA just needs a ranked list of taxa: - log2 fold
changes - correlations with a phenotype - effect sizes from any
statistical model

## Why use ranks instead of thresholds?

Other taxon set enrichment approaches exist, however they use a
threshold to define taxa of interest However, thresholds throw away
information.

In microbiome data especially, where signals are weak and noisy,
arbitrary cut-offs often decide the biology. For example, if you have 10
short chain fatty acid producers all decreased but the false discovery
rate is 0.06 for all of them you probably have some interesting biology
but are not powered to detect it with individual taxa.

Rank-based enrichment keeps all taxa in play and asks whether meaningful
groups drift up or down together.

## Isn’t this just gene set enrichment, but for microbes?

Yes, and that’s the point.

Genomics moved past single-gene thinking years ago. Microbiome analysis
is just catching up.

TaxSEA brings the same logic of set-level biology to microbial
communities.

## Will this work if my DA results are messy or inconsistent?

That’s exactly when TaxSEA is most useful.

If different DA methods give you different taxa but the same kinds of
microbes keep appearing (e.g. oxygen-tolerant, SCFA producers,
disease-associated), TaxSEA will surface that shared signal.

## Have you shown TaxSEA is accurate?

Yes, we provided a rigorous evaluation of accuracy in our original
publication

[Pham et
al. 2025](https://academic.oup.com/bib/article/26/2/bbaf173/8116684)

No black boxes, no heavy modelling assumptions. Just a clean question:
are biologically related microbes shifting together more than expected
by chance?

## When should I not use TaxSEA?

TaxSEA is fast and cheap to run, so in practice it’s almost always worth
trying. That said, there are situations where it’s not the best tool.

TaxSEA tends to work less well when: - Diversity is very low. For
example in early-life microbiomes (e.g. infants), where there are simply
too few taxa for meaningful set-level patterns to emerge.

- The system is highly controlled. n pre-clinical models microbiomes are
  much more consistent between individuals (mice for example are
  coprophagic so communities homogenise). In these settings, classic
  differential abundance often works better.

- Your hypothesis is about a single organism. If you’re essentially
  studying something closer to an infection model, DA is the right tool.

TaxSEA really shines when data are heterogeneous (e.g. adult human gut
microbiomes), and signals are subtle. For example when overall diversity
shifts, but no single taxon stands out as “significant”.

## Can I use TaxSEA with 16S rRNA gene sequencing data?

Technically, yes but it works best with shotgun metagenomic data.

TaxSEA is designed around species- and strain-level resolution. With 16S
data you’re usually limited to higher taxonomic ranks, which means you
lose a lot of biological resolution, and many meaningful traits
(e.g. metabolite production) can’t be assigned confidently. That said,
we’ve included genus-level taxon sets in the database, so TaxSEA will
still run and can be useful for high-level patterns. Just be aware that
the results will be coarser and generally less powerful than when using
shotgun metagenomics.

## Why does TaxSEA use species-level data instead of strain-level?

We deliberately work at the species level to bring some consistency to
very heterogeneous human microbiome data.

Strain-level profiles are even more variable across people than species.
In practice, that makes enrichment-style analyses noisy, slow, and hard
to reproduce. So we make a trade-off: we lose some fine-grained detail,
but we gain stability, speed, and interpretability.

There’s also a practical reason. For many of the traits TaxSEA uses, the
underlying evidence either is conserved at the species level
(e.g. oxygen tolerance, broad metabolic capacity, prior disease
associations), or comes from experiments on one or a few type strains.
So even if we wanted to work fully at strain resolution, the biological
annotation often isn’t there.

Yes, this means we may miss some strain-specific effects in individual
samples. But for population-scale human microbiome studies,
species-level analysis strikes a better balance between biological
meaning and practical usability.
