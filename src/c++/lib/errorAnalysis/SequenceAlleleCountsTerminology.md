# Sequence allele counting terminology

## Terms shared between data gathering and parameter estimation steps:

- **CONTEXT** The pattern which is used to stratify observation counts. Most often this is the contextual reference
sequence pattern around the variant location, but this can also include encoding of the alternate allele type as well.
  - One example of a context which depends on alternate allele type is a homopolymer expansion. In this case the context
consists of a homopolymer sequence in the reference (eg. 20-base "A" homopolymer) and an alternate allele which
expands the homopolymer (a one-base insertion of "A"). Insertions of other bases in the homopolymer tract would not be
enumerated in the same context.

- **CONTEXT INSTANCE** One occurrence of the CONTEXT.

- **OBSERVATION** The sequence patterns observed at a single CONTEXT INSTANCE.


## Terms only relevant to compression of many observations

- **OBSERVATION PATTERN** A sequence pattern which is theoretically observable at zero to many CONTEXT INSTANCEs.

- **OBSERVATION PATTERN COUNT** The number of CONTEXT INSTANCEs at which the OBSERVATION PATTERN was present


## Example

Considering the following simple set of contexts:

- There are 4 CONTEXTs for basecall observations, for positions where the reference is either A, C, G or T
    - Theses contexts are labeled _refA_, _refC_, _refG_ and _refT_

And the example set of reference and reads:

```
Ref:
AACGTAG

Reads:
AAC
GAC
 ATG
  TAT
   GTGA
   GTAG
```

This is decomposed into OBSERVATION PATTERNS as follows:

- _refA_ context:
  - OBSERVATION PATTERN {A,G} is present at two CONTEXT INSTANCEs of the _refA_ CONTEXT, at positions 1 and 6
  - OBSERVATION PATTERN {A,A,A} is present at one CONTEXT INSTANCE of the _refA_ CONTEXT, at position 2
- _refC_ context:
  - OBSERVATION PATTERN {C,C,T,T} is present at one CONTEXT INSTANCE of the _refC_ CONTEXT, at position 3
- _refG_ context:
  - OBSERVATION PATTERN {A,G,G,G} is present at one CONTEXT INSTANCE of the _refG_ CONTEXT, at position 4
  - OBSERVATION PATTERN {A,G} is present at one CONTEXT INSTANCE of the _refG_ CONTEXT, at position 7
- _refT_ context:
  - OBSERVATION PATTERN {T,T,T} is present at one CONTEXT INSTANCE of the _refT_ CONTEXT, at position 5
