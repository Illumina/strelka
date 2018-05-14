# Sequence allele counting terminology

## Terms shared between data gathering and parameter estimation steps:

- **CONTEXT** The pattern which is used to stratify observation counts. Most often this is the contextual reference
sequence pattern around the variant location, but this can also include encoding of the alternate allele type as well.

- **CONTEXT INSTANCE** One occurrence of the CONTEXT.

- **OBSERVATION** The sequence patterns observed at a single CONTEXT INSTANCE.


## Terms only relevant to compression of many observations

- **OBSERVATION PATTERN** A sequence pattern which is theoretically observable at zero to many CONTEXT INSTANCEs.

- **OBSERVATION PATTERN COUNT** The number of CONTEXT INSTANCEs at which the OBSERVATION PATTERN was present


## Example

Considering the following simple set of contexts:

- There are 4 CONTEXTs for basecall observations, for positions where the reference is either A, C, G, T
    - Theses contexts are labeled refA, refC, refG, refT

And the example set of reference and reads:

```
Ref:
AACGTA

Reads:
AAT
GAT
 CTG
  TGT
   GTG
   GTA
```

Positions 1, 2 and 6 are 3 different INSTANCE's of the refA CONTEXT.

OBSERVATION PATTERN (A,G) is present at two CONTEXT INSTANCES or CONTEXT refA, at positions 1 and 6.
