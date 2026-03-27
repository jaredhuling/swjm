# Flag Last Pseudo-Entry per Subject

Returns a logical vector indicating whether each pseudo-entry is the
last entry for its subject. Used for tie-breaking in timeline
construction.

## Usage

``` r
jfm_compute_is_last(pseudo_entries)
```

## Arguments

- pseudo_entries:

  Sorted pseudo data set (from `jfm_wt_death`).

## Value

Logical vector, one per pseudo-entry row.
