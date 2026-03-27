# Compute Exit Times for Pseudo-Entries

For each row in the sorted pseudo-data matrix, computes the exit time:
non-last entries for a subject exit at the next entry's start time; the
last entry exits at Y_i (the subject's observation time).

## Usage

``` r
jfm_compute_exit_times(pseudo_entries, Y)
```

## Arguments

- pseudo_entries:

  Sorted pseudo data set (from `jfm_wt_death`).

- Y:

  Numeric vector of observation times per subject (length n).

## Value

Numeric vector of exit times, one per pseudo-entry row.
