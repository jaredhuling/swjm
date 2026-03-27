# Extract Event Pseudo-Entry Indices

For each event time, extracts the pseudo-entry row index for the subject
who experienced the event. Returns a 0-based integer vector suitable for
passing to C++ score functions.

## Usage

``` r
jfm_event_pseudo_idx(index_matrix, event_subject_ids)
```

## Arguments

- index_matrix:

  The n x ne index matrix (1-based row indices into pseudo_entries),
  either `index_death_matrix` or `index_recurrent_matrix`.

- event_subject_ids:

  Integer vector (length ne) of 1-based subject IDs who experienced each
  event (e.g., td_id or tr_id).

## Value

Integer vector (length ne) of 0-based pseudo-entry row indices.
