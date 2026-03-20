# Counting Process of Recurrent Event

Computes the value of the counting process N_i(t), i.e., the number of
recurrent events that occurred at or before time t.

## Usage

``` r
jfm_counting_process(t, recur_time)
```

## Arguments

- t:

  A scalar time point.

- recur_time:

  A numeric vector of recurrent event times for a subject.

## Value

Integer count of recurrent events at or before time t.
