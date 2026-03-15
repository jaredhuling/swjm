# swjm: Stagewise Variable Selection for Joint Models of Semi-Competing Risks

Implements stagewise regression for penalized variable selection in
joint models of recurrent events and terminal events (semi-competing
risks).

Two model frameworks are supported:

- **Joint Frailty Model (JFM)**: Cox-type model where `alpha` =
  recurrence/readmission coefficients (first p elements of theta),
  `beta` = terminal/death coefficients (second p elements).

- **Joint Scale-Change Model (JSCM)**: AFT-type model where `alpha` =
  recurrence/readmission coefficients (first p elements of theta),
  `beta` = terminal/death coefficients (second p elements).

Three penalty types are available:

- `"coop"`: Cooperative lasso, which encourages shared support between
  the recurrence and terminal event coefficient vectors.

- `"lasso"`: Standard L1 penalty applied coordinate-wise.

- `"group"`: Group lasso, which selects variables jointly across both
  processes.

## Data Format

All functions expect a recurrent-event data frame with the following
columns:

- id:

  Subject identifier (integer).

- t.start:

  Interval start time.

- t.stop:

  Interval stop time.

- event:

  1 for recurrent event, 0 for terminal/censoring row.

- status:

  1 for death, 0 for alive/censored.

- x1, x2, ..., xp:

  Covariate columns.

## Author

**Maintainer**: Huo Lingfeng <huol@example.com>

Authors:

- Huling Jared
