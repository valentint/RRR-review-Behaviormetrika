
<!-- README.md is generated from README.Rmd. Please edit that file -->

# RRR-review-Behaviormetrika

**Replication scripts for the paper [“N. Trendafilov, M. Gallo, V.
Simonacci and V. Todorov: Overview of reduced-rank regression with dense
and sparse coefficients, and a new estimation procedure”]()**

- Clone the repository, resp. unzip the ZIP file in a directory of your
  choice. Two empty sub-directories will be created:

  - `Output` to contain the pdf files of the graphs and
  - `Results` to contain the results of the simulation.

- Open the R project RRR in RStudio, or start R in this directory. This
  will set up the repository as your working directory in R.

- Content of the directory:

  - `Readme`: this file
  - `rrsim.pdf`: description and complete results of the simulation
    study from Section 5.
  - `ccsrrr.R`: contains the R functions `ccsrrr()` and `cv.ccsrrr()`
  - `tobacco.csv`: the Tobacco data set (Izenman, 2008, 6.3.3)
  - `ex_tobacco.R`: an R script to reproduce all examples with the
    Tobacco data set. This script will create the Tables 1, 2 and 3.
  - `dosim2.R`: a script to conduct the simulation study from Section 5.
  - `getResults2.R`: A script to present the results of the simulation
    study (stored in sub-directory `Results`) as Latex tables and plots.
    The plot pdf files are stored in teh sub-directory `Output`.
  - `ex_yeast.R`: a script to reproduce the examples with the Yeast cell
    cycle data set described in Section 6.
  - `timing-microbenchmark.R`: timing of `ccsrrr()` and `secure.path()`
    functions using the package ***`microbenchmark`***.

- Running the simulation.

  - The setup of the simulation is described in the supplementary
    document `rrrsim.pdf`.
  - It can be run separately for each of the four scenarios, each of the
    three error distributions (only the results for the Gaussian errors
    are presented) and for each of the ten methods: `spls`, `srrr`,
    `sarrs`, `rssvd`, `mrce`, `remmap`,
    secure`,`siar`,`sofar`and`ccsrrr\`\` the simulation is run
    separately.
  - That is, set `case=`, `distr=` and `method=`, e.g.`case=1`,
    `distr="norm"` and `method="spls"` and run the script. The result
    will be stored in a file with a name of the following format:
    `method_distr_nXXX_pXX_qXX.csv`, e.g. `spls_norm_n100_p25_q25.csv`.
  - ***Note***: not all methods can be run for all scenarios, because of
    time constraints. All methods run for scenarios 1 and 2, `remmap`
    and `mrce` are skipped for scenario 3 and `remmap`, `mrce` and
    `rssvd` are skipped for scenario 4.
  - Once all methods are run for a given scenario, the presentation of
    the results can be generated. Run `getResults2.R` selecting the
    scenario (1, 2, 3 or 4) and distribution (it defaults to “norm” for
    Gaussian errors). A latex table with the results will be generated
    and five plots will be created in the sub-directory `Results`.

- Running the Tobacco data set examples. The script `ex_tobacco.R` will
  read the file `tobacco.csv` and will generate three tables - also in
  Latex format: Table 1, Table 2 and Table 3.

- Running the Yeast cell cycle data example. The script `ex_yeast.R`
  reads the data from `yeast` from the package ***`spls`*** and runs
  several sparse RRR algorithms on it, including `ccsrrr`. The results -
  how many TFs were selected and how many of these are among the
  confirmed ones - are stored into a table, also in Latex format. The
  selected TFs by `ccsrrr` which are among the confirmed ones are ploted
  into a graph which is stored into the `Output` sub-directory.
