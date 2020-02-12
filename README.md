# Estimating viability and predicting yield in Dilution to Extinction Cultures.


Update to the Button et al. 1993 method for calculating viability of cultures.


## How to set up this codebase

**Step one:** Download this codebase to your machine

```
git clone https://github.com/thrash-lab/viability_test.git
```

**Step two:** Setup a conda environment to run it in

```
conda create -n viability python=3.6 pandas numpy tqdm matplotlib
```

**Step three:** Change to your conda environment
```
conda activate viability
```

**Step four:** Check the code runs
```
python viability_test/viability_test.py -h
```

## How to use this code for your Dilution to Extinction culturing

The pipeline contains three methods:

1. A method to estimate how many positive and pure wells you will get from a DTE experiment for a given taxon
2. A method to evaluate the viability of a culture based on the number of observed wells in a DTE experiment
3. A method to evaluate the viability of a set of culture experiments from a TSV file.

### Estimating yield from a DTE experiment

Let's say you are going to run an experiment whereby you inoculate 500 wells with 2 cells per well from a mixed community. Your taxon of interest makes up 50% of your community (based on 16S rRNA amplicon data). You want to know how many wells of that particular taxon you are likely to observe, assuming the cells are 100% viable:

```
python viability_test/viability_test.py predict_wells -w 500 -i 2 -r 0.05 -v 1
```

Unless you specify `--threads` or `-p`, it will use all the processors on the available machine to peform its bootstrapping.

After a few seconds, you should see something like the following:
```
    You simulated a DTE where you inoculated 500 wells with an average of 2 cells per well.
    Your estimated relative abundance for your taxon of interest was:       50.000%
    Your estimated viability of your taxon of interest was:                 100.000%

    After 9999 simulations:
    The median number of positive (not taxon-specific) wells was:           432 (417-447, 95%CI)
    The median number of wells inoculated with one cell was:                135 (116-155, 95%CI)
    The median number of wells containing your taxon of interest was:       316 (295-337, 95%CI)
    The median number of pure wells for your taxon of interest was:         116 (98-135, 95%CI)
```

So, your experiment will yield between 98 and 135 pure wells for your taxon of interest.


### Evaluating viability from a DTE experiment
Let's say you have performed a DTE experiment and you have observed 5 positive wells for your taxon of interest, instead of the 98-135 wells you were expecting from above, assuming a viability of 100%. You can use the method `estimate_viability` to estimate what viability of your taxon would explain the much lower-than-expected yields:

(note, this function is doing a LOT of bootstrapping, so be patient, or throw it on a large box).

```
python viability_test/viability_test.py estimate_viability -w 500 -i 2 -o 5 -r 0.5
```

Your output should look something like this:
 ```
        You simulated a DTE where you inoculated 500 wells with an average of 2 cells per well.
        A range of decreasing viability was tested using 9999 bootstraps per experiment.
        Your estimated relative abundance for your taxon of interest was:       50.000%
        You observed:                                                           5 wells of interest

        If viability of your taxon had been 100%, you would have expected:      116 (98-135, 95%CI)

        The viability estimates that explain your observed counts is between:   0.900% - 6.500%

```

If your number of observed wells is at least 10% greater than the 95%CI upper limit if viability were 100%, then the program does not attempt to estimate viability, because viability would need to be greater than 100%, which is not possible. In such circumstances, the program will throw the following warning:

```
Your number of observed wells (5) is at least 10% greater than the 95%CI maximum if viability were 100% (2).
Consequently, it makes little sense to try and estimate viability >100%, so I will set min,max values to Inf.
I advise caution in interpreting such extreme value.
```

### Iterating between the two
So now you know that your viability of your culture is between 0.9% to 6.5%. If you were planning a second field campaign to isolate a taxon you really need, you can now use the `predict_wells` function to estimate how many wells you need to inoculate in order to guarantee (within 95% CI) that you isolate that taxon, by setting the viability to the minimum. For this case, we're going to be looking at a LOT of wells, so we can scale back the bootstrapping to `-b 999` to speed things up:

```
python viability_test/viability_test.py -b 999 predict_wells -w 2000 -i 2 -r 0.5 -v 0.009
```

```
    You simulated a DTE where you inoculated 2000 wells with an average of 2 cells per well.
    Your estimated relative abundance for your taxon of interest was:       50.000%
    Your estimated viability of your taxon of interest was:                 0.900%

    After 999 simulations:
    The median number of positive (not taxon-specific) wells was:           1271 (1229-1311, 95%CI)
    The median number of wells inoculated with one cell was:                543 (503-580, 95%CI)
    The median number of wells containing your taxon of interest was:       1263 (1221-1306, 95%CI)
    The median number of pure wells for your taxon of interest was:         6 (2-12, 95%CI)

```

So, if you inoculated 20,000 wells, you would observe 2-12 taxon specific positive wells of growth. That's doable. Now imagine though that your taxon of interest is only 5% of your population.

```
python viability_test/viability_test.py -b 999 predict_wells -w 50000 -i 2 -r 0.05 -v 0.009
```

```

    You simulated a DTE where you inoculated 50000 wells with an average of 2 cells per well.
    Your estimated relative abundance for your taxon of interest was:       5.000%
    Your estimated viability of your taxon of interest was:                 0.900%

    After 999 simulations:
    The median number of positive (not taxon-specific) wells was:           42527 (42371-42682, 95%CI)
    The median number of wells inoculated with one cell was:                13534 (13342-13730, 95%CI)
    The median number of wells containing your taxon of interest was:       4761 (4640-4888, 95%CI)
    The median number of pure wells for your taxon of interest was:         7 (2-12, 95%CI)
```

In that situation, we would need to inoculate 50,000 wells to have a good shot at culturing your target taxon. That would be a good time to start thinking about enrichment culturing or playing around with medium to improve viability!

### Bulk mode for estimating viability

If you have performed a series of experiments, you can place your parameters for estimating viability into a *tab-separated* file. The file must have the following column headers:

1. `inoculum` - the number of cells added to each well
2. `wells` - the number of wells inoculated
3. `num_observed` - the number of observed positive wells
4. `rel_abund` - the relative abundance of your taxon of interest

You can have any other columns in the file as long as it is parseable into a `pandas.DataFrame` using the `read_csv` function.

If our input file `test.tsv` looked like this:

```
ASV     Site    rel_abund       inoculum        wells   num_observed
1000    ARD     0.00075632      2       460     0
1000    ARD2c   0.000264656     2       460     0
5512    CJ2     0.003824898     2       460     5
```

we could simply run:
```
python viability_test/viability_test.py -b 9999 -p 16 bulk -i test.tsv -o test.out
```

The program will take each row in turn and evaluate viability just like in the standard command line version above. The difference is it also outputs the results to a `tsv` file specified in the command. The outputs will be in the same row order as the inputs:

```
ASV     Site    rel_abund       inoculum        wells   num_observed    pure_well_med   pure_well_95pc_low      pure_well_95pc_high     viability_95pc_low      viabiilty_95pc_high
1000    ARD     0.00075632      2       460     0       0.0     0.0     1.0     0.001   0.999
1000    ARD2c   0.000264656     2       460     0       0.0     0.0     0.0     0.001   0.999
5512    CJ2     0.0038248979999999998   2       460     5       0.0     0.0     2.0     inf     inf
```
