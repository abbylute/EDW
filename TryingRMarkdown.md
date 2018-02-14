Trying R Markdown
================
AL
2/14/2018

This is my first markdown document

``` r
tt <- read.table('Temperature_analysis/Tmax_meta.txt')
dim(tt)
```

    ## [1] 6725    3

``` r
head(tt)
```

    ##   elevin   latin     lonin
    ## 1    351 48.5000 -124.0000
    ## 2     12 48.3333 -123.6333
    ## 3     75 48.9333 -123.7500
    ## 4     37 48.5333 -123.3667
    ## 5     26 48.5167 -123.3667
    ## 6      1 48.7167 -123.5500
