# Bounds for each age according to the stratigraphic constraints

This function gives the bounds of each age according to the reduced
network (better computational time) it will return a list of list each
list within have two vector \*index of the ages that preceed the
selected age (`lower`) \*index of the ages that comes after the selected
age (\`upper“) The index are shift by one to includ the study period for
the younger and older :

## Usage

``` r
findbounds(network)
```
