rps11 example output

Read 2 tree file names from list file rps11.txt

        Info        Trees       Description
    65.14748         1000       trees3prime.t
    70.52154         1000       trees5prime.t
    67.83451         2000       average
    46.17921         2000       merged
    21.65530         2000       difference

Required 0.082 total seconds

----- 5-prime info -----------------------------------------------------------------------
     1 -*---  1000    1.00000    0.06497     1  1000     1
     2 -****  1000    1.00000    0.03122     1  1000     1
     3 --***  1000    1.00000    0.13975     1  1000     1
     4 ---*-  1000    1.00000    0.01351     1  1000     1
     5 ----*  1000    1.00000    0.01120     1  1000     1
     6 --*--  1000    1.00000    0.00565     1  1000     1
     7 --*-*   713    0.71300    0.01166     1  1000   202
     8 --**-   147    0.14700    0.00695    10   996   122
     9 ---**   140    0.14000    0.00539     3   994   124
  Tree majrule = [&U] (1:0.03122,2:0.06497,(4:0.01351,(5:0.01120,3:0.00565)0.71300:0.01166)1.00000:0.13975)1.00000;

Lindley Information
  number of taxa: 4
  total prior entropy: 2.70805

Clade --***
  H          = 0.79829
  Hp         = 1.09861
  Hp - H     = 0.30032
  w          = 1.00000
  I          = 0.30032
  I (% max.) = 11.08993

Clade -****
  H          = 0.00000
  Hp         = 1.60944
  Hp - H     = 1.60944
  w          = 1.00000
  I          = 1.60944
  I (% max.) = 59.43161

total I          = 1.90976
total I (% max.) = 70.52154

Lindley information based on conditional clade distribution:
  Posterior entropy:   0.79829
  Prior entropy:       2.70805
  Information gain:    1.90976
  % max. information:  70.52154

Larget coverage:
  Distinct tree topologies: 3
  Posterior coverage:       1.00000

Lindley information based on observed tree topology distribution:
  Posterior entropy:   0.79829
  Prior entropy:       2.70805
  Information gain:    1.90976
  % max. information:  70.52154

----- 3-prime info -----------------------------------------------------------------------
     1 ----*  1000    1.00000    0.00860     1  1000     1
     2 ---*-  1000    1.00000    0.00949     1  1000     1
     3 -****  1000    1.00000    0.01974     1  1000     1
     4 -*---  1000    1.00000    0.07647     1  1000     1
     5 ---**  1000    1.00000    0.07326     1  1000     1
     6 --*--  1000    1.00000    0.01525     1  1000     1
     7 -*-**   602    0.60200    0.00992     1  1000   231
     8 -**--   228    0.22800    0.00501     4   996   172
     9 --***   170    0.17000    0.00493     5   993   146
  Tree majrule = [&U] (1:0.01974,3:0.01525,(2:0.07647,(5:0.00860,4:0.00949)1.00000:0.07326)0.60200:0.00992)1.00000;

Clade -*-**
  H          = 0.00000
  Hp         = 1.09861
  Hp - H     = 1.09861
  w          = 0.60200
  I          = 0.66136
  I (% max.) = 24.42217

Clade --***
  H          = 0.00000
  Hp         = 1.09861
  Hp - H     = 1.09861
  w          = 0.17000
  I          = 0.18676
  I (% max.) = 6.89663

Clade -****
  H          = 0.94382
  Hp         = 1.85992
  Hp - H     = 0.91610
  w          = 1.00000
  I          = 0.91610
  I (% max.) = 33.82868

total I          = 1.76423
total I (% max.) = 65.14748

Lindley information based on conditional clade distribution:
  Posterior entropy:   0.94382
  Prior entropy:       2.70805
  Information gain:    1.76423
  % max. information:  65.14748

Larget coverage:
  Distinct tree topologies: 3
  Posterior coverage:       1.00000

Lindley information based on observed tree topology distribution:
  Posterior entropy:   0.94382
  Prior entropy:       2.70805
  Information gain:    1.76423
  % max. information:  65.14748

----- merged info ------------------------------------------------------------------------

Lindley information based on conditional clade distribution:
  Posterior entropy:   1.45749
  Prior entropy:       2.70805
  Information gain:    1.25056
  % max. information:  46.17921

Larget coverage:
  Distinct tree topologies: 5
  Posterior coverage:       1.00000

Lindley information based on observed tree topology distribution:
  Posterior entropy:   1.45749
  Prior entropy:       2.70805
  Information gain:    1.25056
  % max. information:  46.17921
