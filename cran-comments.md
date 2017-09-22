---
title: "Untitled"
author: "Justin L"
date: "9/20/2017"
output: html_document
---


##Test environments
* local OS X El Capitan, R 3.4.1


## R CMD check results
There was no Erros.

There was 4 WARNINGs:

* checking S3 generic/method consistency
  See section ‘Generic functions and methods’ in the ‘Writing R   Extensions’ manual

* checking Rd cross-references 
  Missing link or links in documentation object 'dr4pl.formula.Rd':
  ‘gof.dr4pl’

  See section 'Cross-references' in the 'Writing R Extensions' manual.

* checking for missing documentation entries
  Undocumented data sets:
  ‘Dirk_data_14’ ‘Dirk_data_15’ ‘Dirk_data_16’ ‘Dirk_data_17’
  ‘Dirk_data_18’ ‘Dirk_data_19’ ‘Dirk_data_20’ ‘Dirk_data_21’
  ‘Dirk_data_22’ ‘Dirk_data_23’ ‘Dirk_data_24’ ‘Dirk_data_25’
  ‘Dirk_data_26’ ‘Dirk_data_27’ ‘Dirk_data_28’ ‘Dirk_data_29’
  ‘Dirk_data_30’ ‘Dirk_data_31’
  All user-level objects in a package should have documentation entries.
  See chapter ‘Writing R documentation files’ in the ‘Writing R
  Extensions’ manual.

  These data sets are NOT intended to be included in the final package

* checking Rd \usage sections 
  Undocumented arguments in documentation object 'GradientFunction'
    ‘x’ ‘y’
  Documented arguments not in \usage in documentation object 'GradientFunction':
   ‘dose’ ‘response’

  Undocumented arguments in documentation object 'Hessian'
    ‘y’

  Undocumented arguments in documentation object 'dr4pl'
   ‘constrained’ ‘grad’ ‘init.parm’ ‘method.init’ ‘method.optim’
   ‘method.robust’ ‘trace’

  Functions with \usage entries need to have the appropriate \alias
  entries, and all their arguments documented.
  The \usage entries must correspond to syntactically valid R code.
  See chapter ‘Writing R documentation files’ in the ‘Writing R
  Extensions’ manual.

There was 1 NOTE:


* checking R code for possible problems 
  print.summary.dr4pl: warning in printCoefmat(object$coefficients,
    P.value = TRUE, has.Pvalue = TRUE): partial argument match of
    'P.value' to 'P.values'
    
## Downstream dependencies
I have also run revdep_check() for a R CMD check on downstream dependencies. There are currently no known ERRORs or Warnings.