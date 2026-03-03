## R CMD check results

0 errors | 0 warnings | 1 note

* This is a revision of version 0.6.1 (not on CRAN), based on helpful comments
  from a CRAN reviewer. Please see the next section for an overview of the
  changes.

## Changes as suggested on submitted version 0.6.1

1. The on.exit() function is called immediately after calls to par().

2. In DESCRIPTION, the use of backticks has been removed for words that are not
   package names. I had tried doing this as a way to prevent warnings about
   spelling (which occur even though the words are listed in inst/WORDLIST).
   For reference, the words are as follows.


