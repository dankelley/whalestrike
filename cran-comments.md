## R CMD check results

0 errors | 0 warnings | 1 note

* This is a revision of version 0.6.0 (not on CRAN), based on helpful comments
  from a CRAN reviewer. See the next section.

## Changes as suggested on previous submission to CRAN

The present version, 0.6.1, addresses comments on the 0.6.0 version that I
submitted a few days ago. I outline the changes briefly below. (At
https://github.com/dankelley/whalestrike/issues/33 I provide more detail, along
with the full comments that I received from the very helpful CRAN reviewer.)

* Add a full paragraph for Description, within the DESCRIPTION file.
* Change http-type reference to doi-type references, in the DESCRIPTION file.
* Update function documentation to state return values.
* Save `par` values before plotting, returning to original values afterwards.
* `app_2025()` saves configuration to a temporary file, not a user-space file.


