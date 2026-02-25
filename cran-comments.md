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
* Document return values from all functions.
* Save `par` values before plotting, returning to original values afterwards.
* `app_2025()` saves configuration to a temporary file, not a user-space file.
* Trim an example in the documentation for `strike()` so it performs in under
  1s on my machine.

## Changes to address NOTES on remote builds

Remote builds (e.g. with `devtools::check_win_release()`) produced some NOTE
messages about spelling errors in the DESCRIPTION file. These are now enclosed
in backticks.  (They are also listed in the inst/WORDLIST file, but that does
not prevent the NOTE messages on the remote builds.)

## A note on spelling in the DESCRIPTION file

The following messages appear in win-builder tests, even though I have enclosed
those words in backwards tick marks in the DESCRIPTION and also included them
in `inst/WORDLIST`. I see no way to get around the NOTE messages, and I hope this is not a problem with CRAN submission.

```
Possibly misspelled words in DESCRIPTION:
  Lethality (10:958, 10:1082)
  lethality (10:276)
  whalestrike (10:1221)
```
