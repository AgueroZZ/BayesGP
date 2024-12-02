## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new release.
* There is no error or warning during the CMD check. The note is due to the size of the complied cpp file.

## Resubmission
This is a resubmission. In this version I have:

* Modified the description field in DESCRIPTION file based on the requirements.
* Used `par(no.readonly = TRUE)` to make sure user settings will not be altered.
* Fixed the missing Rd-tags from the last submission.