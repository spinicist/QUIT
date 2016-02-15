# QUIT Changelog #

## Version 1.1 ##

1. Bug-fixes for qimcdespot.
2. All programs now print the version number with the usage string.
3. Improved multi-threading that will make more consistent use of available
threads.
4. Added 'qisplitsubject' for splitting up images with multiple subjects
scanned at the same time (e.g. ex-vivo scanning).

## Version 1.0 ##

1. I gave up and used a version number for simplicity. The changelog will
only list major changes and new tools, not every small change.
2. Faster qimcdespot - the matrix exponential step has been replaced with an
approximate expression that gives equivalent results but is much faster.
3. Correct echo-time signal equations for SPGR and SSFP. Use 'SPGR_ECHO' or
'SSFP_ECHO' in the input to qimcdespot or qisignal to use them.
4. Expanded qissfpbands with extra output options.
