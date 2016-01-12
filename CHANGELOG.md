# QUIT Changelog #

## Version 1.0 ##

1. I gave up and used a version number for simplicity. The changelog will
only list major changes and new tools, not every small change.
2. Faster qimcdespot - the matrix exponential step has been replaced with an
approximate expression that gives equivalent results but is much faster.
3. Correct echo-time signal equations for SPGR and SSFP. Use 'SPGR_ECHO' or
'SSFP_ECHO' in the input to qimcdespot or qisignal to use them.
4. Expanded qissfpbands with extra output options.