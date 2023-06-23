NEWS-robust
================

<!-- NEWS.md is generated from NEWS.Rmd. Please edit that file -->

    #> Last Update: 2023-06-23 23:17:57.700875

# robust v0.7-2 (2023-06-22)

- tests/robust-Ex.Rout.save: Follow up to changes in robustbase 0.99-0
- src/tmlfor.f: DFLOAT replaced by DBLE

# robust v0.7-1 (2022-07-08)

- Sfloat and Sint are now deprecated, see mail of Prof Riplay from
  08.07.2022 typedef-s inserted in robust.h

# robust v0.7-0 (2022-02-01)

- Replace the call to the deprecated function covMest() by a call to
  CovMest() (package ‘rrcov’), which returns an S4 object

- Code using S.h is converted to use R.h

- BLAS - Fix for USE_FC_LEN_T becoming the default in 4.2.0 (early
  notification)

# robust v0.6-1 (2021-11-16)

- Fixed error on Solaris

# robust v0.6-0 (2021-10-24)

- Fixed an error due to not explicitly including R_ext/Error.h - include
  this file implicitly in compatibility.c.
- Fixed: everything related to ‘covfm’ was moved to package fit.models
  (actually I only commented the exports in NAMESPACE and a couple of
  other references. All the code remains duplicated in robust, also the
  RD files)

# robust v0.5-0 (2020-03-08)

- New submission, maintainer Valentin Todorov
