#!/bin/sh

#should be in Lindquist-Neil-Work/hpcg
git checkout mixed-precision
git rebase --onto master mixed-precision-split
git branch -f mixed-precision-split master

#mixed-precision requires a forced push
git push origin +master +mixed-precision +mixed-precision-split
