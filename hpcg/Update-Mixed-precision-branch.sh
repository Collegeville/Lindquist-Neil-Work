#!/bin/sh

#should be in Lindquist-Neil-Work/hpcg
git checkout mixed-precision
git rebase --onto master mixed-precision-split
git branch -f mixed-precision-split master

#mixed-precision requires a forced push
git push origin +mixed-precision
#mixed-precision-split usually shouldn't need a forced push, so only
#do so if nessaccery
git push origin mixed-precision-split
