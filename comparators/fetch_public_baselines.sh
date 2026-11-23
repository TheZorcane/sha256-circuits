#!/bin/sh
# Download the public SHA-256 circuits used as comparators (network access needed; nothing is executed).
# usage: sh fetch_public_baselines.sh DIR      then: python3 score_public.py DIR
set -e
D=${1:-public_raw}; mkdir -p "$D"
curl -fsSL -o "$D/sha256-input512.slp"        https://raw.githubusercontent.com/usnistgov/Circuits/master/data/slp/sha/sha256-input512.slp
curl -fsSL -o "$D/sha256-iv256+input512.slp"  "https://raw.githubusercontent.com/usnistgov/Circuits/master/data/slp/sha/sha256-iv256%2Binput512.slp"
curl -fsSL -o "$D/sha256Unpadded.txt"         https://raw.githubusercontent.com/matteocam/sha256-circuit/master/sha256Unpadded.txt
curl -fsSL -o "$D/sha256.txt"                 https://nigelsmart.github.io/MPC-Circuits/sha256.txt
ls -l "$D"
