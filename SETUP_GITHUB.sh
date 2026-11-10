#!/usr/bin/env bash
# ============================================================
# Initialise and push the public repository.
#
# Run this from inside caribou-on-hormones/ AFTER:
#   1. Rscript ../13_build_public_dataset.R      (creates the data file)
#   2. cp ../caribou_outputs/public/* data/
#   3. Creating an EMPTY repo named caribou-on-hormones on GitHub
#      (no README, no .gitignore — this folder supplies both)
#
# This is a NEW repository on purpose. The previous repo
# (mammferraz/caribou-on-fcm-ft3) has collection coordinates committed in
# its history: deleting the file there would not remove them, because any
# earlier commit can still be checked out. Archive or delete that repo
# once this one is live.
# ============================================================
set -euo pipefail

REPO_URL="https://github.com/mammferraz/caribou-on-hormones.git"

# ---- Pre-flight: refuse to publish anything with coordinates -----------
# NOTE: these checks run against what git would actually TRACK, not what
# happens to sit on disk. Files matched by .gitignore are never pushed, so
# a leftover outputs/models/ directory on your machine is harmless as long
# as it stays ignored. The checks below verify exactly that.
echo "Pre-flight coordinate check..."

if [ ! -f data/caribou_hormones_analysis_ready.csv ]; then
  echo "  FAIL: data/caribou_hormones_analysis_ready.csv is missing."
  echo "        Run:  Rscript ../13_build_public_dataset.R"
  echo "              cp ../caribou_outputs/public/* data/"
  exit 1
fi

git init -q -b main 2>/dev/null || true
# Everything git would stage, honouring .gitignore
TRACKED=$(git add -An --ignore-missing . 2>/dev/null | sed "s/^add '//;s/'$//")

fail=0

# 1. No coordinate columns in any CSV that would be committed
while IFS= read -r f; do
  case "$f" in *.csv)
    if head -1 "$f" 2>/dev/null | grep -qiE '(^|,)[[:space:]]*(lat|latitude|lon|long|longitude|x_proj|y_proj)[[:space:]]*(,|$)'; then
      echo "  FAIL: $f has a coordinate column."; fail=1
    fi
    # decimal-degree longitudes for this study area (-80 to -95)
    if grep -qE '(^|,)-(8[0-9]|9[0-5])\.[0-9]{4,}' "$f" 2>/dev/null; then
      echo "  FAIL: $f contains decimal-degree longitude values."; fail=1
    fi ;;
  esac
done <<< "$TRACKED"

# 2. No serialised model objects. gam() stores its model frame inside the
#    fitted object, so s(x_proj, y_proj) models carry every coordinate.
if echo "$TRACKED" | grep -qi '\.rds$'; then
  echo "  FAIL: .rds files would be committed. mgcv objects embed the"
  echo "        coordinate model frame. Add them to .gitignore or delete."
  fail=1
fi

# 3. No spatial layers or rasters
if echo "$TRACKED" | grep -qiE '\.(geojson|gpkg|shp|shx|dbf|prj|tif|tiff)$'; then
  echo "  FAIL: spatial layers would be committed — these geolocate the animals."
  fail=1
fi

# 4. The raw dataset must never appear
if echo "$TRACKED" | grep -q 'All_data_updated.csv'; then
  echo "  FAIL: All_data_updated.csv would be committed — raw file WITH coordinates."
  fail=1
fi

[ "$fail" -eq 0 ] || { echo; echo "Aborted. Nothing initialised."; exit 1; }

# Report anything present-but-ignored, so it is a conscious decision
IGNORED=$(git status --porcelain --ignored 2>/dev/null | grep '^!!' | sed 's/^!! //')
if [ -n "$IGNORED" ]; then
  echo "  Present on disk but correctly IGNORED (will not be pushed):"
  echo "$IGNORED" | sed 's/^/    /'
fi

echo "  All checks passed."
echo

# ---- Initialise and push ----------------------------------------------
git add .

echo
echo "Files staged for the FIRST public commit:"
git status --short
echo
read -r -p "Review the list above. Proceed with commit and push? [y/N] " ok
[ "$ok" = "y" ] || { echo "Aborted. Nothing pushed."; exit 1; }

git commit -m "Initial public release: analysis code and de-identified dataset

Faecal cortisol and triiodothyronine metabolite analysis for boreal
woodland caribou, northern Ontario (2021-2024, n = 1,380).

Collection coordinates are withheld throughout: boreal caribou are listed
as Threatened under Canada's Species at Risk Act. Spatially derived
covariates are published in their place so the reported models remain
reproducible. See README.md."

git remote add origin "$REPO_URL"
git push -u origin main

echo
echo "Done. Next steps:"
echo "  1. Archive or DELETE mammferraz/caribou-on-fcm-ft3. Its history still"
echo "     contains the collection coordinates, so archiving leaves them"
echo "     publicly readable — deleting is the safer choice."
echo "  2. Keep this repo PRIVATE until the paper is accepted, then flip it"
echo "     to public: Settings > General > Danger Zone > Change visibility."
echo "     NOTE: flipping to public exposes the ENTIRE commit history, not"
echo "     just the current files. Private is not a staging area."
