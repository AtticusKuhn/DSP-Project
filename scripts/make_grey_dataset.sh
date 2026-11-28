#!/usr/bin/env bash
# Build a greyscale experimentation dataset from never-compressed colour PNGs,
# following the pipeline in dataset/README.txt.
in_dir=$1; out_dir=$2

qualities=(50 60 70 75 80 85 90 95)

mkdir -p "$out_dir"

for f in "$in_dir"/train*.png; do
  name=$(basename "$f" .png)
  grey="$out_dir/${name}_grey.png"
  pngtopnm "$f" | ppmtopgm | pnmtopng > "$grey"
  for q in "${qualities[@]}"; do
    pngtopnm "$grey" | cjpeg -quality "$q" | djpeg | pnmtopng > "$out_dir/${name}_q$(printf %02d "$q").png"
  done
  echo "done: $name"
done
