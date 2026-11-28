#!/usr/bin/env bash
# Build a greyscale double-JPEG-compression dataset from never-compressed
# colour PNGs, following the pipeline in dataset/README.txt.
in_dir=$1
out_dir=$2
qualities=(50 60 70 75 80 85 90 95)
mkdir -p "$out_dir"
tmp=$(mktemp -d); trap 'rm -rf "$tmp"' EXIT
for f in "$in_dir"/train*.png; do
  name=$(basename "$f" .png)
  pngtopnm "$f" | ppmtopgm > "$tmp/grey.pgm"
  for q1 in "${qualities[@]}"; do
    cjpeg -quality "$q1" "$tmp/grey.pgm" | djpeg > "$tmp/once.pgm"
    for q2 in "${qualities[@]}"; do
      cjpeg -quality "$q2" "$tmp/once.pgm" | djpeg | pnmtopng > "$out_dir/${name}_q$(printf %02d "$q1")_q$(printf %02d "$q2").png"
    done
  done
  echo "done: $name"
done
