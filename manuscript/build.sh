#!/usr/bin/env bash
set -euo pipefail
cd -- "$(dirname -- "$0")"
documents=(methodology-8.3.5 methodology-8.3.5-supplement)
for document in "${documents[@]}"; do
  pdflatex -interaction=nonstopmode -halt-on-error "$document.tex"
  bibtex "$document"
done
for pass in 1 2 3; do
  for document in "${documents[@]}"; do
    pdflatex -interaction=nonstopmode -halt-on-error "$document.tex"
  done
done
