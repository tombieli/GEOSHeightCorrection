#!/bin/bash
docker run --rm -v `pwd`:/data pandoc/core --from=markdown --citeproc --bibliography bibliography.bib --csl https://www.zotero.org/styles/ieee --metadata link-citations=true --to=gfm --output=README.md README.src.md
mv README.md ..
