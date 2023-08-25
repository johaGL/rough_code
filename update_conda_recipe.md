I found in https://bioconda.github.io/contributor/guidelines.html?highlight=bioconductor#hashes  how to get the "shasum": 
`wget -O- $URL | shasum -a 256`
so for our release 0.1.1 , I ran in the terminal: 
`wget -O- https://pypi.io/packages/source/d/dimet/dimet-0.1.1.tar.gz | shasum -a 256`
And I've got : bbff4a390d5e2f3d55410e46155025b688a1bf0da79361d1b688dbd2e82ec487



