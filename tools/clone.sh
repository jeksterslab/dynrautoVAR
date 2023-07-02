#!/bin/bash

git clone git@github.com:jeksterslab/dynrautoVAR.git
rm -rf "$PWD.git"
mv dynrautoVAR/.git "$PWD"
rm -rf dynrautoVAR
