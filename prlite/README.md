# prlite description

A few support functions common to patrec, emdw and adcv

## Dependencies
 `sudo apt-get install libssl-dev`

## Installing prlite

Since prlite depends on `gLinear` you need to install it first. Follow
the instructions in its `README.md` file. We normally create a
`~/devel` directory and install both `gLinear` and `prlite` there.

## To compile

Your prlite directory should have a `build` subdirectory in it. If it does
not, create one. Then:

	`cd build`
	`cmake ../`
	`make prlite`

This should create `libprlite.so`
