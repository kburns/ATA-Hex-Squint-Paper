#! /bin/bash

if [ x"$1" = x-d ] ; then
    doit=
    shift
else
    doit=:
fi

if [ x"$1" = x-r ] ; then
    resid=true
    shift
else
    resid=false
fi

lower=$1
upper=$2

if [ x"$upper" = x ] ; then
    echo >&2 "error: specify lower and upper numbers"
    exit 1
fi

h=$HOME/ATA-Hex-Squint-Paper/sqsim
v=../tiny.uv
cf=$HOME/sw/carmafiller/carmafiller
imgargs="cellsize=10arcsec wprojplanes=128 npix=2048 operation=csclean niter=1000"
cfargs="polmode=1 snumbase=1"

for i in $(seq $lower $upper) ; do
    i=$(printf '%03d' $i)
    echo $i
    $h/sqsim.py $v $i.uv >$i.st 2>$i.log || exit 1
    $cf vis=$i.uv ms=$i.ms $cfargs &>> $i.log || exit 1
    mfsimager $i.ms image=$i $imgargs &>> $i.log || exit 1
    if $resid ; then
	$h/sqfit.py $i.st $i.restored $i.fr.fits >$i.fk || exit 1
    else
	$h/sqfit.py $i.st $i.restored >$i.fk || exit 1
    fi
    #$doit rm -rf $i.uv $i.st $i.ms $i.model $i.residual $i.restored $i.log
    $doit rm -rf $i.uv $i.ms $i.model $i.residual $i.log
done
