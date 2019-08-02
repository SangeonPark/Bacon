#!/bin/bash

#/eos/cms/store/cmst3/user/pharris/production/00/WJetsToLNu_TuneZ2Star_8TeV-madgraph-tarball_Bacon/
eosdir=$1
gen=$2
name=$3
mkdir $name
eosdir=/eos/cms/store/cmst3/user/pharris/production/00/$eosdir
for x in `/afs/cern.ch/project/eos/installation/0.3.15/bin/eos.select ls  $eosdir `; do
  echo $x
  runBacon 10000000 root://eoscms.cern.ch/$eosdir/$x $gen
  mv Output.root $name/$x
done
