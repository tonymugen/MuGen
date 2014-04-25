#! /usr/bin/perl -w
use strict;

#
#	Run the snpP program posterior samples
#

my $nuVal = shift @ARGV;
chomp $nuVal;

my $flKind = shift @ARGV; # flKind should not have a trailing underscore, but should have the leading one if need be
chomp $flKind;

my $her = shift @ARGV;
chomp $her;

my $snpNum = shift @ARGV;
chomp $snpNum;

my $dPhen = shift @ARGV;
chomp $dPhen;

my $nChn = shift @ARGV;
chomp $nChn;

my $nMcmc = shift @ARGV;
chomp $nMcmc;

my $cDir = shift @ARGV;
chomp $cDir;

####
my $shScript = <<'END_SH';
#!/bin/sh
#PBS -q v4
#PBS -l nodes=1,walltime=15:00:00
#PBS -N snpPrun_DIRKIND_NUV
#PBS -A jgm45_0001
#PBS -j oe

echo "PBS_JOBID: $PBS_JOBID"
echo "PBS_JOBNAME: $PBS_JOBNAME"
echo "HOSTNAME: $HOSTNAME"

cd $TMPDIR

cp $HOME/HBqgen/snpP .
cp -v $HOME/sharedData/DIR/saved_chainsHER/betSNPKIND*.gbin .
cp -v $HOME/sharedData/DIR/saved_chainsHER/mhlSNPKIND*.gbin .

./snpP -n NUV -l MNM -c NCN -s NSNP -d DNM -t 8 -f KIND > trkDIRKINDoutNUV.txt &

# Need to keep running until all background jobs are finished, otherwise
# the batch system will kill background jobs when this script exits.
while [ $(jobs -r | wc -l) -gt 0 ]; do
    sleep 120
    if [ -e trkDIRKINDoutNUV.txt ]; then cp trkDIRKINDoutNUV.txt $HOME/DIR/; fi
done

cp -v pVal* $HOME/sharedData/DIR/saved_chainsHER/

rm *.gbin


END_SH
####

my $shFile = 'snpPrun'.'_'.$cDir.$flKind.'_'.$nuVal.'.sh';
my $locSh = $shScript;
$her = '' if ($her eq 'N');

$locSh =~ s/HER/$her/g;
$locSh =~ s/NCN/$nChn/g;
$locSh =~ s/MNM/$nMcmc/g;
$locSh =~ s/NSNP/$snpNum/g;
$locSh =~ s/DNM/$dPhen/g;
$locSh =~ s/DIR/$cDir/g;
$locSh =~ s/KIND/$flKind/g;
$locSh =~ s/NUV/$nuVal/g;
open (SHOUT, ">$shFile");
print SHOUT $locSh;
close SHOUT;
system("nsub $shFile");

