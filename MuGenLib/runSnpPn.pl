#! /usr/bin/perl -w
use strict;

#
#	Run the snpP program posterior samples
#

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

my $nlSw = shift @ARGV;
chomp $nlSw;

####
my $shScript = <<'END_SH';
#!/bin/sh
#PBS -q v4
#PBS -l nodes=1,walltime=10:00:00
#PBS -N snpPnRunHER_PDRKINDNLM
#PBS -A jgm45_0001
#PBS -j oe

echo "PBS_JOBID: $PBS_JOBID"
echo "PBS_JOBNAME: $PBS_JOBNAME"
echo "HOSTNAME: $HOSTNAME"
echo "TMPDIR: $TMPDIR"

if [ ! -d $TMPDIR ]; then mkdir $TMPDIR; fi

cd $TMPDIR

cp $HOME/HBqgen/snpPn .
cp -v $HOME/sharedData/PDR/saved_chainsHER/betNLMSNPKIND*.gbin .
cp -v $HOME/sharedData/PDR/saved_chainsHER/plsNLMSNPKIND*.gbin .

./snpPn -l MNM -c NCN -s NSNP -d DNM -t 8 -f KIND NLF > trk_PDRKIND_NLMHERoutN.txt &

# Need to keep running until all background jobs are finished, otherwise
# the batch system will kill background jobs when this script exits.
while [ $(jobs -r | wc -l) -gt 0 ]; do
    sleep 120
    if [ -e trk_PDRKIND_NLMHERoutN.txt ]; then cp trk_PDRKIND_NLMHERoutN.txt $HOME/PDR/; fi
done

cp -v pValNLMSNPKIND*.gbin $HOME/sharedData/PDR/saved_chainsHER/

rm betNLMSNPKIND*.gbin
rm plsNLMSNPKIND*.gbin
rm trk_PDRKIND_NLMHERoutN.txt


END_SH
####

my $nlm = '';
$nlm = 'Nl' if ($nlSw eq 'n');
my $nlf = '';
$nlf = '-N' if ($nlSw eq 'n');

$her = '' if ($her eq 'N');

my $shFile = 'snpPnRun'.$her.'_'.$cDir.$flKind.$nlm.'.sh';
my $locSh = $shScript;

$locSh =~ s/HER/$her/g;
$locSh =~ s/NCN/$nChn/g;
$locSh =~ s/MNM/$nMcmc/g;
$locSh =~ s/NSNP/$snpNum/g;
$locSh =~ s/DNM/$dPhen/g;
$locSh =~ s/PDR/$cDir/g;
$locSh =~ s/KIND/$flKind/g;
$locSh =~ s/NLM/$nlm/g;
$locSh =~ s/NLF/$nlf/g;

open (SHOUT, ">$shFile");
print SHOUT $locSh;
close SHOUT;
system("nsub $shFile");

