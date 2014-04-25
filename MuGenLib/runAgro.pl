#! /usr/bin/perl -w
use strict;

#
#	Run the agroMVGWA program on the Stuttgart Agronomic data
#

my $nuValRep = shift @ARGV;
chomp $nuValRep;

my $nuValG = shift @ARGV;
chomp $nuValG;

my $snpNam = shift @ARGV;
chomp $snpNam;

my $snpShrt = $snpNam;
$snpShrt =~ s/(tr)|(p(\d{1,2}))//g;
####
my $shScript = <<'END_SH';
#!/bin/sh
#PBS -q v4
#PBS -l nodes=1,walltime=15:00:00
#PBS -N agroMMGWArun_KIND_NUR_NUG_NUM
#PBS -A jgm45_0001
#PBS -j oe

echo "PBS_JOBID: $PBS_JOBID"
echo "PBS_JOBNAME: $PBS_JOBNAME"
echo "HOSTNAME: $HOSTNAME"

if [ ! -d $TMPDIR ]; then mkdir $TMPDIR; fi
cd $TMPDIR

cp $HOME/HBqgen/agroMVGWA .
cp -v $HOME/sharedData/agro/data/agro*.gbin .
cp -v $HOME/sharedData/agro/data/SNPKIND.gbin .
cp -v $HOME/sharedData/agro/data/SNPSMKNDXtX.gbin .

./agroMVGWA -b 50000 -s 40000 -t 20 -T 400 -n NUR -g NUG -c NUM -f KIND -C 8 > trkKINDoutNUR_NUG_NUM.txt &

# Need to keep running until all background jobs are finished, otherwise
# the batch system will kill background jobs when this script exits.
while [ $(jobs -r | wc -l) -gt 0 ]; do
    sleep 120
    if [ -e trkKINDoutNUR_NUG_NUM.txt ]; then cp trkKINDoutNUR_NUG_NUM.txt $HOME/agro/; fi
done

cp -v Sig* $HOME/sharedData/agro/saved_chains/
cp -v LN* $HOME/sharedData/agro/saved_chains/
cp -v bet* $HOME/sharedData/agro/saved_chains/
cp -v pls* $HOME/sharedData/agro/saved_chains/

rm *.gbin


END_SH
####

for my $chNum (1..5){
	my $shFile = 'agroGWA'.$snpNam.'_'.$chNum.'_'.$nuValRep.'_'.$nuValG.'.sh';
	my $locSh = $shScript;
	$locSh =~ s/NUM/$chNum/g;
	$locSh =~ s/KIND/$snpNam/g;
	$locSh =~ s/SMKND/$snpShrt/g;
	$locSh =~ s/NUR/$nuValRep/g;
	$locSh =~ s/NUG/$nuValG/g;
	open (SHOUT, ">$shFile");
	print SHOUT $locSh;
	close SHOUT;
	system("nsub $shFile");
}

