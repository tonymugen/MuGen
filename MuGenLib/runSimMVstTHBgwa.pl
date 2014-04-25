#! /usr/bin/perl -w
use strict;

#
#	Run the MVstTHBgwa program on sets of simulated data
#

my $simBeg = shift @ARGV;
chomp $simBeg;
my $simEnd = shift @ARGV;
chomp $simEnd;

die "simEnd has to be larger than simBeg\n" if ($simBeg > $simEnd);

my $nuRval = shift @ARGV;
chomp $nuRval;

my $nuGval = shift @ARGV;
chomp $nuGval;

my $herVal = shift @ARGV;
chomp $herVal;

####
my $shScript = <<'END_SH';
#!/bin/sh
#PBS -q v4
#PBS -l nodes=1,walltime=5:00:00
#PBS -N MMGWArunSMN_NUM_NUR_NUGHER
#PBS -A jgm45_0001
#PBS -j oe

echo "PBS_JOBID: $PBS_JOBID"
echo "PBS_JOBNAME: $PBS_JOBNAME"
echo "HOSTNAME: $HOSTNAME"
echo "TMPDIR: $TMPDIR"

if [ ! -d $TMPDIR ]; then mkdir $TMPDIR; fi

cd $TMPDIR

cp $HOME/HBqgen/MVstTHBgwa .
cp $HOME/HBqgen/KeigVal.gbin .
cp $HOME/HBqgen/LNind.gbin .
cp $HOME/HBqgen/R.gbin .
cp $HOME/sharedData/GWAS_sim/dataGSLHER/YoutHERSMN.gbin .
cp $HOME/sharedData/GWAS_sim/additional_data/MESA_XtX.gbin .
cp $HOME/sharedData/GWAS_sim/additional_data/MESAsnpTR.gbin .

./MVstTHBgwa -b 2000 -s 20000 -t 10 -T 200 -c NUM -C 8 -n NUR -g NUG -S 669958 -f outHERSMN > trkHERSMNoutNUR_NUG_NUM.txt &

# Need to keep running until all background jobs are finished, otherwise
# the batch system will kill background jobs when this script exits.
while [ $(jobs -r | wc -l) -gt 0 ]; do
    sleep 120
    if [ -e trkHERSMNoutNUR_NUG_NUM.txt ]; then cp trkHERSMNoutNUR_NUG_NUM.txt $HOME/GWAS_sim/; fi
done

cp -v Sig* $HOME/sharedData/GWAS_sim/saved_chainsHER/
cp -v bet* $HOME/sharedData/GWAS_sim/saved_chainsHER/
cp -v mhl* $HOME/sharedData/GWAS_sim/saved_chainsHER/

rm *.gbin


END_SH
####


for my $sim ($simBeg .. $simEnd){
	my $nJobs = `qstat | grep -c ajg67`;
	chomp $nJobs;
	
	while ($nJobs >= 40){
		sleep 300;
		$nJobs = `qstat | grep -c ajg67`;
		chomp $nJobs;
	}

	for my $chNum (1..5){
		my $shFile = 'MMGWA'.$sim.'_'.$chNum.'_'.$nuRval.'_'.$nuGval.$herVal.'.sh';
		my $locSh = $shScript;
		$locSh =~ s/NUM/$chNum/g;
		$locSh =~ s/HER/$herVal/g;
		$locSh =~ s/SMN/$sim/g;
		$locSh =~ s/NUR/$nuRval/g;
		$locSh =~ s/NUG/$nuGval/g;
		open (SHOUT, ">$shFile");
		print SHOUT $locSh;
		close SHOUT;
		system("nsub $shFile");
	}
}

