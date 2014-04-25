#! /usr/bin/perl -w
use strict;

#
#	Run the wrsaMVGWA program on the Janelle's wilds RSA data
#

my $nuRval = shift @ARGV;
chomp $nuRval;

my $nuGval = shift @ARGV;
chomp $nuGval;

my $snpKind = shift @ARGV;
chomp $snpKind;

my $dayNam = shift @ARGV;
chomp $dayNam;

####
my $shScript = <<'END_SH';
#!/bin/sh
#PBS -q v4
#PBS -l nodes=1,walltime=5:00:00
#PBS -N wrsaMMGWArun_KIND_DAY_NUR_NUG_NUM
#PBS -A jgm45_0001
#PBS -j oe

echo "PBS_JOBID: $PBS_JOBID"
echo "PBS_JOBNAME: $PBS_JOBNAME"
echo "HOSTNAME: $HOSTNAME"

cd $TMPDIR

cp $HOME/HBqgen/wrsaMVGWA .
cp -v $HOME/sharedData/rootArchitecture/wilds/data/*WRSADAY.gbin .
cp $HOME/sharedData/rootArchitecture/wilds/data/wSNPD700KKIND.gbin .
cp $HOME/sharedData/rootArchitecture/wilds/data/wSNPD700KXtX.gbin .
cp $HOME/sharedData/rootArchitecture/wilds/data/wSNPD700Keval.gbin .
cp $HOME/sharedData/rootArchitecture/wilds/data/wSNPD700KRmat.gbin .

./wrsaMVGWA -b 10000 -s 20000 -t 10 -T 200 -n NUR -g NUG -c NUM -f DAY -F KIND -C 8 > trkDAYKIND_NUR_NUG_NUM.txt &

# Need to keep running until all background jobs are finished, otherwise
# the batch system will kill background jobs when this script exits.
while [ $(jobs -r | wc -l) -gt 0 ]; do
    sleep 120
    if [ -e trkDAYKIND_NUR_NUG_NUM.txt ]; then cp trkDAYKIND_NUR_NUG_NUM.txt $HOME/rootArchitecture/MVcpp/; fi
done

cp -v bet* $HOME/sharedData/rootArchitecture/wilds/saved_chains/
cp -v mhl* $HOME/sharedData/rootArchitecture/wilds/saved_chains/

if [ $(ls | grep -c 'Sig') -gt 0 ]; then cp -v Sig* $HOME/sharedData/rootArchitecture/wilds/saved_chains/; fi
if [ $(ls | grep -c 'LN') -gt 0 ]; then cp -v LN* $HOME/sharedData/rootArchitecture/wilds/saved_chains/; fi


rm *.gbin


END_SH
####


for my $chNum (1..5){
	my $shFile = 'wrsaGWA'.$snpKind.'_'.$dayNam.'_'.$chNum.'_'.$nuRval.'_'.$nuGval.'.sh';
	my $locSh = $shScript;
	$locSh =~ s/NUM/$chNum/g;
	$locSh =~ s/DAY/$dayNam/g;
	$locSh =~ s/KIND/$snpKind/g;
	$locSh =~ s/NUG/$nuGval/g;
	$locSh =~ s/NUR/$nuRval/g;
	open (SHOUT, ">$shFile");
	print SHOUT $locSh;
	close SHOUT;
	system("nsub $shFile");
}

