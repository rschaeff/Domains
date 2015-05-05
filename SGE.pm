package Domains::SGE;
require Exporter;

use warnings;
use strict;

use Data::Dumper;
use Carp;

our @ISA = qw(Exporter);
our @EXPORT = (
		"&qsub",
		"&qsub_fast",
		"&job_create", 
		"&qstat_wait_list", 
		"&qstat_wait",
		"&throttled_qsub",
		);


sub qsub_fast { 
	my $sub = 'qsub_fast';
	my $job_file = shift @_;
	system("qsub $job_file");
	#print `qsub  $job_file`;
	return 1;
}
sub qsub { 
	my $sub = 'qsub';

	my $job_file = shift @_;
	#if (!-f $job_file){ 
	#	die "ERROR! $sub: $job_file not found\n";
	#}

	my $output = `qsub -terse $job_file`; #Use -terse for better parsing.

	my $job_id;
	if ($output =~ /(\d+)/) { 
		$job_id = $1;
		return $job_id;
	}else{
		print "ERROR $sub: $output\n";
		return 0;
	}
}
sub throttled_qsub { 
	my $sub = 'throttled_qsub';

	#use Term::ProgressBar;
	my ($job_aref, $limit) = @_;

	#my $progress = Term::ProgressBar->new({count => scalar(@$job_aref), ETA => 'linear'}); 

	my @lines = split(/\n/, `qstat`);
	my $wc = scalar(@lines);
	my $i = 0;
	my $jobs = scalar(@$job_aref);
	JOB:
	while ($i <= scalar(@$job_aref)) { 

		$wc = `qstat | wc`;
		@lines = split(/\n/, `qstat`);
		$wc = scalar(@lines);
		if ($wc < $limit) { 
			my $j = 0;
			while ($j < 1000) { 
				if ($i < $jobs) { 
					print "qsub $i $$job_aref[$i]\n";
					qsub_fast($$job_aref[$i]);
				}else{
					last JOB;
				}
				#$progress->update($i);
				$i++;
				$j++;
			}
		}
	}

	while ($wc > 5) { 
		@lines = split(/\n/, `qstat`);
		$wc = scalar(@lines);
	}
	return ;
}
sub job_create { 
	my $sub = 'job_create';
	my ($fn, $lns) = @_;

	open (OUT, ">$fn") or die "ERROR! Could not open $fn for writing:$!\n";
		print OUT "#!/bin/bash\n";
		print OUT "#\$ -cwd\n";
		print OUT "#\$ -j y \n";
		print OUT "#\$ -S /bin/bash\n";
		print OUT "#\$ -M dustin.schaeffer\@gmail.com\n";

	if (ref $lns eq 'ARRAY') { 
		foreach my $ln (@$lns) { 
			print OUT "$ln\n";
		}
	}else{
		print OUT "$lns\n";
	}
	close OUT;
}
sub qstat_wait { 

	my $sub = 'qstat_wait';

	my ($qid) = @_;
	if (!$qid || $qid !~ /\d+/) { 
		die "ERROR! $sub: $qid is not a valid job id\n";
	}

	my $qstring = `qstat -j $qid 2> /dev/null`;
	if (!$qstring) { 
		return 0;
	}else{
		return 1;
	}
}
sub qstat_wait_list { 
	my $sub = 'qstat_wait_list';

	my ($qid_list_aref) = @_;

	my %qid_lookup;
	foreach my $qid (@$qid_list_aref) { 
		$qid_lookup{$qid}++;
	}

	my %not_complete;
	my $string = `qstat`;
	my @lines = split(/\n/, $string);
	
	foreach my $line (@lines) { 
		if ($line !~ /^\d+/) { next } 
		my @F = split(/\s+/, $line);
		#print "$line\n";
		my $qid = $F[0];
		#print "$qid\n";
		if (!$qid || $qid !~ /\d+/) { 
			#die "ERROR! $sub: $qid is not a valid job id\n";
			croak "ERROR! $sub: $qid is not a valid job id\n";
		}

		#my $qstring = `qstat -j $qid 2> /dev/null`;
		#print "DEBUG $sub: $qid $qstring\n";
		if ($qid_lookup{$qid}) { 
			print "$qid running..\n";
			$not_complete{$qid}++;
		}	
	}

	foreach my $qid (@$qid_list_aref) { 
		if ($not_complete{$qid}) { 
			print "NOT COMPLETE\n";
			return 1;
		}

	}
	print "COMPLETE!\n";
	return 0;
}
