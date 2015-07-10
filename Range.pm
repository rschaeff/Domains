package Domains::Range;
require Exporter;

use warnings;
use strict;

use Data::Dumper;

use Carp;

our @ISA = qw(Exporter);
our @EXPORT = (	"&multi_chain_pdb_range_expand",
		"&pdb_range_expand",
		"&multi_chain_pdb_rangify",
		"&pdb_rangify",
		"&rangify",
		"&multi_chain_rangify",
		"&multi_chain_range_expand",
		"&multi_chain_range_include",
		"&multi_chain_range_exclude",
		"&range_expand", 
		"&range_include", 
		"&range_exclude",
		"&region_boundary_byres",
		"&struct_region",
		"&multi_chain_struct_region",
		"&scopify_range",
		"&scop_range_split",
		"&multi_chain_region_coverage",
		"&multi_chain_region_coverage2",
		"&region_coverage",
		"&multi_chain_residue_coverage",
		"&residue_coverage",
		"&union_range_boundary",
		"&multi_chain_ungap_range",
		"&ungap_range",
		"&ungap_range_aref"
		);

my $DEBUG = 0;

#return the PDB seq_id for a set of residues given in 'normalized structure' reference
sub struct_region { 	
	my $sub = 'struct_region';
	my ($aref, $struct_aref) = @_;
	if (ref($aref) ne "ARRAY" || ref($struct_aref) ne "ARRAY") { 
		print "WARNING! $sub: Bad inputs\n";
		return 0;
	}
	my @return;
	my $s1 = scalar(@$aref);
	my $s2 = scalar(@$struct_aref);
	if ($DEBUG > 1) { 
		#printf "DEBUG s1:%i s2:%i\n", $s1, $s2;
		print "DEBUG: Struct_reference ID ". rangify(@$aref) . "\n";
		print "DEBUG: Struct_seq_reference ID " . rangify(@$struct_aref) ."\n";
	}
	my $warning = 0;
	for (my $i = 0; $i < scalar(@$aref); $i++) {
		if (defined($$struct_aref[$$aref[$i]-1])) { 
			push(@return, $$struct_aref[$$aref[$i]-1]) ;
		}else{
			#print Dumper($struct_aref);
			#print "WARNING: struct_region: residue not defined or not structured-> $i $$aref[$i]\n";
			$warning++;
		}
	}  
	if ($warning > 20) { 
		print "WARNING! $sub: $warning unstructured residues for " . scalar(@$aref) . " total residues\n";
	}
	if ($DEBUG > 1 ) { 
		print "DEBUG: return " . rangify(@return) . "\n";
	}
	return \@return;
}

sub multi_chain_struct_region { 
	my $sub = 'multi_chain_struct_region';
	my ($aref, $struct_aref, $chain_aref) = @_;
	my @return;
	my @chain_return;
	my $s1 = scalar(@$aref);
	my $s2 = scalar(@$struct_aref);
	my $warning = 0;

	if (scalar(@$struct_aref) != scalar(@$chain_aref)) { 
		printf "ERROR! $sub: Inputs not the same size: %i != %i\n", scalar(@$struct_aref), scalar(@$chain_aref);
	}

	if ($DEBUG > 1) { 
		#printf "DEBUG s1:%i s2:%i\n", $s1, $s2;
		print "DEBUG: Struct_reference ID ". rangify(@$aref) . "\n";
		print "DEBUG: Struct_seq_reference ID " . rangify(@$struct_aref) ."\n";
	}
	for (my $i = 0; $i < scalar(@$aref); $i++) {
		if (defined($$struct_aref[$$aref[$i]-1])) { 
			push(@return, $$struct_aref[$$aref[$i]-1]) ;
			push(@chain_return, $$chain_aref[$$aref[$i]-1]);
		}else{
			#print Dumper($struct_aref);
			#print "WARNING: $sub: residue not defined or not structured-> $i $$aref[$i]\n";
			$warning++; 
		}
	}  
	if ($warning > 1)  {
		print "WARNING! $sub: $warning residues not defined or structured\n";
	}
	if ($DEBUG > 1 ) { 
		print "DEBUG: return " . rangify(@return) . "\n";
	}
	return (\@return, \@chain_return, $warning);
}


#Split a scop range (i.e. [A:1-100, A:150-200]) into its components.
sub scop_range_split { 
	my $sub = 'scop_range_split';
	my ($scop_range_str, $no_ins) = @_;

	my @segs;
	my %chains;
	my @chains;
	if ($scop_range_str =~ /[A-Za-z0-9]+:\-?\d+[A-Z]?\-\d+[A-Z]?/) { 
		while ($scop_range_str =~ /([A-Za-z0-9]+):(\-?\d+)([A-Z]?)\-(\-?\d+)([A-Z]?)/g) { 
			
			my $chain = $1;
			my $seg_start	= $2;
			my $seg_start_ins_code	= $3;
			if ($seg_start_ins_code) { 
				print "WARNING! $sub: Insertion code $seg_start_ins_code on $scop_range_str\n";
				if ($no_ins) { 
					#?
				}else{
					$seg_start = $seg_start . $seg_start_ins_code
				}
			}
			my $seg_end	= $4;
			my $seg_end_ins_code	= $5;
			if ($seg_end_ins_code) { 
				print "WARNING! $sub: Insertion code $seg_end_ins_code on $scop_range_str\n";
				if ($no_ins) { 
				}else{
					$seg_end = $seg_end . $seg_end_ins_code;
				}
			}
			#print "DEBUG! $sub: $chain $seg_start $seg_end\n";
			push (@segs, "$seg_start-$seg_end");
			push (@chains, "$chain");
			$chains{$chain}++;
		}
	}elsif($scop_range_str =~ /([A-Za-z0-9]+):(\-?\d+[A-Z]?)/) { 
		push (@chains, $1);
		push (@segs, $2);
		$chains{$1}++;
	}

	my $range_str = join(",", @segs);
	my $chain_str;
	if (scalar(keys(%chains)) > 1) { 
		#print "WARNING! $sub: multi-chain range parse for $scop_range_str, take heed!\n";

		my @sort_chains = sort { $a cmp $b } keys %chains;
		$chain_str = join(",", @chains);
	}else{
		$chain_str = $chains[0];
	}
	if (!defined($range_str) || !defined($chain_str)) { 
		die "$sub: Failed expansion on $scop_range_str\n";
	}
	return ($range_str, $chain_str); 
}


sub pdb_range_expand { 
	my $sub = 'pdb_range_expand';

	my ($pdb_range_str, $pdbnum_aref, $pdb, $chain) = @_;

	if (!$pdb_range_str || !$pdbnum_aref) { 
		die "ERROR $sub: inputs $pdb_range_str $pdbnum_aref\n";
	}

	my %pdbnum_to_seqid;
	for (my $i = 1; $i < scalar(@$pdbnum_aref); $i++) { 
		if (defined($$pdbnum_aref[$i]) && $$pdbnum_aref[$i] =~ /\d+/ && !defined($pdbnum_to_seqid{$$pdbnum_aref[$i]})) { 
			$pdbnum_to_seqid{$$pdbnum_aref[$i]} = $i;
		}
	}

	my @seqid;
	while ($pdb_range_str =~ /(\-?\d+[A-Z]?)\-(\-?\d+[A-Z]?)/g) { 

		my $start 	= $1;
		my $end 	= $2;
		if ($DEBUG) { 
			print "DEBUG $sub: $pdb_range_str $start $end $pdbnum_to_seqid{$start} $pdbnum_to_seqid{$end}\n";
		}

		if ($start eq $end) { next } 
		if (!$pdbnum_to_seqid{$start}) { 
			#print "WARNING! $sub: PDB number $start not found, skipping segment\n";
			#next;
			print "WARNING! $sub: PDB number $start not found, attempting fast forward\n";
			for (my $i = $start; $i < $end; $i++) { 
				print "$sub FFORWARD $i\n";
				if ($pdbnum_to_seqid{$i}) { 
					$start = $i;
					last;
				}
			}
			if (!$pdbnum_to_seqid{$start}) { 
				print "ERROR! $sub: Forward failed, skipping segement\n";
				next;
			}
			#return 0;
		}
		if (!$pdbnum_to_seqid{$end}) { 
			print "WARNING! $sub PDB num $end not found, attempting rewind\n";
			for (my $i = $end; $i > $start; $i--) { 
				print "$sub REWIND $i\n";
				if ($pdbnum_to_seqid{$i}) { 
					$end = $i;
					last;
				}
			}
			if (!$pdbnum_to_seqid{$end}) { 
				print "WARNING! $sub: PDB number $end not found, skipping segement\n";
				next;
			}
			#return 0;
		}
		for (my $i = $pdbnum_to_seqid{$start}; $i < $pdbnum_to_seqid{$end} + 1; $i++) { 
			push (@seqid, $i);
		}
	}
	return \@seqid;
}

sub scopify_range { 
	my $sub = 'scopify_range';

	my ($range_str, $chain) = @_;

	if ($DEBUG) { 
		print "DEBUG $sub: range_str $range_str chain $chain\n";
	}

	my @segs;
	while ($range_str =~ /(\-?\d+[A-Z]?\-\-?\d+[A-Z]?)/g) { 
		my $seg = $1;
		my $scop_seg_str .= "$chain:$seg";
		push (@segs, $scop_seg_str);
	}

	my $str = join(",", @segs);
	#my $scop_str = "[$str]";
	return $str;
}

		

sub pdb_rangify { 
	my $sub = 'pdb_rangify';

	my ($seqid_aref, $pdbnum_aref) = @_;

	if ($DEBUG) { 
		print "DEBUG $sub: seqid_aref $seqid_aref, $pdbnum_aref\n";
	}

	my $seqid_range_str = rangify(@$seqid_aref);
	if ($DEBUG ) { 
		print "DEBUG $sub: seqid range str $seqid_range_str\n";
	}

	my @pdbnum_range_segs;
	while ($seqid_range_str =~ /((\-?\d+)[A-Z]?\-(\d+)[A-Z]?)/g) { 

		my $seg			= $1;
		my $seg_start_seqid	= $2;
		my $seg_end_seqid	= $3;

		if (defined($$pdbnum_aref[$seg_start_seqid]) && defined($$pdbnum_aref[$seg_end_seqid])) { 
			if ($DEBUG) { 
				print "DEBUG $sub: start pdb num $$pdbnum_aref[$seg_start_seqid] end pdb num $$pdbnum_aref[$seg_end_seqid]\n";
			}
			push (@pdbnum_range_segs, "$$pdbnum_aref[$seg_start_seqid]-$$pdbnum_aref[$seg_end_seqid]");
		}else{
			if (!defined($$pdbnum_aref[$seg_start_seqid]) && defined($$pdbnum_aref[$seg_end_seqid])) { 
				print "WARNING! $sub: No pdb numbers for seg $seg\n";
			}elsif(!defined($$pdbnum_aref[$seg_start_seqid])) { 
				print "WARNING! $sub: No start pdbnum for seg $seg\n";
			}elsif(!defined($$pdbnum_aref[$seg_end_seqid])) { 
				print "WARNING! $sub: No end pdbnum for seg $seg\n";
			}else{
				die "ERROR! $sub: Catastrophic regexp failure on $seqid_range_str\n";
			}
		}
	}
	my $pdbnum_range_str = join(",", @pdbnum_range_segs);

	if ($DEBUG) { 
		print "DEBUG: $sub: pdbnum range str $pdbnum_range_str\n";
	}
	return $pdbnum_range_str;
} 

sub multi_chain_pdb_rangify { 
	my $sub = 'multi_chain_pdb_rangify';
	my ($pos_aref, $pdbnum_href, $chain_aref) = @_;

	if (!$pos_aref) { 
		print "WARNING! $sub: pos_aref undefined, returning...\n";
		return 0;
	}
	if (!$chain_aref) { 
		print "WARNING! $sub: chain_aref undefined, returning...\n";
		return 0;
	}
	if (!$pdbnum_href) { 
		print "WARNING! $sub: pdbnum_href undefined, returning...\n";
		return 0;
	}

	if (scalar(@$pos_aref) != scalar(@$chain_aref)) { 
		my $a = scalar(@$pos_aref);
		my $b = scalar(@$chain_aref);
		die "ERROR $sub: Inputs not same size $a != $b\n";
	}

	my $i = 0;
	my $range;
	my $on = 0;
	while ($i < scalar(@$pos_aref)) {
		if ($i == 0) { 
			if (!defined($$pdbnum_href{$$chain_aref[0]}{$$pos_aref[0]})) { 
				#print "$sub: $$chain_aref[0] $$pos_aref[0]?\n";
				#print Dumper($pos_aref);
				#print Dumper($pdbnum_href);
				#die;
				return 0;
			}
			$range .= "$$chain_aref[0]:$$pdbnum_href{$$chain_aref[0]}{$$pos_aref[0]}";
			if ($$pos_aref[$i+1] == $$pos_aref[$i] + 1 && $$chain_aref[$i+1] eq $$chain_aref[$i]) { 
				$range .= "-";
				$on = 1;
			}else{
				$range .= ",";
				$on = 0;
			}
		}elsif ($i < scalar(@$pos_aref) -1) {
			if ($on) { 
				if ($$pos_aref[$i-1] +1 == $$pos_aref[$i] && $$chain_aref[$i-1] eq $$chain_aref[$i]) { 
					$i++;
					next;
				}elsif ($$pos_aref[$i-1] + 1 < $$pos_aref[$i] || $$chain_aref[$i-1] ne $$chain_aref[$i]) { 
					$range .= $$pdbnum_href{$$chain_aref[$i-1]}{$$pos_aref[$i-1]}  . ",";
					if ($$pos_aref[$i+1] eq $$pos_aref[$i]+1) { 
						$range .= "$$chain_aref[$i]:$$pdbnum_href{$$chain_aref[$i]}{$$pos_aref[$i]}-";
					}elsif ($$pos_aref[$i+1] > $$pos_aref[$i] + 1 || $$chain_aref[$i+1] ne $$chain_aref[$i]) { 
						if (defined $$pdbnum_href{$$chain_aref[$i]}{$$pos_aref[$i]} ) { 
							$range .= "$$chain_aref[$i]:$$pdbnum_href{$$chain_aref[$i]}{$$pos_aref[$i]},";
						}elsif (defined $$pdbnum_href{$$chain_aref[$i-1]}{$$pos_aref[$i-1]}) { 
							$range .= "$$chain_aref[$i]:$$pdbnum_href{$$chain_aref[$i-1]}{$$pos_aref[$i-1]},";
						}else{
							die "ERROR! $sub: No pdbnum found for $$chain_aref[$i] $$pos_aref[$i]\n";
						}
						$on = 0; 
					}elsif ($$pos_aref[$i+1] == $$pos_aref[$i]) { 
						$i++;
						next;
					}else{

						printf "ERROR! $sub: i-1->(%i %i %s) i->(%i %i %s) i+1->(%i %i %s)\n", $i-1, $$pos_aref[$i-1], $$chain_aref[$i-1], $i, $$pos_aref[$i], $$chain_aref[$i], $i+1, $$pos_aref[$i+1], $$chain_aref[$i+1];
					}
				}else{
					#die "?$hora_id: $pos[$i-1] $pos[$i]\n";
					print "WARNING: possible ins? $$pos_aref[$i-1] $$pos_aref[$i] $$chain_aref[$i-1] $$chain_aref[$i]\n";
					$i++;
					next;
				}
			}else{
				if ($$pos_aref[$i+1] eq $$pos_aref[$i] + 1) { 
					$range .= "$$chain_aref[$i]:$$pdbnum_href{$$chain_aref[$i]}{$$pos_aref[$i]}-";
					$on = 1;
				}elsif($$pos_aref[$i+1] > $$pos_aref[$i] + 1 || $$chain_aref[$i+1] ne $$chain_aref[$i]) { 
					if (defined $$pdbnum_href{$$chain_aref[$i]}{$$pos_aref[$i]} ) { 
						$range .= "$$chain_aref[$i]:$$pdbnum_href{$$chain_aref[$i]}{$$pos_aref[$i]},";
					}elsif (defined $$pdbnum_href{$$chain_aref[$i-1]}{$$pos_aref[$i-1]}) { 
						$range .= "$$chain_aref[$i-1]:$$pdbnum_href{$$chain_aref[$i-1]}{$$pos_aref[$i-1]},";
					}else{
						die "ERROR! $sub: No pdbnum found for $$chain_aref[$i] $$pos_aref[$i]\n";
					}
					$on = 0;
				}elsif($$pos_aref[$i+1] == $$pos_aref[$i] ) { 
					$i++;
					next;
				}else{ 
					die "?!\n";
				}
			}
		}else{ 
			$range .= $$pdbnum_href{$$chain_aref[-1]}{$$pos_aref[-1]};
			#$range .= $$pdbnum_aref[-1];
		}
		$i++;
	}
	return $range;
}



sub rangify { 
	my $sub = 'rangify';
	my @pos = @_;

	my $i = 0;
	my $range;
	my $on = 0;
	
	if (scalar(@pos) == 0) { 
		return 0;
		#die "ERROR! $sub: Empty input array\n";
	}

	if (scalar(@pos) == 1) { 
		return 0;
		die "ERROR! $sub: Singleton input array\n";
	}

	while ($i < scalar(@pos)) {
		if ($i == 0) { 
			$range .= $pos[0] ;
			if ($pos[$i+1] == $pos[$i] + 1) { 
				$range .= "-";
				$on = 1;
			}else{
				$range .= ",";
				$on = 0;
			}
		}elsif ($i < scalar(@pos) -1) {
			if ($on) { 
				if ($pos[$i-1] +1 == $pos[$i]) { 
					$i++;
					next;
				}elsif ($pos[$i-1] + 1 < $pos[$i]) { 
					$range .= $pos[$i-1]  . ",";
					if ($pos[$i+1] eq $pos[$i]+1) { 
						$range .= $pos[$i] . "-"; #Drops sometimes
					}elsif ($pos[$i+1] > $pos[$i] + 1) { 
						$range .= $pos[$i] . ",";
						$on = 0;
					}elsif ($pos[$i+1] == $pos[$i]) { 
						$i++;
						next;
					}else{
						die "ERROR! $sub: i-1:" . ($i-1) ." $pos[$i-1] i: $i $pos[$i] i+1: ". ($i+1) ." $pos[$i+1]\n";
					}
				}else{
					#die "?$hora_id: $pos[$i-1] $pos[$i]\n";
					print "WARNING: possible ins? $pos[$i-1] $pos[$i]\n";
					$i++;
					next;
				}
			}else{
				if ($pos[$i+1] eq $pos[$i] + 1) { 
					$range .= $pos[$i] . "-";
					$on = 1;
				}elsif($pos[$i+1] > $pos[$i] + 1) { 
					$range .= $pos[$i] .",";
					$on = 0;
				}elsif($pos[$i+1] == $pos[$i] ) { 
					$i++;
					next;
				}else{ 
					die "?!\n";
				}
			}
		}else{ 
			if ($on) { 
				if ($pos[$i-1]+1 == $pos[$i]) { 
					$range .= $pos[$#pos];
				}elsif($pos[$i-1]+1 < $pos[$i]) { 
					$range .= $pos[$i-1] . "," . $pos[$#pos];
				}else{
					printf "WARNING! Range terminus: %i %i\n", $pos[$i-1]+1, $pos[$i];
				}
			}else{
				$range .=  $pos[$#pos];
			}

		}
		$i++;
	}
	return $range;
}
sub multi_chain_rangify { 
	my $sub = 'multi_chain_rangify';
	my ($pos_aref, $chain_aref) = @_;

	if (scalar(@$pos_aref) != scalar(@$chain_aref)) { 
		my $a = scalar(@$pos_aref);
		my $b = scalar(@$chain_aref);
		die "ERROR $sub: Inputs not same size $a != $b \n";
	}

	my $i = 0;
	my $range;
	my $on = 0;
	while ($i < scalar(@$pos_aref)) {
		if ($i == 0) { 
			if (!defined($$chain_aref[0])) { return 0 } 
			$range .= "$$chain_aref[0]:$$pos_aref[0]";
			if ($$pos_aref[$i+1] == $$pos_aref[$i] + 1 && $$chain_aref[$i+1] eq $$chain_aref[$i]) { 
				$range .= "-";
				$on = 1;
			}else{
				$range .= ",";
				$on = 0;
			}
		}elsif ($i < scalar(@$pos_aref) -1) {
			if ($on) { 
				if ($$pos_aref[$i-1] +1 == $$pos_aref[$i] && $$chain_aref[$i-1] eq $$chain_aref[$i]) { 
					$i++;
					next;
				}elsif ($$pos_aref[$i-1] + 1 < $$pos_aref[$i] || $$chain_aref[$i-1] ne $$chain_aref[$i]) { 
					$range .= $$pos_aref[$i-1]  . ",";
					if ($$pos_aref[$i+1] eq $$pos_aref[$i]+1) { 
						$range .= "$$chain_aref[$i]:$$pos_aref[$i]-";
					}elsif ($$pos_aref[$i+1] > $$pos_aref[$i] + 1 || $$chain_aref[$i+1] ne $$chain_aref[$i]) { 
						$range .= "$$chain_aref[$i]:$$pos_aref[$i],";
						$on = 0;
					}elsif ($$pos_aref[$i+1] == $$pos_aref[$i]) { 
						$i++;
						next;
					}else{
						printf "ERROR! $sub: i-1->(%i %i %s) i->(%i %i %s) i+1->(%i %i %s)\n", $i-1, $$pos_aref[$i-1], $$chain_aref[$i-1], $i, $$pos_aref[$i], $$chain_aref[$i], $i+1, $$pos_aref[$i+1], $$chain_aref[$i+1];
						return 0 ;
						#die "ERROR! $sub: i-1:" . $i-1 ."$$pos_aref[$i-1] i: $i $$pos_aref[$i] i+1: ". $i+1 ." $$pos_aref[$i+1]\nj\n";
					}
				}else{
					#die "?$hora_id: $pos[$i-1] $pos[$i]\n";
					print "WARNING: possible ins? $$pos_aref[$i-1] $$chain_aref[$i-1]  $$pos_aref[$i] $$chain_aref[$i]\n";
					$i++;
					next;
				}
			}else{
				if ($$pos_aref[$i+1] eq $$pos_aref[$i] + 1) { 
					$range .= "$$chain_aref[$i]:$$pos_aref[$i]-";
					$on = 1;
				}elsif($$pos_aref[$i+1] > $$pos_aref[$i] + 1 || $$chain_aref[$i+1] ne $$chain_aref[$i]) { 
					$range .= $$pos_aref[$i] .",";
					$on = 0;
				}elsif($$pos_aref[$i+1] == $$pos_aref[$i] ) { 
					$i++;
					next;
				}else{ 
					die "?!\n";
				}
			}
		}else{ 
			$range .= $$pos_aref[-1];
		}
		$i++;
	}
	return $range;
}
sub multi_chain_pdb_range_expand { 
	my $sub = 'multi_chain_pdb_range_expand';

	my ($range_str, $seqid_aref, $chain_aref, $pdbnum_href) = @_;
	$range_str =~ /\w:\-?\d+\-\d+/ or return;

	my %pdbnum_to_seqid;
	for (my $i = 0; $i < scalar(@$seqid_aref); $i++) { 
		my $seqid = $$seqid_aref[$i];
		my $chain = $$chain_aref[$i];
		if (exists($$pdbnum_href{$chain}{$seqid})) { 
			my $pdbnum = $$pdbnum_href{$chain}{$seqid};
			$pdbnum_to_seqid{$chain}{$pdbnum} = $seqid;
		}
	}

		

	my @range;
	my @chain;
	my @segs = split(/,/, $range_str);
	foreach my $seg (@segs) { 
		if ($seg =~ /(\w):(\-?\d+\w?)\-(\-?\d+\w?)/) { 
			my $chain = $1;
			my $start = $2;
			my $end	= $3;

			if ($start eq $end) { next } 
			if (!$pdbnum_to_seqid{$chain}{$start}) { 
				print "WARNING! $sub: PDB number $start not found,  attempting fast forward\n";
				for (my $i = $start; $i < $end; $i++) { 
					print "$sub FFORWARD $i\n";
					if ($pdbnum_to_seqid{$chain}{$i}) { 
						$start = $i;
						last;
					}
				}
				if (!$pdbnum_to_seqid{$chain}{$start}) { 
					print "ERROR! $sub: Forward failed, skipping segment\n";
					next;
				}
			}

			if (!$pdbnum_to_seqid{$chain}{$end}) { 
				print "WARNING! $sub: PDB number $end not found, attempting rewind\n";
				
				for (my $i = $end; $i > $start; $i--) { 
					print "$sub REWIND $i\n";
					if ($pdbnum_to_seqid{$chain}{$i}){ 
						$end = $i;
						last;
					}
				}
				if (!$pdbnum_to_seqid{$chain}{$end}) { 
					print "WARNING! $sub: PDB number $end not found, skipping segment\n";
					next;
				}
			}
			for (my $i = $pdbnum_to_seqid{$chain}{$start}; $i < $pdbnum_to_seqid{$chain}{$end}+1; $i++) { 
				push (@range, $i);
				push (@chain, $chain);
			}
		}else{ 
			die "ERROR! $sub: seg $seg\n";
		}
	}
	return (\@range, \@chain);
	
}
sub multi_chain_range_expand { 
	my $sub = 'multi_chain_range_expand';

	my ($range_str) = @_;
	$range_str =~ /\w:\-?\d+\-\d+/ or return;

	my @range;
	my @chain;
	my @segs = split(/,/, $range_str);
	foreach my $seg (@segs) { 
		if ($seg =~ /(\w+):(\-?\d+)\-(\-?\d+)/) { 
			my $chain = $1;
			my $start = $2;
			my $end	= $3;
			for (my $i = $start; $i < $end +1; $i++) { 
				push (@range, $i);
				push (@chain, $chain);
			}
		}elsif($seg =~ /(\w+):(\d+)/) { 
			push (@range, $2);
			push (@chain, $1);
		}else{ 
			print "WARNING! $sub: seg $seg\n";
		}
	}
	return (\@range, \@chain);
}
sub range_expand_obsolete { 
	my $sub = 'range_expand_obsolete';
	my ($range_str) = @_;
	
	$range_str =~ /\d+/ or return;
	my @range;
	my @segs = split(/,/, $range_str); 
	foreach my $seg (@segs) { 
		if ($seg =~ /(\-?\d+)\-(\-?\d+)/) { 
			my $start = $1;
			my $end = $2;
			for (my $i = $start; $i < $end + 1; $i ++) { 
				push (@range, $i);
			}
		}elsif($seg =~ /\d+/) { 
			print "??$seg\n";
			push (@range,$seg);
		}else{
			#die "ERROR! $sub: Bad segment regexp: $seg\n";
			print "WARNING! $sub: Bad segment regexp: $seg $range_str\n";
			return 0;
		}
	}
	return \@range;
}

sub range_expand { 
	my $sub = 'range_expand';
	my ($range_str) = @_;
	
	$range_str =~ /\-?\d+/ or return;
	my @range;
	my @segs = split(/,/, $range_str); 
	my %seen;
	foreach my $seg (@segs) { 
		if ($seg =~ /(\-?\d+)\-(\-?\d+)/) { 
			my $start = $1;
			my $end = $2;
			for (my $i = $start; $i < $end + 1; $i ++) { 
				if (!$seen{$i}) { 
					push (@range, $i);
				}
				$seen{$i}++;
			}
		}elsif($seg =~ /\d+/) { 
			push (@range,$seg);
		}else{
			#die "ERROR! $sub: Bad segment regexp: $seg\n";
			print "WARNING! $sub: Bad segment regexp: $seg $range_str\n";
			return 0;
		}
	}
	@range = sort {$a <=> $b} @range;
	return \@range;
}

sub multi_chain_range_include { 
	my ($range1_aref, $chain1_aref, $range2_aref, $chain2_aref, $sort_href) = @_;

	if (scalar(@$range1_aref) != scalar(@$chain1_aref)) { die } 
	if (scalar(@$range2_aref) != scalar(@$chain2_aref)) { die } 

	my %range2_lookup;
	my %seen;
	for (my $i = 0; $i < scalar(@$range2_aref); $i++) { 
		my $range2 = $$range2_aref[$i];
		my $chain2 = $$chain2_aref[$i];
		$seen{$chain2}{$range2}++;
	}

	for (my $i = 0; $i < scalar(@$range1_aref); $i++) { 
		my $range1 = $$range1_aref[$i];
		my $chain1 = $$range1_aref[$i];
		if (!$seen{$chain1}{$range1}) { 
			push (@$range2_aref, $$range1_aref[$i]);
			push (@$chain2_aref, $$chain1_aref[$i]);
		}

	}
	my @sort_array;
	for (my $i = 0; $i<scalar(@$range2_aref); $i++) { 
		$sort_array[$i]{chain}	= $$chain2_aref[$i];
		$sort_array[$i]{seqid}	= $$range2_aref[$i];
	}

	@sort_array = sort {    
						my $anum = $a->{seqid}; my $bnum = $b->{seqid}; 
						my $achain = $a->{chain}; my $bchain = $b->{chain};
						$$sort_href{$achain}{$anum} <=> $$sort_href{$bchain}{$bnum};
						} @sort_array;

	undef ($range2_aref);
	undef ($chain2_aref);
	for (my $i = 0; $i < scalar(@sort_array); $i++) { 
		push (@$range2_aref, $sort_array[$i]{seqid});
		push (@$chain2_aref, $sort_array[$i]{chain});
	}

	return 1;
}

sub range_include { 
	my ($range1_aref, $range2_aref) = @_;
	for (my $i = 0; $i < scalar(@$range1_aref); $i++) { 	
		if (! grep {$_ == $$range1_aref[$i]} @$range2_aref) {
			push (@$range2_aref, $$range1_aref[$i]);
		}
	}
	@$range2_aref = sort {$a <=> $b} @$range2_aref;
	return 1;
}

sub multi_chain_range_exclude { 
	my $sub = 'multi_chain_range_exclude';
	my ($range1_aref, $chain1_aref, $range2_aref, $chain2_aref, $sort_href) = @_;
	if ($DEBUG) {
		printf "DEBUG $sub: r1_rangify %s r2_rangify %s\n", multi_chain_rangify($range1_aref, $chain1_aref), multi_chain_rangify($range2_aref, $chain2_aref);
		printf "DEBUG $sub: %i %i %i %i\n", scalar(@$range1_aref), scalar(@$chain1_aref), scalar(@$range2_aref), scalar(@$chain2_aref);
	}
	if (scalar(@$range1_aref) != scalar(@$chain1_aref)) { die } 
	if (scalar(@$range2_aref) != scalar(@$chain2_aref)) { die } 

	my %seen;
	for (my $i = 0; $i < scalar(@$range1_aref); $i++) { 
		if (!defined($$range1_aref[$i])) { die } 
		if (!defined($$chain1_aref[$i])) { die } 

		$seen{$$chain1_aref[$i]}{$$range1_aref[$i]}++;
		#print "SEEN $$range1_aref[$i] $$chain1_aref[$i]\n";
	}
	for (my $i = 0; $i < scalar(@$range2_aref); $i++) { 
		if ($seen{$$chain2_aref[$i]}{$$range2_aref[$i]} ) { 
			if ($DEBUG > 1) { 
				print "DEBUG $sub: DELETE $i $$range2_aref[$i] $$chain2_aref[$i]\n";
			}
			#print "DELETE $i $$range2_aref[$i] $$chain2_aref[$i]\n";
			delete($$range2_aref[$i]);
			delete($$chain2_aref[$i]);
		}
	}
	my @sort_array;
	my $pos = 0;
	printf "%i %i\n", scalar(@$chain2_aref), scalar(@$range2_aref);
	for (my $i = 0; $i<scalar(@$chain2_aref); $i++) { 
		if (defined($$chain2_aref[$i])) { 
			$sort_array[$pos]{chain}	= $$chain2_aref[$i];
			$sort_array[$pos]{seqid}	= $$range2_aref[$i];
			#print "? $pos $$chain2_aref[$i] $$range2_aref[$i]\n";
			$pos++;
		}
	}

	@sort_array = sort {    my $anum = $a->{seqid}; 	my $bnum = $b->{seqid}; 
				my $achain = $a->{chain}; 	my $bchain = $b->{chain};
				if (!$$sort_href{$achain}{$anum}) { 
					print "WARNING $sub: No value for ach $achain anum $anum\n"
				}
				if (!$$sort_href{$bchain}{$bnum}) { 
					print "WARNING $sub: No value for ach $bchain anum $bnum\n"
				}
				$$sort_href{$achain}{$anum} <=> $$sort_href{$bchain}{$bnum} 
				} @sort_array;

	undef (@$range2_aref);
	undef (@$chain2_aref);
	for (my $i = 0; $i < scalar(@sort_array); $i++) { 
		push (@$range2_aref, $sort_array[$i]{seqid});
		push (@$chain2_aref, $sort_array[$i]{chain});
	}

	my $i = 0;
	while (defined($$range2_aref[$i])) {
		$i++;
	}
	print "DEBUG $sub: i $i\n";
	splice(@$range2_aref, $i);
	splice(@$chain2_aref, $i);

	if ($DEBUG) { 
		printf "DEBUG $sub: r1_rangify %s r2_rangify %s r1_scalar %i r2_scalar %i\n", multi_chain_rangify($range1_aref, $chain1_aref), multi_chain_rangify($range2_aref, $chain2_aref), scalar(@$range1_aref), scalar(@$range2_aref);
	}
	return 1;

}



sub range_exclude { 
	my $sub = 'range_exclude';
	my ($range1_aref, $range2_aref) = @_;
	if ($DEBUG) { 
		printf "DEBUG $sub: r1_rangify %s r2_rangify %s\n", rangify(@$range1_aref), rangify(@$range2_aref);
	}
	my %seen;
	for (my $i = 0; $i < scalar(@$range1_aref); $i++) { 
		$seen{$$range1_aref[$i]}++;
	}
	for (my $i = 0; $i < scalar(@$range2_aref); $i++) { 
		if ($seen{$$range2_aref[$i]}) { 
			if ($DEBUG > 1) { 
				print "DEBUG $sub: DELETE $i $$range2_aref[$i]\n";
			}
			delete($$range2_aref[$i]);
		}
	}
	@$range2_aref = sort{
			my ($anum, $bnum);
			if (!$a) { $anum = 9999 } else { $anum = $a }
			if (!$b) { $bnum = 9999 } else { $bnum = $b }
			$anum <=> $bnum} @$range2_aref;
	my $i = 0;
	while (defined($$range2_aref[$i])) {
		$i++;
	}
	splice(@$range2_aref, $i);
	if ($DEBUG) { 
		printf "DEBUG $sub: r1_rangify %s r2_rangify %s\n", rangify(@$range1_aref), rangify(@$range2_aref);
	}
	return 1;
}

sub region_coverage { 
	my $sub = 'region_coverage';
	my ($range1_aref, $range2_aref) = @_;
	
	if (!$range1_aref || !$range2_aref) { 
		print "WARNING! $sub: undefined range? skipping...\n";
		return 0;
	}
	
	#Calculate intersection of range arrays, express as coverage of range1. Do averaging outside of subroutine
	my $length1 = scalar(@$range1_aref);

	my %isect;
	my %union;
	foreach my $pos1 (@$range1_aref) { $union{$pos1} = 1}
	foreach my $pos2 (@$range2_aref) { 
		if ($union{$pos2} ) { $isect{$pos2} = 1}
		$union{$pos2} = 1;
	}
	my @isect = sort {$a<=>$b} keys %isect;
	if ($length1 > 0) { 
		my $coverage = scalar(@isect) / $length1;
		if ($DEBUG) { 
			printf "DEBUG $sub: %i %i %0.2f\n", scalar(@isect), $length1, scalar(@isect) / $length1;
		}
		return $coverage;
	}else{
		return 0;
	}
}

sub residue_coverage {
	my ($range1_aref, $range2_aref) = @_;


	my $length1 = scalar(@$range1_aref);
	
	my %isect;
	my %union;
	foreach my $pos1 (@$range1_aref) { $union{$pos1} = 1}
	foreach my $pos2 (@$range2_aref) { 
		if ($union{$pos2} ) { $isect{$pos2} = 1}
		$union{$pos2} = 1;
	}
	my @isect = sort {$a<=>$b} keys %isect;
	if ($length1 > 0) { 
		#my $coverage = scalar(@isect) / $length1;
		my $res_count = scalar(@isect);
		return $res_count;
	}else{
		return 0;
	}
}

sub multi_chain_residue_coverage { 
	my ($range1_aref, $chain1_aref, $range2_aref, $chain2_aref) = @_;

	my $length1 = scalar(@$range1_aref);

	my %isect;
	my %union;
	for (my $i = 0; $i < scalar(@$range1_aref); $i++) { 
		my $pos1 = $$range1_aref[$i];
		my $chain1 = $$chain1_aref[$i];
		my $comp_key = $pos1 . $chain1;

		$union{$comp_key} = 1;
	}
	for (my $i = 0; $i < scalar(@$range2_aref); $i++) { 
		my $pos2 = $$range2_aref[$i];
		my $chain2 = $$chain2_aref[$i];
		my $comp_key = $pos2 . $chain2;
		if ($union{$comp_key}) { $isect{$comp_key} = 1 } 
		$union{$comp_key} = 1;
	}
	my @isect = sort {$a cmp $b} keys %isect;
	if ($length1 > 0) { 
		my $res_count = scalar(@isect);
		return $res_count;
	}else{
		return 0;
	}
}

sub multi_chain_region_coverage { 
	my $sub = 'multi_chain_region_coverage';
	my ($range1_aref, $chain1_aref, $range2_aref, $chain2_aref) = @_;
	my $length1 = scalar(@$range1_aref);

	if ($DEBUG) { 
		printf "DEBUG top: $sub %i %i %i %i\n", scalar(@$range1_aref), scalar(@$chain1_aref), scalar(@$range2_aref), scalar(@$chain2_aref);
	}

	if (scalar(@$range1_aref) != scalar(@$chain1_aref)) { croak "ERROR! $sub: range 1/chain 1 array not same size\n" } 
	if (scalar(@$range2_aref) != scalar(@$chain2_aref)) { croak "ERROR! $sub: range 2/chain 2 array not same size\n" } 

	my %isect;
	my %union;
	for (my $i = 0; $i < scalar(@$range1_aref); $i++) { 
		my $pos1 = $$range1_aref[$i];
		my $chain1 = $$chain1_aref[$i];
		my $comp_key = $pos1 . ":" . $chain1;

		$union{$comp_key} = 1;
	}
	for (my $i = 0; $i < scalar(@$range2_aref); $i++) { 
		my $pos2 = $$range2_aref[$i];
		my $chain2 = $$chain2_aref[$i];
		if (!defined($pos2) || !defined($chain2)) { #0 truth
			print "WARNING! $sub: pos chain mismatch for $i $pos2 $chain2\n";
			return 0;
		} 


		my $comp_key = $pos2 .":". $chain2;
		if ($union{$comp_key}) { $isect{$comp_key} = 1 } 
		$union{$comp_key} = 1;
	}
	my @isect = sort {$a cmp $b} keys %isect;
	if ($length1 > 0) { 
		if ($DEBUG ) { 
			printf "DEBUG: $sub %i %i %2.f\n", scalar(@isect), $length1, scalar(@isect)/$length1;
		}
		my $res_coverage = scalar(@isect) / $length1;
		return $res_coverage;
	}else{
		print "WARNING! $sub: lengt1 == 0?\n";
		return 0;
	}
	
}
sub multi_chain_region_coverage2 { 
	my $sub = 'multi_chain_region_coverage2';
	my ($range1_aref, $chain1_aref, $range2_aref, $chain2_aref) = @_;
	my $length1 = scalar(@$range1_aref);
	my $length2 = scalar(@$range2_aref);

	if ($DEBUG) { 
		printf "DEBUG top: $sub %i %i %i %i\n", scalar(@$range1_aref), scalar(@$chain1_aref), scalar(@$range2_aref), scalar(@$chain2_aref);
	}

	if (scalar(@$range1_aref) != scalar(@$chain1_aref)) { croak "ERROR! $sub: range 1/chain 1 array not same size\n" } 
	if (scalar(@$range2_aref) != scalar(@$chain2_aref)) { croak "ERROR! $sub: range 2/chain 2 array not same size\n" } 

	my %isect;
	my %union;
	for (my $i = 0; $i < scalar(@$range1_aref); $i++) { 
		my $pos1 = $$range1_aref[$i];
		my $chain1 = $$chain1_aref[$i];
		my $comp_key = $pos1 . ":" . $chain1;

		$union{$comp_key} = 1;
	}
	for (my $i = 0; $i < scalar(@$range2_aref); $i++) { 
		my $pos2 = $$range2_aref[$i];
		my $chain2 = $$chain2_aref[$i];
		if (!defined($pos2) || !defined($chain2)) { #0 truth
			print "WARNING! $sub: pos chain mismatch for $i $pos2 $chain2\n";
			return 0;
		} 


		my $comp_key = $pos2 .":". $chain2;
		if ($union{$comp_key}) { $isect{$comp_key} = 1 } 
		$union{$comp_key} = 1;
	}
	my @isect = sort {$a cmp $b} keys %isect;
	if ($length1 > 0) { 
		if ($DEBUG ) { 
			printf "DEBUG: $sub %i %i %2.f\n", scalar(@isect), $length1, scalar(@isect)/$length1;
		}
		my $res_coverage1 = scalar(@isect) / $length1;
		my $res_coverage2 = scalar(@isect) / $length2;
		return ($res_coverage1, $res_coverage2);
	}else{
		print "WARNING! $sub: lengt1 == 0?\n";
		return 0;
	}
	
}
sub region_boundary_byres { 
	my $sub = 'region_boundary_byres';
	my ($range1_str, $range2_str, $boundary_cutoff) = @_;

	my @range1_segs;
	my $i = 0;
	if ($range1_str =~ /\-?\d+\-\d+/) { 
		while ($range1_str =~ /(-?\d+)\-(\d+)/g) { 
			my $seg_start = $1;
			my $seg_end = $2;
			$range1_segs[$i]{start}	= $seg_start;
			$range1_segs[$i]{end}	 = $seg_end;
			$i++;
		}
	}else{
		die "ERROR! $sub: range1 string failed regep:$!\n";
	}

	my @range2_segs;
	$i = 0;
	if ($range2_str =~ /\-?\d+\-\d+/) { 
		while ($range2_str =~ /(-?\d+)\-(\d+)/g) { 
			my $seg_start = $1;
			my $seg_end = $2;
			$range2_segs[$i]{start}	= $seg_start;
			$range2_segs[$i]{end}	 = $seg_end;
			$i++;
		}
	}else{
		die "ERROR! $sub: range2 string failed regexcp:$!\n";
	}

	if ( scalar(@range1_segs) != scalar(@range2_segs)) { 
		return 0;
	}

	for (my $i = 0; $i < scalar(@range1_segs); $i++) { 
		if ( abs ($range1_segs[$i]{start} - $range2_segs[$i]{start}) > $boundary_cutoff) { 
			return 0;
		}
		if ( abs($range1_segs[$i]{end} - $range2_segs[$i]{end}) > $boundary_cutoff) { 
			return 0;
		}
	}
	return 1;

}
sub union_range_boundary {
	my ($range1_aref, $range2_aref) = @_;
	my %union;
	
	foreach my $s1 (@$range1_aref) { $union{$s1}++ } 
	foreach my $s2 (@$range2_aref) { $union{$s2}++ }

	my @union = sort {$a <=> $b} keys %union;

	my $first = $union[0];
	my $last  = $union[-1];
	return \@union;
}  
sub multi_chain_ungap_range { 
	my $sub = 'multi_chain_ungap_range';
	my ($range_str, $GAP_TOL) = @_;

	if ($GAP_TOL < 1) { 
		die "ERROR! $sub: gap tolerance not defined\n";
	}

	my @segs = split(/\,/, $range_str);

	my $i = 0;
	while ($i < scalar(@segs) -1) { 
		if ($segs[$i] =~ /(\w):(\-?\d+)\-(\-?\d+)/) { 
			my $chain1 	= $1;
			my $start1	= $2;
			my $end1 	= $3;

			if ($segs[$i+1] =~ /(\w):(\-?\d+)\-(\-?\d+)/) { 
				my $chain2 	= $1;
				my $start2	= $2;
				my $end2	= $3;

				if ($chain1 eq $chain2 && abs($start2 - $end1) < $GAP_TOL) { 
					splice(@segs, $i, 2, "$chain1:$start1-$end2");
				}else{
					$i++;
				}
			}
			elsif($segs[$i+1] =~ /(\w):(\d+)/) { 
				my $solo_chain2 = $1;
				my $solo2 = $2;
				if ($chain1 eq $solo_chain2 && abs($solo2 - $end1) < $GAP_TOL) { 
					splice(@segs, $i, 2, "$solo_chain2:$start1-$solo2")
				}else{	
					$i++;
				}
			}elsif($segs[$i+1] =~ /(\d+)/) { 
				#print "WARNING! $sub: Junky singleton seg $segs[$i+1], skipping\n";
				splice(@segs, $i, 2, "$segs[$i]");
			}else{
				print "WARNING! $sub: $segs[$i+1], huh?!\n";
				return 0;
			}
		}
		elsif($segs[$i] =~ /(\w):(\d+)/) { 
			my $solo_chain1 = $1;
			my $solo1 = $2;
			if ($segs[$i+1] =~ /(\w):(\-?\d+)\-(\-?\d+)/) { 
				my $chain2 = $1;
				my $start2 = $2;
				my $end2 = $3;
				if ($solo_chain1 eq $chain2 && abs($solo1 - $end2) < $GAP_TOL) { 
					splice(@segs, $i, 2, "$chain2:$solo1-$end2");
				}else{
					$i++;
				}
			}
			elsif($segs[$i+1] =~ /(\w):(\d+)/) { 
				my $solo_chain2 = $1;
				my $solo2 = $2;
				if ($solo_chain2 eq $solo2 && abs($solo2 - $solo1) < $GAP_TOL) { 
					splice(@segs, $i, 2, "$solo_chain2:$solo1-$solo2");
				}else{
					$i++;
				}
			}elsif($segs[$i+1] =~ /(\d+)/) { 
				#print "WARNING! $sub: Junky singleton seg $segs[$i+1], removing\n";
				splice(@segs, $i, 2, "$segs[$i]");
			}else{
				print "$sub: huh?! $segs[$i+1]\n";
				return 0;
			}
		}else{
			print "EH? $segs[$i]\n";
			$i++;
		}

	}
	$range_str = join(",", @segs);
	return $range_str;
}
 
sub ungap_range_aref { 
	my $sub = 'ungap_range_aref';
	my ($range_aref, $GAP_TOL) = @_;
	my $range = rangify($range_aref);
	$range = ungap_range($range, $GAP_TOL);
	$range_aref = range_expand($range);
}

sub ungap_range { 
	my $sub = 'ungap_range';
	my ($range_str, $GAP_TOL) = @_;

	if ($GAP_TOL < 1) { 
		confess "ERROR! $sub: gap tolerance not defined\n";
	}
	my @segs = split(/\,/, $range_str);

	
	my $i = 0;
	while ($i < scalar(@segs) - 1) { 
		if ($segs[$i] =~ /(\-?\d+)\-(\-?\d+)/) { 
			my $start1 = $1;
			my $end1 = $2;
			if ($segs[$i+1] =~ /(\-?\d+)\-(\-?\d+)/) { 
				my $start2 = $1;
				my $end2 = $2;
				if (abs($start2 - $end1) < $GAP_TOL) { 
					splice(@segs, $i, 2, "$start1-$end2");
				}else{
					$i++;
				}
			}
			elsif($segs[$i+1] =~ /(\d+)/) { 
				my $solo2 = $1;
				if (abs($solo2 - $end1) < $GAP_TOL) { 
					splice(@segs, $i, 2, "$start1-$solo2");
				}else{
					$i++;
				}
			}else{
				die "ungap_range: huh?! $range_str $segs[$i+1]\n";
			}
		}
		elsif($segs[$i] =~ /(\d+)/) { 
			my $solo1 = $1;
			if ($segs[$i+1] =~ /(\-?\d+)\-(\-?\d+)/) { 
				my $start2 = $1;
				my $end2 = $2;
				if (abs($solo1 - $start2) < $GAP_TOL) { 
					splice(@segs, $i, 2, "$solo1-$end2");
				}else{
					$i++;
				}
			}
			elsif($segs[$i+1] =~ /(\d+)/) { 
				my $solo2 = $1;
				if (abs($solo2 - $solo1) < $GAP_TOL) { 
					splice(@segs, $i, 2, "$solo1-$solo2");
				}else{
					$i++;
				}
			}else{
				die "ungap_range: huh?! $range_str $segs[$i+1]\n";
			}
		}else{
			die "ungap_range: huh? $range_str $segs[$i]\n";
		}
	}		
	$range_str = join(",", @segs);
	return $range_str;
}
1;
