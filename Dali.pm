package Domains::Dali;
require Exporter;

use warnings;
use strict;

use XML::LibXML;
use Domains::Range;

our @ISA =  qw(Exporter);

our @EXPORT =qw(dali_summ iset seq_grab blosum62_href_load adorn_rep_stats);

my $DEBUG = 0;

my $LD_LIBRARY_PATH = '/usr1/local/lib';
$ENV{LD_LIBRARY_PATH} .= $LD_LIBRARY_PATH;

our @aa = qw(ALA CYS ASP GLU PHE GLY HIS ILE LYS LEU MET ASN PRO GLN ARG SER THR VAL TRP TYR);
our %blosum_score_matrix = blosum62_href_load();

my $ABSOLUTE_MINIMUM = 3.0;

sub adorn_rep_stats { 
	my ($rep_stats_href, $range_cache_fn) = @_;

	open (my $fh, "<", $range_cache_fn) or die "ERROR! Could not open ref range cache $range_cache_fn for reading:$!\n";
	foreach my $ln (<$fh>) { 
		my @F = split(/\s+/, $ln);
		my $uid 		= $F[0];
		my $ecod_domain_id 	= $F[1];
		my $seqid_range 	= $F[2];
		my $pdb			= $F[3];
		my $chain 		= $F[4];
		if ($chain eq '.') { next } 
		my $short_uid		= substr($uid, 2, 5);
		if ($$rep_stats_href{$uid}) { 
			$$rep_stats_href{$uid}{range} = $seqid_range;
			$$rep_stats_href{$uid}{ecod_domain_id} = $ecod_domain_id;
			$$rep_stats_href{$uid}{pdb} = $pdb;
			$$rep_stats_href{$uid}{chain} = $chain;
		}
	}
}

sub dali_summ { 
	my $sub = 'dali_summ';

	my ($dali_files_aref, $dali_summ_fn, $query_href, $rep_stats) = @_;
	
	print "$sub: $dali_summ_fn\n";

	if ($DEBUG) { 
		printf "DEBUG: %i files in input, output to $dali_summ_fn\n", scalar(@$dali_files_aref);
	}

	my $dali_summ_xml_doc = XML::LibXML->createDocument();
	my $dali_summ_root_node	= $dali_summ_xml_doc->createElement('dali_summ_doc');
	$dali_summ_xml_doc->setDocumentElement($dali_summ_root_node);

	my $dali_hit_list_node	= $dali_summ_xml_doc->createElement('dali_hits');
	$dali_summ_root_node->appendChild($dali_hit_list_node);

	my %uid_lookup;
	foreach my $key (keys %$rep_stats) { 
		if (! exists $$rep_stats{$key}{ecod_domain_id} ) { 
			warn "WHAT? $key\n";
			next;
		}
		my $ecod_domain_id = $$rep_stats{$key}{ecod_domain_id};
		$uid_lookup{$ecod_domain_id} = $key;
	}


	my $i = 0;
	my $j = 0;
	my @hits;
	FILE:
	for (my $i = 0; $i < scalar(@$dali_files_aref); $i++) { 
		my $dali_file = $$dali_files_aref[$i];
		if ($DEBUG) { 
			print "DEBUG: $i $dali_file\n";
		}
		if (!-f $dali_file) { 
			print "WARNING! $sub: $dali_file not found, skipping...\n";
			next;
		}

		open (IN, $dali_file) or die "ERROR! $sub: Could not open $dali_file for reading:$!\n";

		my $ecod_domain_id;
		if ( $dali_file =~ /\.(e\w{4}.\d+)/) { 
			$ecod_domain_id = $1;
		}

		my ($hit_pdb, $hit_chain, $hit_uid);
		if ($dali_file =~ /\/e([0-9]\w{3})(\.)/) { #Go digging in ref for range to get parse
		#	print "$dali_file\n";
			#print "WARNING! MC PARSE not implemented go do it.\n";
			next;
		}elsif($dali_file =~ /(e([0-9]\w{3})(\w+))\./) {
			#$hit_pdb = $1;
			#$hit_chain = $2;
			$hit_uid = $uid_lookup{$ecod_domain_id};
		}elsif($dali_file =~ /(\d{9})/) { 
			$hit_uid = $1;
			$ecod_domain_id = $$rep_stats{$hit_uid}{ecod_domain_id};
		}else{
			die "So confused...\n";
		}
		$hit_pdb 	= $$rep_stats{$hit_uid}{pdb};
		$hit_chain 	= $$rep_stats{$hit_uid}{chain};
			
		#my ($hit_seqid_aref, $hit_struct_seqid_aref, $hit_pdbnum_aref, $hit_asym_id) = pdbml_seq_parse($hit_pdb, $hit_chain);
		my ($hit_seqid_aref, $hit_struct_seqid_aref, $hit_asym_id);
		if (exists $$rep_stats{$ecod_domain_id}) { #This should never happen 
			 $hit_seqid_aref 		= $$rep_stats{$ecod_domain_id}{seqid_aref};
			 $hit_struct_seqid_aref 	= $$rep_stats{$ecod_domain_id}{struct_seqid_aref};
			 $hit_asym_id			= $$rep_stats{$ecod_domain_id}{asym_id};
		}elsif(exists $$rep_stats{$hit_uid}) { 
			 $hit_seqid_aref 		= $$rep_stats{$hit_uid}{seqid_aref};
			 $hit_struct_seqid_aref 	= $$rep_stats{$hit_uid}{struct_seqid_aref};
			 $hit_asym_id			= $$rep_stats{$hit_uid}{asym_id};
		}else{
			die "ERROR! Could not deduce rep_stats structure\n";
		}
		if (!$hit_seqid_aref) { 
			print "WARNING! $hit_pdb $hit_chain $ecod_domain_id obsolete, skipping...\n";
			next;
		}
		my (@matrix, $z, $rmsd, $id);
		my ($mob_range, $ref_range);
		my (@full_mob_range, @full_ref_range);
		my $al_count = 0;

		my $summ_exists = 0;

		my $nmatch = 0;
		my $sum_scores1 = 0;
		my $sum_scores2 = 0;
		my $sum_scores12 = 0;

		my %blosum_ref_count;
		my %blosum_mob_count;
		for (my $i = 0; $i < scalar(@aa); $i++) { 
			$blosum_ref_count{$aa[$i]} = 0;
			$blosum_mob_count{$aa[$i]} = 0;
		}

		while (my $ln = <IN>) { 
			if ($ln =~ /-protein/) { 
				$al_count++;
				if ($al_count == 6) { last } 
			}

			if ($ln =~ /-matrix/) { 
				$ln =~ s/"//g;
				my @F = split(/\s+/, $ln);
				$F[3]	=~ /U\(([1-3])\,\.\)/ or die "ERROR $sub: matrix regexp failure!\n";
				my $mat_i = $1;
				for my $mat_j (1, 2, 3, 4) { 
					$matrix[$mat_i][$mat_j] = $F[($mat_j+3)];
				}
			}
			if ($ln =~ /-transrotate/) { 
				my @F = split(/\s+/, $ln);
				$matrix[4][1]	= 0;
				$matrix[4][2]	= 1;
				$matrix[4][3]	= 2;
				$matrix[4][4]	= 3;
			}
			if ($ln =~ /-ranges/) { 
				#These are pdb residue numbers, not pos inx
				my $ref_start 	= substr($ln, 25, 4);
				$ref_start 	=~ s/\s+//g;
				my $ref_end	= substr($ln, 31, 4);
				$ref_end 	=~ s/\s+//g;

				my $mob_start	= substr($ln, 40, 4);
				$mob_start	=~ s/\s+//g;
				my $mob_end	= substr($ln, 46, 4);
				$mob_end	=~ s/\s+//g;

				$ref_range	= "$ref_start-$ref_end";
				$mob_range	= "$mob_start-$mob_end";
				#AH HA - query_reg is already in seq_id reference.
				my $ref_range_aref = struct_region(range_expand($ref_range), $$query_href{struct_seqid_aref});
				my $mob_range_aref = struct_region(range_expand($mob_range), $hit_struct_seqid_aref);

				push (@full_ref_range, @$ref_range_aref);
				push (@full_mob_range, @$mob_range_aref);

				#Calculate BLOSUM score from bits matrix 

				my $ref_seq_href = seq_grab($ref_range_aref, $$query_href{seq});

				my $mob_seq_href;
				if( exists $$rep_stats{$hit_uid}{seq}){
					$mob_seq_href = seq_grab($mob_range_aref, $$rep_stats{$hit_uid}{seq});
					#printf "keys mobseq %i\n", keys %$mob_seq_href;
				}else{
					die "??";
				}
				die if keys %$mob_seq_href == 0;


				for (my $i = 0; $i < scalar(@$ref_range_aref); $i++) {  

					my $ref_aa = $$ref_seq_href{$$ref_range_aref[$i]};
					my $mob_aa = $$mob_seq_href{$$mob_range_aref[$i]};
					if (!$mob_aa or !$ref_aa) { 
						print "Failed hp: $hit_pdb ma $mob_aa mra $$mob_range_aref[$i] ra $ref_aa rra $$ref_range_aref[$i]\n";
						next;
					}

					$nmatch++;
					$blosum_ref_count{$ref_aa}++;
					$blosum_mob_count{$mob_aa}++;
					if (!$blosum_score_matrix{$ref_aa}{$ref_aa}) { 
						#die "ERROR! No BLOSUM score for $ref_aa $hit_pdb $$ref_range_aref[$i]\n";
						#print "WARNING! No BLOSUM score for $ref_aa $hit_pdb $$ref_range_aref[$i]\n";
						next;
					}
					if (!$blosum_score_matrix{$mob_aa}{$mob_aa}) { 
						#die "ERROR! No BLOSUM score for $mob_aa $hit_pdb $$mob_range_aref[$i]\n";
						#print "WARNING! No BLOSUM score for $mob_aa $hit_pdb $$mob_range_aref[$i]\n";
						next;
					}
					if (!$blosum_score_matrix{$mob_aa}{$ref_aa}) { 
						#die "ERROR! No BlOSUM score for $ref_aa/$mob_aa\n";
						#print "WARNING! No BLOSUM score for $ref_aa $mob_aa\n";
						next;
					}
					$sum_scores1	+= $blosum_score_matrix{$ref_aa}{$ref_aa};
					$sum_scores2	+= $blosum_score_matrix{$mob_aa}{$mob_aa};
					$sum_scores12	+= $blosum_score_matrix{$ref_aa}{$mob_aa};

					#printf "%i: %s %s\n", $nmatch, $ref_aa, $mob_aa;
				}

				#END SCORE
				if ($DEBUG > 1) { 
					print "DEBUG: $sub: $dali_file $ref_range $mob_range\n";
				}
			}
				

			if ($ln =~ /-summar/) { 
				$summ_exists = 1;
				$ln =~ s/^\s+//;
				my @F = split(/\s+/, $ln);
				$z 	= $F[3];
				if ($z < $ABSOLUTE_MINIMUM) {next } 
				$rmsd	= $F[4];
				$id	= $F[7];
				#if ($DEBUG > 1) { 
				#	print "DEBUG $sub: $z $rmsd $id $dali_file \n";
				#}

				if (scalar(@full_mob_range) == 0 || scalar(@full_ref_range) == 0 ) {  next } 
				$mob_range = rangify(@full_mob_range);
				$ref_range = rangify(@full_ref_range);

				#if (!$ecod_rep_range{$ecod_domain_id}) { next } 
				my ($domain_range_aref, $domain_chain_aref);
				if (exists $$rep_stats{$ecod_domain_id}{range}) { 
					($domain_range_aref, $domain_chain_aref) = multi_chain_range_expand($$rep_stats{$ecod_domain_id}{range});
				}elsif (exists $$rep_stats{$hit_uid}{range}) { 
					($domain_range_aref, $domain_chain_aref) = multi_chain_range_expand($$rep_stats{$hit_uid}{range});
				}else{
					next;
				}
				my $domain_res = scalar(@$domain_range_aref);
				my $mob_res = scalar(@{ range_expand($mob_range) });

				my $coverage = sprintf "%.2f", $mob_res/$domain_res;

				#BLOSUM FINAL

				my $rand_sum = 0.0;
				for (my $i = 0; $i < scalar(@aa); $i++) { 
					for (my $j = 0; $j < scalar(@aa); $j++) { 
						$rand_sum += ($blosum_ref_count{$aa[$i]}/$nmatch) * ( $blosum_mob_count{$aa[$j]}/$nmatch)*$blosum_score_matrix{$aa[$i]}{$aa[$j]};
					}
				}

				$rand_sum *= $nmatch;

				my $norm_blosum_hit_score = 0;
				if ( ($sum_scores1+$sum_scores2) - $rand_sum != 0) { 
					$norm_blosum_hit_score = ($sum_scores12 - $rand_sum) / (0.5 * ($sum_scores1 + $sum_scores2) - $rand_sum);
				}else{
					$norm_blosum_hit_score = 'NaN';
				}

				if ($DEBUG == 1) { 
					print "DEBUG $sub; $j $z $rmsd $id $norm_blosum_hit_score $dali_file\n";
				}

				$hits[$j]{mob_range} = $mob_range;
				$hits[$j]{ref_range} = $ref_range;
				$hits[$j]{z_score} 	= $z;
				$hits[$j]{rmsd}		= $rmsd;
				$hits[$j]{id}		= $id;
				$hits[$j]{dali_fn} 	= $dali_file;
				$hits[$j]{coverage} 	= $coverage;
				$hits[$j]{ecod_domain_id}	= $ecod_domain_id;
				$hits[$j]{blosum_dali_score}	= $norm_blosum_hit_score;
				if ($hit_uid) { 
					$hits[$j]{uid} =  $hit_uid;
				}
				$j++;

				undef(@full_mob_range);
				undef(@full_ref_range);

			}

		}
		close IN;


#		if (scalar(@full_mob_range) == 0 || scalar(@full_ref_range) == 0 ) { next } 
#		$mob_range = rangify(@full_mob_range);
#		$ref_range = rangify(@full_ref_range);

#		if ($summ_exists > 0) { 
#			if ($DEBUG == 1) { 
#				print "DEBUG $sub: $j $z $rmsd $id $dali_file \n";
#			}
#
#			$dali_file =~ /(e\w{4}.\d+)/;
#			my $ecod_domain_id = $1;
#
#			my ($domain_range_aref, $domain_chain_aref) = multi_chain_range_expand($ecod_rep_range{$ecod_domain_id});
#
#			my $domain_res = scalar(@$domain_range_aref);
#			my $mob_res = scalar(@{ range_expand($mob_range) });
#
#
#			my $coverage = sprintf "%.2f", $mob_res/$domain_res;
#			
#
#			$hits[$j]{mob_range} = $mob_range;
#			$hits[$j]{ref_range} = $ref_range;
#			$hits[$j]{z_score} 	= $z;
#			$hits[$j]{rmsd}		= $rmsd;
#			$hits[$j]{id}		= $id;
#			$hits[$j]{dali_fn} 	= $dali_file;
#			$hits[$j]{coverage} 	= $coverage;
#			$hits[$j]{ecod_domain_id}	= $ecod_domain_id;
#			$j++;
#
#		}
	}

	@hits = sort {$b->{z_score} <=> $a->{z_score}} @hits;

	for (my $i = 0; $i < scalar(@hits); $i++) { 

		my $z	 = $hits[$i]{z_score};
		my $rmsd = $hits[$i]{rmsd};
		my $id	 = $hits[$i]{id};
		my $blsm_norm = $hits[$i]{blosum_dali_score};

		my $mob_range = $hits[$i]{mob_range};
		my $ref_range = $hits[$i]{ref_range};

		my $dali_file = $hits[$i]{dali_fn};
		my $ecod_domain_id = $hits[$i]{ecod_domain_id};
		my $coverage = $hits[$i]{coverage};

		if ($DEBUG) { 
			print "DEBUG $sub: $i $id $dali_file $z $rmsd $id $mob_range $ref_range\n";
		}
		my $dali_hit_node	= $dali_summ_xml_doc->createElement('hit');

		my $dali_file_node	= $dali_summ_xml_doc->createElement('dali_file');
		$dali_file_node->appendTextNode($dali_file);
		$dali_hit_node->appendChild($dali_file_node);

		$dali_hit_node->setAttribute('z-score', $z);
		$dali_hit_node->setAttribute('RMSD', $rmsd);
		$dali_hit_node->setAttribute('identity', $id);
		$dali_hit_node->setAttribute('ecod_domain_id', $ecod_domain_id);
		if ($hits[$i]{uid}) { 
			$dali_hit_node->setAttribute('uid', $hits[$i]{uid});
		}
		$dali_hit_node->setAttribute('coverage', $coverage);
		$dali_hit_node->setAttribute('blosum_norm_dali_score', $blsm_norm);

		#This should be seqid range if your reference pdb set is seqid indexed
		my $dali_query_range_node	= $dali_summ_xml_doc->createElement('query_reg');
		#$dali_query_range_node->appendTextNode($mob_range);
		$dali_query_range_node->appendTextNode($ref_range);
		$dali_hit_node->appendChild($dali_query_range_node);

		my $dali_hit_range_node		= $dali_summ_xml_doc->createElement('hit_reg');
		#$dali_hit_range_node->appendTextNode($ref_range);
		$dali_hit_range_node->appendTextNode($mob_range);
		$dali_hit_node->appendChild($dali_hit_range_node);

		my $pdb_range = pdb_rangify(range_expand($ref_range), $$query_href{pdbnum_aref});

		my $dali_query_pdb_range_node = $dali_summ_xml_doc->createElement('query_pdb_range');
		$dali_query_pdb_range_node->appendTextNode($pdb_range);
		$dali_hit_node->appendChild($dali_query_pdb_range_node);

		$dali_hit_list_node->appendChild($dali_hit_node);

	
	}
	my $doc_string = $dali_summ_xml_doc->toString(1);
	open (OUT, ">$dali_summ_fn") or die "ERORR! Could not open $dali_summ_fn\n";
	print OUT $doc_string;
	close OUT;
}

sub blosum62_href_load { 

	my $sub = 'blosum62_href_load';

	my @blosum_score_matrix;
	my %return;

	my @cheek_aa = qw(ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL);

	@blosum_score_matrix = (
	["1.9646"],
	["-0.7068","2.7367"],
	["-0.7654","-0.2199","2.8266"],
	["-0.8767","-0.8029","0.6358","2.8871"],
	["-0.2043","-1.6946","-1.3299","-1.7300","4.2911"],
	["-0.4020","0.4914","0.0008","-0.1567","-1.4509","2.6426"],
	["-0.4319","-0.0577","-0.1340","0.7552","-1.8062","0.9273","2.4514"],
	["0.0798","-1.1521","-0.2114","-0.6568","-1.2502","-0.8926","-1.0551","2.7816"],
	["-0.8126","-0.1249","0.2892","-0.5595","-1.4939","0.2240","-0.0588","-1.0204","3.7555"],
	["-0.6609","-1.4951","-1.6085","-1.5606","-0.6138","-1.3848","-1.5972","-1.8624","-1.6158","1.9993"],
	["-0.7323","-1.0773","-1.6895","-1.8028","-0.6387","-1.0670","-1.4232","-1.8135","-1.3934","0.7608","1.9247"],
	["-0.3670","1.0544","-0.0895","-0.3509","-1.5182","0.6363","0.3877","-0.7640","-0.3605","-1.3351","-1.2234","2.2523"],
	["-0.4676","-0.6836","-1.0754","-1.5293","-0.7099","-0.2105","-0.9990","-1.3383","-0.7756","0.5634","0.9959","-0.6774","2.6963"],
	["-1.1050","-1.3932","-1.4970","-1.7419","-1.1877","-1.5822","-1.5962","-1.5537","-0.6171","-0.0804","0.2074","-1.5393","0.0063","3.0230"],
	["-0.4071","-1.0543","-1.0002","-0.7401","-1.3976","-0.6410","-0.5581","-1.0668","-1.0805","-1.3783","-1.4300","-0.5068","-1.2382","-1.7986","3.6823"],
	["0.5579","-0.3824","0.3005","-0.1305","-0.4375","-0.0506","-0.0735","-0.1462","-0.4408","-1.1741","-1.2213","-0.1017","-0.7404","-1.1845","-0.4045","1.9422"],
	["-0.0227","-0.5612","-0.0230","-0.5254","-0.4333","-0.3377","-0.4316","-0.7877","-0.8429","-0.3588","-0.5987","-0.3348","-0.3331","-1.0538","-0.5376","0.6906","2.2727"],
	["-1.2634","-1.3397","-1.8480","-2.1072","-1.1521","-0.9732","-1.4177","-1.2457","-1.1711","-1.2903","-0.8159","-1.4782","-0.7124","0.4588","-1.8271","-1.3759","-1.2145","5.2520"],
	["-0.8820","-0.8469","-1.0409","-1.5325","-1.2036","-0.7105","-1.0102","-1.5199","0.8463","-0.6657","-0.5310","-0.9100","-0.4974","1.4696","-1.4599","-0.8429","-0.8030","1.0771","3.2975"],
	["-0.0947","-1.2513","-1.4382","-1.5713","-0.4038","-1.0992","-1.2211","-1.5694","-1.5587","1.2735","0.3942","-1.1312","0.3436","-0.4245","-1.1744","-0.8231","-0.0278","-1.4171","-0.6038","1.8845"]
	);


	for (my $i = 0; $i < scalar(@blosum_score_matrix); $i++) { 
		for (my $j = 0; $j < scalar(@blosum_score_matrix); $j++) { 
			if (!$blosum_score_matrix[$i][$j]) { next } 

			my $aa1 = $cheek_aa[$i];
			my $aa2 = $cheek_aa[$j];

			$return{$aa1}{$aa2} = $blosum_score_matrix[$i][$j];
			$return{$aa2}{$aa1} = $blosum_score_matrix[$i][$j];
		}
	}

	return %return;
}

sub seq_grab { 
	my $sub = 'seq_grab';

	my ($range_aref, $seq_href) = @_;

	my %sub_seq;
	foreach my $seqid (@$range_aref) { 
		if ($$seq_href{$seqid}) { 
			$sub_seq{$seqid} = $$seq_href{$seqid};
		}
	}
	return \%sub_seq;
}

sub isect { 
	my $sub = 'isect';

	my ($aref1, $aref2) = @_;

	my %seen;
	foreach my $i (@$aref1) { 
		$seen{$i}++;
	}
	my @isect;
	foreach my $i (@$aref2) { 
		if ($seen{$i}) { 
			push (@isect, $i);
		}
	}
	return \@isect;
}
	

