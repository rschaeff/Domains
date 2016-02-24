package Domains::Partition;
require Exporter;

use warnings;
use strict;

use Data::Dumper;

use Date::Format;
use Date::Parse;

use Storable;
use Cwd;

use Carp;
use File::Copy;
use File::Remove qw(remove);

use Domains::Range;
use Domains::PDBml;
use Domains::Dali;
use Domains::SGE;
use XML::Grishin;

use ECOD::Reference;

use Getopt::Long;

umask(0002);
my $UID =  1219; #rschaeff
my $GID = 902; #ecod


our @ISA = qw(Exporter);
our @EXPORT = (
		"&domain_partition",
		#"&domain_partition_asm",
		"&dali_summary",
		"&reference_cache_load",
		"&generate_singleton_run_list",
		"&generate_domain_summary_files",
		"&generate_struct_search_jobs",
		"&job_list_generate_fasta",
		"&run_list_maintain",
		"&job_list_maintain",
		"&jobify_dali",
		"&immediate_dali",
		'$pdb_dir',
		'$pdbml_dir',
		'$status_dir',
		'&buildali_hhblits',
		'&blastp_job_file',
		'&chblastp_job_file',
		'&hh_gen',
		'&hh_run',
		'&hh_parse',
		'&register_repair',
		'&reference_library_load',
		'&read_self_comps',
		'&struct_search_dali_query_gen',
		'&generate_query_pdb',
		'&bsumm',
		'&query_glob',
		'&domains',
		'&optimize',
		'&peptide_filter',
		'&peptide_collate',
		'&load_partition_vars',
		'&find_hh_dali_domains',
		'$DOMAIN_PREFIX',
		'$DOMAIN_DATA_DIR',
		'$GAP_TOL',
		);



#These shouldn't be defined in a package
my $DEBUG = 1;
my $FORCE_OVERWRITE = 0;
my $FORCE_HH_OVERWRITE = 0;

#generate_struct_search
my $GENERATE_NEW_FILES = 1;
my $PSIPRED_EXE 	= '/usr1/psipred/bin';
my $PSIPRED_DATA	= '/usr1/psipred/data';
my $DALI_EXE 		= '/home/rschaeff/bin/DaliLite_nohtml';
my $HHBLITS_EXE		= '/home/rschaeff/hhsuite-2.0.15-linux-x86_64/bin/hhblits';
my $HHMAKE_EXE		= '/home/rschaeff/hhsuite-2.0.15-linux-x86_64/bin/hhmake';
my $HHSEARCH_EXE	= '/home/rschaeff/hhsuite-2.0.15-linux-x86_64/bin/hhsearch';
my $HH_LIB		= '/home/rschaeff/hhsuite-2.0.15-linux-x86_64/lib/hh';
my $SINGLE_BOUNDARY_OPTIMIZE = '/data/ecod/weekly_updates/weeks/bin/single_boundary_optimize.pl';
my $SINGLE_PEPTIDE_FILTER    = '/data/ecod/weekly_updates/weeks/bin/single_xml_peptide_prefilter.pl';

#my $NR20_TMP_DB		= '/home/rschaeff_1/side/wtf_hh/nr20_12Aug11';
my $NR20_TMP_DB			= '/local_scratch/rschaeff/nr20_12Aug11';
my $PDB_CACHE			= '/home/rschaeff_1/pdb_junk_bin';

my $PERFORM_NREP_CONCISE = 1 ;
my $POOR_ONLY = 0;
my $NEW_COVERAGE_THRESHOLD = 0.9;
my $OLD_COVERAGE_THRESHOLD = 0.05;
my $DALI_SIGNFICANCE_THRESHOLD = 5.0;
my $HIT_COVERAGE_THRESHOLD = 0.7;
my $NO_PSIBLAST = 1; #Is PSIBLAST providing ANY value over HHSEARCH? It is unclear.

#job_list_asess
my $INCLUSION_THRESHOLD = 0.7;
my $FILTER_VERSION = 6; #peptide filter
my $CC_FILTER_VERSION = 1;
my $SCREEN_REPS = 1;	#Screen for rep95 nodes 
my $SCREEN_MC_DOM = 1; #Screen for MC domain components

our $DOMAIN_PREFIX;
our $GAP_TOL;


my $DOMAIN_DATA_DIR = '/data/ecod/domain_data';
if (!-d $DOMAIN_DATA_DIR) { 
	die "ERROR! domain data dir $DOMAIN_DATA_DIR not found\n";
}
my $pdb_dir = "/usr2/pdb/data/structures/divided/pdb/";
if (!-d $pdb_dir) { 
	die "ERROR! pdb dir $pdb_dir not found\n";
}
my $pdbml_dir = "/usr2/pdb/data/structures/divided/XML-noatom/";
if (!-d $pdbml_dir) { 
	die "ERROR! pdbml dir $pdbml_dir not found\n";
}

my $status_dir = '/usr2/pdb/data/status';
if (!-d $status_dir) { 
	die "ERROR! status directory $status_dir not found\n";
}

sub load_partition_vars { 
	$DOMAIN_PREFIX = 'domains_v13';
	$GAP_TOL = 20;
	$DOMAIN_DATA_DIR;
}

sub struct_search_dali_query_gen { 
	my $sub = 'struct_search_dali_query_gen';

	my ($job_list_xml_fn, $ref) = @_;
	load_references();
	if (!$ref) {
		$ref = $LATEST_REFERENCE;
	}
	my $ecod_xml_doc 	= xml_open($REF_XML{$ref});
	my $job_xml_doc 	= xml_open($job_list_xml_fn);
	my $reference 		= $ref;

	#rep find 
	my $rep_domains = $ecod_xml_doc->findnodes(qq{//domain[\@manual_rep='true']});
	if ($DEBUG) { 
		printf "DEBUG: Found %i rep domain in ref $reference\n", $rep_domains->size();
	}

	my @ecod_reps;
	my %uid_lookup;
	foreach my $rep_domain_node ($rep_domains->get_nodelist() ) { 

		my ($uid, $ecod_domain_id) = get_ids($rep_domain_node);
		my $short_uid = substr($uid,2,5);
		 
		next if isObsolete($rep_domain_node);
	
		if ($rep_domain_node->findvalue('structure/@structure_obsolete') eq 'true') { next } 

		if ($ecod_domain_id =~ /e(\w(\w{2})\w)/) { 
			my $ecod_pdb 	= $1;
			my $ecod_two 	= $2;
			my $path 	= "$DOMAIN_DATA_DIR/$short_uid/$uid/$uid.seqres.pdb";
			if (!-f $path) { 
				print "WARNING! NO CLEAN DOMAIN PDB FILE FOR $ecod_domain_id ($path)\n";
			}else{
				push (@ecod_reps, $ecod_domain_id);
				$uid_lookup{$ecod_domain_id} = $uid;
			}
		}else{
			print "WARNING! malformed domain_id $ecod_domain_id, skipping\n";
		}
	}
	my $job_nodes = find_job_nodes($job_xml_doc);

	if ($DEBUG) { 
		printf "DEBUG: Found %i struct_search jobs\n", $job_nodes->size();
	}

	my ($job_dump_dir, $job_list_dir) = get_job_dirs_from_job_xml($job_xml_doc);

	foreach my $job_node ($job_nodes->get_nodelist()) { 

		my $job_id	= $job_node->findvalue('@id');
		my $reference	= $job_node->findvalue('reference');
		my ($pdb, $chain) = get_pdb_chain($job_node);

		if ($reference ne $ref) { 
			print "WARNING! job reference $reference does not match input refrences $LATEST_REFERENCE\n";
		}
		my $mode	= $job_node->findvalue('mode');

		if ($mode ne 'struct_search') { next } 	

		my $pdb_chain = $pdb . "_" . $chain;
		print "pc: $pdb_chain\n";
		#Abandon if DALI summ or domains exists, use only for run completion not overwrite

		my $recurse = 0;
		my $recurse_range;
		if ($job_node->findvalue('@type') eq 'recurse'){ 
			my $seqid_range = $job_node->findvalue('seqid_range');
			$recurse_range = $seqid_range;
			$seqid_range =~ s/\,/_/g;
			$pdb_chain .=  "_$seqid_range";
			$recurse++;
		}
		
		my $dali_summ_fn = "$job_dump_dir/$pdb_chain/$pdb_chain.$reference.dali_summ.xml"; 
		if (-f $dali_summ_fn) { 
			print "WARNING! Dali summary file found for $pdb_chain, skipping...\n";
			next;
		}

		#my $peptide_filter;

		my $query_pdb_chain_fn = "$job_dump_dir/$pdb_chain/$pdb_chain.pdb";

		if (!-f $query_pdb_chain_fn || $FORCE_OVERWRITE) { #generate query pdb 
			print "gen: $query_pdb_chain_fn\n";
			generate_query_pdb($query_pdb_chain_fn, $pdb, $chain);
		}else{
			if ($DEBUG) { 
				print "DEBUG: $query_pdb_chain_fn exists!\n";
			}
		}

		#DaliLite

		#Generate query brk file 
		my $DaliLite_exe  = $DALI_EXE;
		if (!-f $DaliLite_exe) { 
			die "ERROR! Could not find DaliLite executable $DaliLite_exe\n";
		}

		#Goofball DALI
		my $dali_tmp_dir = "$job_dump_dir/$pdb_chain/dali_tmp";
		if (!-d $dali_tmp_dir) { 
			mkdir($dali_tmp_dir);
		}

		foreach my $ecod_rep (@ecod_reps) { 
			my $uid = $uid_lookup{$ecod_rep};
			my $short_uid = substr($uid,2,5);

			if ($DEBUG > 1) { 
				print "DALI jobify $query_pdb_chain_fn $DOMAIN_DATA_DIR/$short_uid/$uid/$uid.seqres.pdb\n";
			}
			my $dali_output_fn = "q.$pdb_chain.$ecod_rep.$reference.dali";

			jobify_dali($dali_tmp_dir, "../../$pdb_chain.pdb" , "$DOMAIN_DATA_DIR/$short_uid/$uid/$uid.seqres.pdb", $dali_output_fn);
		}

	}
}
sub jobify_dali { 
	my $sub = 'jobify_dali';
	
	my ($dali_dir, $pdb1_fn, $pdb2_fn, $dali_output_fn) = @_;

	if (!-f $DALI_EXE) { 
		die "ERROR! $sub: DALI exe $DALI_EXE not found\n";
	}
	my $job_file = $dali_output_fn;
	$job_file =~ s/dali/job/;

	$pdb1_fn =~ /\/([\d\w\_\-\.]+)\.pdb/;
	my $name1 = $1;

	$pdb2_fn =~ /\/([\d\w\.]+)\.pdb/;
	my $name2 = $1;
	my $mixup = $name1 . "_" . $name2;

	open (my $fh, ">", "$dali_dir/$job_file") or die "ERROR! $sub: Could not open $dali_dir/$job_file for writing:$!\n";

	print $fh "#!/bin/bash\n";
	print $fh "#\$ -cwd \n";
	print $fh "#\$ -j y\n";
	print $fh "#\$ -M dustin.schaeffer\@gmail.com\n";
	print $fh "#\$ -v LD_LIBRARY_PATH\n";
	print $fh "mkdir $dali_dir/$mixup\n"; 
	print $fh "cd $dali_dir/$mixup\n";
	print $fh "$DALI_EXE -pairwise $pdb1_fn $pdb2_fn > /dev/null\n";
	my $DALI_OUTPUT = 'mol1?.result';
	print $fh "mv $DALI_OUTPUT ../$dali_output_fn\n";
	print $fh "rm -f dali.lock\n";
	print $fh "cd ..\n";
	print $fh "rm -f $dali_dir/$mixup/*\n";
	print $fh "rmdir $dali_dir/$mixup\n";
	close OUT;

	return "$dali_dir/$job_file";
}

sub immediate_dali { 
	my $sub = 'immediate_dali';

	my ($dali_dir, $pdb1_fn, $pdb2_fn, $dali_output_fn) = @_;

	if (!-f $DALI_EXE) { 
		die "ERROR! $sub: DALI exe $DALI_EXE not found\n";
	}
	my $job_file = $dali_output_fn;
	$job_file =~ s/dali/job/;

	$pdb1_fn =~ /\/([\w\_\-]+)\.pdb/;
	my $name1 = $1;

	$pdb2_fn =~ /\/([\w\.]+)\.pdb/;
	my $name2 = $1;
	my $mixup = $name1 . "_" . $name2;

	my $curdir = getcwd();
	mkdir("$dali_dir/$mixup");
	chdir("$dali_dir/$mixup");
	system("$DALI_EXE -pairwise $pdb1_fn $pdb2_fn > /dev/null\n");
	my $DALI_OUTPUT = 'mol1?.result';
	#move($DALI_OUTPUT, "../$dali_output_fn");
	system("mv $DALI_OUTPUT ../$dali_output_fn");
	remove("dali.lock");
	chdir("..");
	remove("$dali_dir/$mixup/*");
	rmdir("$dali_dir/$mixup");
	chdir($curdir);

}

sub existsMode { 
	return $_[0]->exists(qq{run_list_job_xml_file[\@mode="$_[1]"]}) ? 1 : 0  
}
sub generate_query_pdb { 
	my ($out_fn, $pdb, $chain, $input_seqid_aref) = @_;


	my ($seqid_aref, $struct_seqid_aref, $pdbnum_aref, $asym_id) = pdbml_seq_parse($pdb, $chain);
	if (!$seqid_aref) { 
		print "WARNING! PDB $pdb deprecated, skipping...\n";
		return 0;
	}
	my $coord_aref;
	if ($input_seqid_aref) { 
		($coord_aref) = pdbml_coord_parse_pull($pdb, $input_seqid_aref, $asym_id);
	}else{
		($coord_aref) = pdbml_coord_parse_pull($pdb, $seqid_aref, $asym_id);
	}

	#is the output coord ref sane?
	my $coord_size = scalar(@$coord_aref); #How large is it? Needs to exceed 30 residues for Dali to work
	my $ca_only = 0;
	my %atoms;
	for (my $i = 0; $i < scalar(@$coord_aref); $i++) { 
		my $atom = $$coord_aref[$i]{label_atom_id};
		$atoms{$atom}++;
	}

	if ($DEBUG) { 
		printf "$pdb, $chain, keys %i\n", scalar(keys %atoms);
	}
	if (scalar(keys %atoms) == 1) { 
		 print "$pdb, $chain CAONLY\n";
		 return 0;
	}

	pdbml2pdb_cartn($out_fn, $coord_aref, "seq");


}
sub generate_struct_search_jobs { 
	my $sub = 'generate_struct_search_jobs';

	my $job_assess 		= 'seq_iter';
	my $job_generate 	= 'struct_search';

	load_partition_vars();

	my ($run_list_xml_doc) = @_;

	DP_NODE:
	foreach my $dp_node (find_domain_parse_nodes($run_list_xml_doc)) { 
	
		if (existsMode($dp_node, $job_generate)) { 
			next;
		}

		my $run_id		= $dp_node->findvalue('@run_id');
		my $week_label		= $dp_node->findvalue('week_label');
		my $run_list_dir 	= $dp_node->findvalue('run_list_dir');
		
		if ($DEBUG)  { 
			print "DEBUG: main: $run_id $week_label $run_list_dir\n";
		}


		if (existsMode($dp_node, $job_assess) && !existsMode($dp_node, $job_generate)) { 

			my $assess_file = get_run_list_job_xml_file($dp_node, $job_assess); #Shouldn't this be $job_assess
			my $generate_file;

			#my $generate_file_XPath = qq{run_list_job_xml_file[\@mode="$job_generate"]};
			
			if ($assess_file =~ /job/) { 
				$generate_file = $assess_file;
				$generate_file =~ s/job/$job_generate.job/;
			}else{
				die "REGEXP substitute fails on $assess_file for generate_file generation\n";
			}

			my ($generate_xml_doc, $generate_root_node) = xml_create('job_set_top');

			create_job_dir_nodes_for_job_xml($generate_xml_doc, $generate_root_node, $run_list_dir);
				
			xml_write($generate_xml_doc, "$run_list_dir/$generate_file");

			my $dp_job_xml_file_node = $run_list_xml_doc->createElement('run_list_job_xml_file');
			$dp_job_xml_file_node->setAttribute('mode', $job_generate);
			$dp_job_xml_file_node->appendTextNode($generate_file);

			$dp_node->appendChild($dp_job_xml_file_node);
			

			job_list_assess($run_list_dir, $assess_file, $generate_file, $job_assess, $job_generate, $DOMAIN_PREFIX);

		}else{
			print "WARNING! No $job_assess node found for $run_id $week_label\n";
		}
	}
}
sub job_list_assess { 
	my $sub = 'job_list_assess';

	my ($run_list_dir, $assess_file, $generate_file, $job_assess, $job_generate, $domain_prefix) = @_;

	if ($DEBUG) { 
		print "DEBUG $sub: top $assess_file $generate_file\n";
	}

	if (! -f "$run_list_dir/$assess_file") { 
		die "ERROR! $sub: assess file $assess_file not found\n";
	}

	if (! -f "$run_list_dir/$generate_file") { 
		die "ERROR! $sub: generate file $generate_file not found\n";
	}

	my $assess_xml_doc 	= xml_open("$run_list_dir/$assess_file");
	my $generate_xml_doc 	= xml_open("$run_list_dir/$generate_file");

	#process assess jobs
	my %known;
	foreach my $as_node (find_job_nodes($generate_xml_doc)) { 
		my $job_id 	= $as_node->findvalue('@id');

		my ($query_pdb, $query_chain) = get_pdb_chain($as_node);
		my $reference 	= $as_node->findvalue('reference');
		my $mode		= $as_node->findvalue('mode');

		$known{$query_pdb}{$query_chain}{$reference}{$mode}++;
	}

	my @gen;
	my $gen_id = 1;
	my $gen_i = 0;
	my $reps;

	if ($assess_xml_doc->exists('//@rep95')) { 
		$reps = 1;
	}
	foreach my $aj_node (find_job_nodes($assess_xml_doc)) { 

		my $job_id 	= $aj_node->findvalue('@id');
		my $rep95	= $aj_node->findvalue('@rep95');
		my $mc_dom	= $aj_node->findvalue('@sw_mc_complex');

		my ($query_pdb, $query_chain) = get_pdb_chain($aj_node);	

		my $reference 	= $aj_node->findvalue('reference');
		my $mode	= $aj_node->findvalue('mode');

		if ($SCREEN_REPS && !$rep95 && $reps) { 
			#print "WARNING! Not rep $query_pdb $query_chain\n";
			next;
		} 
		if ($SCREEN_MC_DOM && $mc_dom) { 
			#print "WARNING! MC_DOM componenent $query_pdb $query_chain\n";
			next;
		}


		if (isPeptide($aj_node)) { 
			if ($DEBUG) { 
				print "DEBUG: $sub: $job_id $query_pdb $query_chain peptide filter apply, skipping...\n";
			}
			next;
		}else{
			#print "WARNING! $sub: No peptide filter node found for $job_id $query_pdb $query_chain... proceeding anyways\n";
		}
		
		if (isCoiledCoil($aj_node)) { 
			if ($DEBUG) { 
				print "DEBUG: $sub $job_id $query_pdb $query_chain coiled coil filter apply, skipping...\n";
			}
			next;
		}else{
			#print "WARNING! $sub: No coiled coil filter node for $job_id $query_pdb $query_chain...\n";
		}
	
		if ($mode eq $job_assess) { 
			if (job_assess($run_list_dir, $query_pdb, $query_chain, $domain_prefix, $reference, $job_assess)) { 
				#SUCCESS
			}elsif(!$known{$query_pdb}{$query_chain}{$reference}{$mode}) { 
				#Generate struct_search node
				$gen[$gen_i]{id} 	= $gen_id;
				$gen[$gen_i]{query_pdb}	= $query_pdb;
				$gen[$gen_i]{query_chain}	= $query_chain;
				$gen[$gen_i]{reference}	= $reference;
				$gen[$gen_i]{mode} = $job_generate;
				
				$gen_i++;
				$gen_id++;
			}else{
				die "So confused\n";
			}
		}
	}

	my $as_job_list_node = $generate_xml_doc->findnodes('//job_list')->get_node(1);

	for (my $i = 0; $i < $gen_i; $i++) { 

		my $as_node	= create_job($generate_xml_doc, \%{$gen[$i]});
		if ($SCREEN_REPS) { 
			$as_node->setAttribute('rep95', 'true');
		}

		$as_job_list_node->appendChild($as_node);
	}

	xml_write($generate_xml_doc, "$run_list_dir/$generate_file");

}
sub job_assess { 
	my $sub = 'job_assess'; 
	
	my ($run_list_dir, $query_pdb, $query_chain, $domain_prefix, $reference, $job_assess) = @_;

	if ($DEBUG) { 
		print "DEBUG: $sub: top $run_list_dir, $query_pdb, $query_chain, $domain_prefix, $reference, $job_assess\n";
	}

	my $pdb_chain = $query_pdb . "_" . $query_chain;

	my $fp = "$run_list_dir/ecod_dump/$pdb_chain/$domain_prefix.$pdb_chain.$reference.xml";
	if (-f $fp) { 
		#Found file!

		my $domain_xml_doc = xml_open($fp);
		
		my $chain_domain_coverage = $domain_xml_doc->findvalue('//chain_domain_coverage');
		if ($domain_xml_doc->exists('//optimized_chain_domain_coverage')) { 
			my $unused_residues = $domain_xml_doc->findvalue('//optimized_chain_domain_coverage/@unused_res');

			if ($unused_residues > 0) { 
				return 0
			}else{
				return 1
			}
		}elsif ($chain_domain_coverage == 0) {
			#Empty domain file case
			return 0;
		}elsif ($chain_domain_coverage == 1) { 
			#Full assignment case
			return 1;
		}else{
			print "WARNING! $sub: Could not find optimized ranges in $fp\n";
			return 1;
		}

	}else{
		print "WARNING! $sub: Could not find $fp\n";
		return 1;
	}
}
sub generate_singleton_run_list { 
	my $sub = 'generate_singleton_run_list';

	my ($job_list_file, $run_label) = @_;

	my $dir;
	if (!$run_label) { 
		if ($job_list_file =~ /(.*)\/(repair.develop\d+)/) { 
			$dir 	= $1;
			$run_label = $2;
		}else{
			die "ERROR! Could not surmise run_label from job_list_file (assume looks like repair.developN)\n";
		}
	}else{
		$dir = ".";
	}
	my ($run_list_xml_doc, $root_node) = xml_create('domain_parse_run_list');

	my $domain_parse_run_node = $run_list_xml_doc->createElement('domain_parse_run');
	$domain_parse_run_node->setAttribute('run_id', 1);
	$domain_parse_run_node->appendTextChild('week_label', $run_label);

	$root_node->appendChild($domain_parse_run_node);

	my $abs_job_list_fn = File::Spec->rel2abs($job_list_file);
	$abs_job_list_fn =~ /(.*)\/([\w\.]+)/;

	my $run_list_dir = $1;
	my $run_list_job_xml_file = $2;

	$domain_parse_run_node->appendTextChild('run_list_dir', $run_list_dir);

	my $run_list_job_file_node	= $run_list_xml_doc->createElement('run_list_job_xml_file');
	$run_list_job_file_node->setAttribute('mode', 'seq_iter');
	$run_list_job_file_node->appendTextNode($run_list_job_xml_file);

	$domain_parse_run_node->appendChild($run_list_job_file_node);

	my $doc_string = $run_list_xml_doc->toString(1);

	my $out_fn = "$dir/run_list.$run_label.xml";
	xml_write($run_list_xml_doc, $out_fn);
	return $out_fn;

}

sub generate_domain_summary_files { 
	my $sub = 'generate_domain_summary_files';

	my ($domain_run_list) = @_;
	
	my $run_list_summary_doc = xml_open($domain_run_list);

	my $run_list_XPath = '//run_list_summary_doc/run_list_summary/run_list';

	my $weekly_update_dir = $run_list_summary_doc->findvalue('//weekly_update_dir');
	if (! -d $weekly_update_dir) { die "ERROR! Could not find weekly update dir $weekly_update_dir\n" } 

	my @run_data;
	my $i = 0;

	foreach my $node (find_run_list_nodes($run_list_summary_doc)) { 

		my $id 			= $node->findvalue('@id');
		my $process		= $node->findvalue('@process');
		my $mode		= $node->findvalue('mode');
		my $run_list_file	= $node->findvalue('run_list_file');
		my $run_list_label	= $node->findvalue('run_list_label');
		my $run_list_desc	= $node->findvalue('run_list_description');
		my $run_list_reference	= $node->findvalue('reference');

		if ($DEBUG) { 
			print "DEBUG $id $process $run_list_label $run_list_reference\n";
		}

		#Need to think about multiple versioning for later...
		my $domain_prefix	= $node->findvalue('domain_prefix');

		$run_data[$i]{id}		= $id;
		$run_data[$i]{process}		= $process;
		$run_data[$i]{mode}		= $mode;
		$run_data[$i]{run_list_file}	= $run_list_file;
		$run_data[$i]{run_list_label}	= $run_list_label;
		$run_data[$i]{run_list_desc}	= $run_list_desc;
		$run_data[$i]{domain_prefix}	= $domain_prefix;
		$run_data[$i]{reference}	= $run_list_reference;
		
		$i++;
	}

	for (my $i = 0; $i < scalar(@run_data); $i++) { 

		my $id			= $run_data[$i]{id};
		my $run_list_file 	= $run_data[$i]{run_list_file};
		my $run_list_label 	= $run_data[$i]{run_list_label};

		my $process 		= $run_data[$i]{process}; #Boolean

		if ($DEBUG > 3) { 
			print "DEBUG main: $id $run_list_label $run_list_file\n";
		}

		if ($process eq 'true') { 
			my $run_list_domain_summary_file = run_process_v2(\%{$run_data[$i]}, $weekly_update_dir);
			$run_data[$i]{domain_summary_file} = $run_list_domain_summary_file;
			unlink("$weekly_update_dir/run_list.latest.domain_summary.xml");
			symlink($run_list_domain_summary_file, "$weekly_update_dir/run_list.latest.domain_summary.xml");
			if (! -l "$weekly_update_dir/run_list.latest.domain_summary.xml") { die "Latest symlink not generated\n"; }
		}

	}

	my $run_list_nodes = $run_list_summary_doc->findnodes($run_list_XPath);
	foreach my $node ($run_list_nodes->get_nodelist() ) { 

		my $id 	= $node->findvalue('@id');

		if ($DEBUG) { 
			print "$id\n";
		}

		my @runs = grep {$_->{id} == $id} @run_data;
		if ($DEBUG ) { 
			foreach my $run (@runs) { 
				print "DEBUG: $id $$run{id} $$run{run_list_file}\n"; 
			}
		}
		if (scalar(@runs) > 1) { 
			die "ERROR! multiple runs with id $id?\n";
		}
		if ($runs[0]{domain_summary_file}) { 
			if ($node->exists('domain_summary_file')) { 
				foreach my $ds_node ($node->findnodes('domain_summary_file')->get_nodelist() ) { 
					$node->removeChild($ds_node);
				}
			}
			my $ds_node = $run_list_summary_doc->createElement('domain_summary_file');
			$ds_node->appendTextNode($runs[0]{domain_summary_file});
			$node->appendChild($ds_node);
		}
	}
	xml_write($run_list_summary_doc, $domain_run_list);
	system("xml_pp $domain_run_list > tmp");
	if (-s "tmp" > 0) { 
	 system("mv tmp $domain_run_list");
	}
}
sub run_process_v2 { 
		my $sub = 'run_process_v2';

		my ($run_list_href, $output_directory) = @_;

		my $id	= $$run_list_href{id};
		my $run_list_file	= $$run_list_href{run_list_file};
		my $run_list_label	= $$run_list_href{run_list_label};
		my $domain_prefix	= $$run_list_href{domain_prefix};
		my $run_list_mode	= $$run_list_href{mode};			
		my $reference		= $$run_list_href{reference};

		my $run_list_desc 	= $$run_list_href{run_list_desc};

		if ($DEBUG) { 
			print "DEBUG $sub: $id $run_list_label $run_list_file\n";
		}

		if (!-f $run_list_file) { 
			print "WARNING! $run_list_file not found, skipping...\n";
			return 0;
		}

		my $run_list_xml_doc = xml_open($run_list_file);

		#Create the document that holds all the data for a run
		my ($domain_summary_doc, $domain_summary_top_node) = xml_create('domain_parse_summary');

		#Copy the run_list_summary head data to the domain summary file
		my $ds_run_list_file_node = $domain_summary_doc->createElement('domain_summary_run_list_file');
		$ds_run_list_file_node->appendTextNode($run_list_file);

		my $ds_run_list_label_node = $domain_summary_doc->createElement('domain_summary_run_list_label');
		$ds_run_list_label_node->appendTextNode($run_list_label);

		my $ds_run_list_desc_node = $domain_summary_doc->createElement('domain_summary_run_list_description');
		$ds_run_list_desc_node->appendTextNode($run_list_desc);

		#Head node for chain list
		my $ds_chain_list_node = $domain_summary_doc->createElement('pdb_chain_list');

		$domain_summary_top_node->appendChild($ds_run_list_file_node);
		$domain_summary_top_node->appendChild($ds_run_list_label_node);
		$domain_summary_top_node->appendChild($ds_run_list_desc_node);
		$domain_summary_top_node->appendChild($ds_chain_list_node);

		#Go read the run_list summary, the associated hora job list summaries, and the associate domain xmls and build
		my $run_list_XPath = '//domain_parse_run_list/domain_parse_run';

		#ref library read
		my $ref_domain_list = reference_domain_cache_load($reference);

		if ($run_list_xml_doc->exists($run_list_XPath)) { 
			my $run_list_nodes = $run_list_xml_doc->findnodes($run_list_XPath);

			foreach my $run_node ($run_list_nodes->get_nodelist() ) { 
				my $run_list_dir	= $run_node->findvalue('run_list_dir');
				my $run_id 		= $run_node->findvalue('@run_id');
				my $week_label 		= $run_node->findvalue('week_label');

				if ($DEBUG) { 
					print "DEBUG $sub: $run_id $run_list_dir\n";
				}
				if (!-d $run_list_dir) { 
					print "WARNING! $sub: id $id dir $run_list_dir not found, skipping...\n";
					next;
				}
				foreach my $jx_file_node ($run_node->findnodes('run_list_job_xml_file')->get_nodelist() ) { 
					my $run_mode = $jx_file_node->findvalue('@mode');
					if ($DEBUG ) { 
						print "DEBUG $sub: $run_id $run_mode\n";
					}
					#domain_mode_process_v2($domain_summary_doc, $ref_xml_doc, $run_mode, $reference, $run_node, $run_list_dir, $run_list_label, $domain_prefix);
					domain_mode_process_v2($domain_summary_doc, $ref_domain_list, $run_mode, $reference, $run_node, $run_list_dir, $run_list_label, $domain_prefix);
				}
				#Merge domain parse jobs into a single category, seq_iter only, struct_search only, and mixed. Relabel domain IDs 
				#domain_merge($domain_summary_doc, $week_label);
				#Only use all of one type
				domain_merge_v2($domain_summary_doc, $reference, $week_label);
			}
			add_week_list($domain_summary_doc);
			
			#Postprocess ranges, renumber domains, label
			#domain_simple_postprocess($domain_summary_doc);

		}else{
			print "WARNING! $run_list_file doesn't have XPath $run_list_XPath, skipping...\n";
			return 0;
		}
		my $domain_summary_file_name = "$output_directory/$run_list_label.$domain_prefix.domain_summary.xml";

		xml_write($domain_summary_doc, $domain_summary_file_name);
		return $domain_summary_file_name;
}
sub hasDomains {
	$_[0]->exists('.//domain');
}
sub domain_merge_v2 { 
		my $sub = 'domain_merge_v2';

		my ($domain_summary_doc, $reference, $week_label) = @_;

		my $pdb_chain_XPath = qq{//pdb_chain[\@week_label='$week_label']};

		if ($DEBUG) { 
			print "DEBUG $sub: $reference $week_label\n";
		}

		foreach my $pc_node ($domain_summary_doc->findnodes($pdb_chain_XPath)->get_nodelist() ) { 

			my $pdb		= $pc_node->findvalue('@pdb');
			my $chain	= $pc_node->findvalue('@chain');

			my @types;
			my %dp_domains;
			foreach my $dp_node ($pc_node->findnodes('domain_parse')) { 

				my $type = $dp_node->findvalue('@type');

				next if !hasDomains($dp_node);

				push (@types, $type);
				
				foreach my $dom_node ($dp_node->findnodes('job/domain_list/domain')) { 
					push (@{$dp_domains{$type}}, $dom_node);
				}

			}

			my $merge_node	= $domain_summary_doc->createElement('domain_parse');
			$merge_node->setAttribute('type', 'merge');

			my $job_node	= $domain_summary_doc->createElement('job');
			$job_node->setAttribute('mode', 'merge');
			$job_node->setAttribute('reference', $reference);

			$merge_node->appendChild($job_node);

			my $dom_list_node	= $domain_summary_doc->createElement('domain_list');
			$job_node->appendChild($dom_list_node);
			
			my $good_types;
			if (scalar(@types) == 1) { 
				foreach my $dom_node (@{$dp_domains{$types[0]}}) { 
					my $new_dom_node = $dom_node->cloneNode(1);
					$dom_list_node->appendChild($new_dom_node);
				}
				#my $dp_node 	= $pc_node->findnodes(qq{domain_parse[\@type="$types[0]"]})->get_node(1);

				my $dp_cdc_XPath = qq{domain_parse[\@type="$types[0]"]/chain_domain_coverage};
				my $dp_cdc_node;

				my $dp_ocdc_XPath = qq{domain_parse[\@type="$types[0]"]/optimized_chain_domain_coverage};
				my $dp_ocdc_node;

				if ($pc_node->exists($dp_cdc_XPath)) { 
					$dp_cdc_node = $pc_node->findnodes($dp_cdc_XPath)->get_node(1);
				}else{
					$dp_cdc_node = $domain_summary_doc->createElement('chain_domain_coverage');
					$dp_cdc_node->appendTextNode('0');
				}

				if ($pc_node->exists($dp_ocdc_XPath)) { 
					$dp_ocdc_node = $pc_node->findnodes($dp_ocdc_XPath)->get_node(1);
				}else{
					$dp_ocdc_node = $domain_summary_doc->createElement('optimized_chain_domain_coverage');
					$dp_ocdc_node->setAttribute('used_res', '0');
				}
				my $new_cdc_node = $dp_cdc_node->cloneNode(1);
				my $new_ocdc_node = $dp_ocdc_node->cloneNode(1);

				$merge_node->appendChild($new_cdc_node);
				$merge_node->appendChild($new_ocdc_node);
			}else{
				#Take either all the seq_iter or all the struct_search nodes	
				my @good_types;
				foreach my $type (@types) { 
					#my $dp_node 	= $pc_node->findnodes(qq{domain_parse[\@type="$type"]})->get_node(1);
					#my $cdc 	= $dp_node->findvalue('chain_domain_coverage');	

					my $dp_node 	= $pc_node->findnodes(qq{domain_parse[\@type="$type"]})->get_node(1);
					if ($dp_node->exists('optimized_chain_domain_coverage')) { 
						my $ocdc = $dp_node->findvalue('optimized_chain_domain_coverage/@unused_res');
						if ($ocdc && $ocdc <= 25 ) {  # ==0 to <= 25 Thursday, May 09 2013
							push (@good_types, $type);
						}
					}
				}
				
				my $default1 = 'struct_search';
				my $default2 = 'seq_iter';
				my $use_type;
				if (scalar(@good_types) == 2)  { 
					$use_type = $default1;
				}elsif(scalar(@good_types) == 1) { 
					$use_type = $good_types[0];
				}else{
					$use_type = $default2;
				}

				#my $dp_node 	= $pc_node->findnodes(qq{domain_parse[\@type="$use_type"]})->get_node(1);

				my $dp_cdc_XPath = qq{domain_parse[\@type="$use_type"]/chain_domain_coverage};
				my $dp_cdc_node;

				my $dp_ocdc_XPath = qq{domain_parse[\@type="$use_type"]/optimized_chain_domain_coverage};
				my $dp_ocdc_node;

				if ($pc_node->exists($dp_cdc_XPath)) { 
					$dp_cdc_node = $pc_node->findnodes($dp_cdc_XPath)->get_node(1);
					
				}else{
					$dp_cdc_node = $domain_summary_doc->createElement('chain_domain_coverage');
					$dp_cdc_node->appendTextNode('0');
				}
				my $new_cdc_node = $dp_cdc_node->cloneNode(1);

				if ($pc_node->exists($dp_ocdc_XPath)) { 
					$dp_ocdc_node = $pc_node->findnodes($dp_ocdc_XPath)->get_node(1);
				}else{
					$dp_ocdc_node = $domain_summary_doc->createElement('optimized_chain_domain_coverage');
					$dp_ocdc_node->setAttribute('used_res', 0);
					$dp_ocdc_node->appendTextNode('0');
				}
				my $new_ocdc_node = $dp_ocdc_node->cloneNode(1);

				foreach my $dom_node ( @{$dp_domains{$use_type}}) { 

					my $new_dom_node = $dom_node->cloneNode(1);
					$dom_list_node->appendChild($new_dom_node);
				}
				$merge_node->appendChild($new_cdc_node);
				$merge_node->appendChild($new_ocdc_node);
			}
			$pc_node->appendChild($merge_node);
		}
			
		$pdb_chain_XPath = qq{//pdb_chains[\@week_label='$week_label']};

		foreach my $pc_node ($domain_summary_doc->findnodes($pdb_chain_XPath)->get_nodelist() ) { 

		my $pdb		= $pc_node->findvalue('@pdb');
		my $chain	= $pc_node->findvalue('@chain');

		my @types;
		my %dp_domains;
		foreach my $dp_node ($pc_node->findnodes('domain_parse')) { 

			my $type = $dp_node->findvalue('@type');

			if (!$dp_node->exists('job_asm/domain_list/domain')) { next } 

			push (@types, $type);
			
			foreach my $dom_node ($dp_node->findnodes('job_asm/domain_list/domain')) { 
				push (@{$dp_domains{$type}}, $dom_node);
			}

		}

		my $merge_node	= $domain_summary_doc->createElement('domain_parse');
		$merge_node->setAttribute('type', 'merge');

		my $job_node	= $domain_summary_doc->createElement('job');
		$job_node->setAttribute('mode', 'merge');
		$job_node->setAttribute('reference', $reference);

		$merge_node->appendChild($job_node);

		my $dom_list_node	= $domain_summary_doc->createElement('domain_list');
		$job_node->appendChild($dom_list_node);
		

		my $good_types;
		if (scalar(@types) == 1) { 
			foreach my $dom_node (@{$dp_domains{$types[0]}}) { 
				my $new_dom_node = $dom_node->cloneNode(1);
				$dom_list_node->appendChild($new_dom_node);
			}
			#my $dp_node 	= $pc_node->findnodes(qq{domain_parse[\@type="$types[0]"]})->get_node(1);

			my $dp_cdc_XPath = qq{domain_parse[\@type="$types[0]"]/chain_domain_coverage};
			my $dp_cdc_node;

			my $dp_ocdc_XPath = qq{domain_parse[\@type="$types[0]"]/optimized_chain_domain_coverage};
			my $dp_ocdc_node;

			if ($pc_node->exists($dp_cdc_XPath)) { 
				$dp_cdc_node = $pc_node->findnodes($dp_cdc_XPath)->get_node(1);
			}else{
				$dp_cdc_node = $domain_summary_doc->createElement('chain_domain_coverage');
				$dp_cdc_node->appendTextNode('0');
			}

			if ($pc_node->exists($dp_ocdc_XPath)) { 
				$dp_ocdc_node = $pc_node->findnodes($dp_ocdc_XPath)->get_node(1);
			}else{
				$dp_ocdc_node = $domain_summary_doc->createElement('optimized_chain_domain_coverage');
				$dp_ocdc_node->setAttribute('used_res', '0');
			}
			my $new_cdc_node = $dp_cdc_node->cloneNode(1);
			my $new_ocdc_node = $dp_ocdc_node->cloneNode(1);

			$merge_node->appendChild($new_cdc_node);
			$merge_node->appendChild($new_ocdc_node);
		}else{
			#Take either all the seq_iter or all the struct_search nodes	
			my @good_types;
			foreach my $type (@types) { 
				#my $dp_node 	= $pc_node->findnodes(qq{domain_parse[\@type="$type"]})->get_node(1);
				#my $cdc 	= $dp_node->findvalue('chain_domain_coverage');	

				my $dp_node 	= $pc_node->findnodes(qq{domain_parse[\@type="$type"]})->get_node(1);
				if ($dp_node->exists('optimized_chain_domain_coverage')) { 
					my $ocdc = $dp_node->findvalue('optimized_chain_domain_coverage/@unused_res');
					if ($ocdc && $ocdc <= 25 ) {  # ==0 to <= 25 Thursday, May 09 2013
						push (@good_types, $type);
					}
				}
			}
			
			my $default1 = 'struct_search';
			my $default2 = 'seq_iter';
			my $use_type;
			if (scalar(@good_types) == 2)  { 
				$use_type = $default1;
			}elsif(scalar(@good_types) == 1) { 
				$use_type = $good_types[0];
			}else{
				$use_type = $default2;
			}

			#my $dp_node 	= $pc_node->findnodes(qq{domain_parse[\@type="$use_type"]})->get_node(1);

			my $dp_cdc_XPath = qq{domain_parse[\@type="$use_type"]/chain_domain_coverage};
			my $dp_cdc_node;

			my $dp_ocdc_XPath = qq{domain_parse[\@type="$use_type"]/optimized_chain_domain_coverage};
			my $dp_ocdc_node;

			if ($pc_node->exists($dp_cdc_XPath)) { 
				$dp_cdc_node = $pc_node->findnodes($dp_cdc_XPath)->get_node(1);
				
			}else{
				$dp_cdc_node = $domain_summary_doc->createElement('chain_domain_coverage');
				$dp_cdc_node->appendTextNode('0');
			}
			my $new_cdc_node = $dp_cdc_node->cloneNode(1);

			if ($pc_node->exists($dp_ocdc_XPath)) { 
				$dp_ocdc_node = $pc_node->findnodes($dp_ocdc_XPath)->get_node(1);
			}else{
				$dp_ocdc_node = $domain_summary_doc->createElement('optimized_chain_domain_coverage');
				$dp_ocdc_node->setAttribute('used_res', 0);
				$dp_ocdc_node->appendTextNode('0');
			}
			my $new_ocdc_node = $dp_ocdc_node->cloneNode(1);

			foreach my $dom_node ( @{$dp_domains{$use_type}}) { 

				my $new_dom_node = $dom_node->cloneNode(1);
				$dom_list_node->appendChild($new_dom_node);
			}
			$merge_node->appendChild($new_cdc_node);
			$merge_node->appendChild($new_ocdc_node);
		}
		$pc_node->appendChild($merge_node);
	}
		
}
			

sub reference_library_load { 
	my $sub = 'reference_library_load';

	my ($reference) = @_;
	
	my $ref_loc;
	if ($REF_XML{$reference} && -f $REF_XML{$reference}) { 
		$ref_loc = $REF_XML{$reference};
	}else{
		die "ERROR! $sub: no ref $reference\n";
	}

	my $ref_xml_doc = xml_open($ref_loc);

	return $ref_xml_doc;
}

sub reference_domain_cache_load{ 
	my $sub = 'reference_domain_cache_load';

	my ($reference) = @_;
	
	my $cache_loc;
	if ($REF_RANGE_CACHE{$reference} && -f $REF_RANGE_CACHE{$reference}) { 
		$cache_loc = $REF_RANGE_CACHE{$reference};
	}

	open (my $fh, "<", $cache_loc) or die "ERROR! Could not open $cache_loc for reading:$!\n;";

	my %range_cache;
	while (my $ln = <$fh>) { 
		$ln =~ /^#/ and next;
		my @F = split(/\s+/, $ln);
		if (scalar(@F) != 5) { next } 
		my $uid	= $F[0];
		my $ecod_domain_id = $F[1];
		my $range = $F[2];
		my $pdb = $F[3];
		my $chain = $F[4];

		$range_cache{$ecod_domain_id} = $range;
	}
	
	return \%range_cache;
}
sub add_week_list { 

	my $sub = 'add_week_list';
	
	my ($domain_summary_doc) = @_;

	my %week_label;

	my $dom_summary_doc_top_node = $domain_summary_doc->findnodes('//domain_parse_summary')->get_node(1);

	if($domain_summary_doc->exists('//pdb_chain[@week_label]')) { 
		my $week_label_nodes	= $domain_summary_doc->findnodes('//pdb_chain[@week_label]');

		foreach my $node ($week_label_nodes->get_nodelist() ) { 
			my $week_label = $node->findvalue('@week_label');
			$week_label{$week_label}++;
		}
		my $ds_week_list_node = $domain_summary_doc->createElement('week_list');

		my @keys = sort {$a cmp $b} keys %week_label;
		for (my $i = 0; $i < (scalar(@keys)); $i++) { 
			my $key = $keys[$i];
			my $id = $i+1;
			my $week_label_node = $domain_summary_doc->createElement('week');

			$week_label_node->setAttribute('id', $id);
			$week_label_node->setAttribute('label', $key);

			$ds_week_list_node->appendChild($week_label_node);
		}
		$dom_summary_doc_top_node->appendChild($ds_week_list_node);
	}else{
		print "WARNING! $sub: No week label attributes found\n";
		return 0;
	}
	return 1;
}
			
sub rectify_range_overlap { 
	my $sub = 'rectify_range_overlap';

	my ($range_aref1, $range_aref2) = @_;

	if ($DEBUG) { 
		my $range1 = rangify(@$range_aref1);
		my $range2 = rangify(@$range_aref2);
		print "DEBUG $sub: r1 $range1 r2 $range2\n";
	}


	my @overlap;
	for (my $i = 0; $i < scalar(@$range_aref1); $i++) { 
		for (my $j = 0; $j < scalar(@$range_aref2); $j++) { 
			if ($$range_aref1[$i] == $$range_aref2[$j]) { 
				push (@overlap, $$range_aref1[$i]);
			}
		}
	}
	my $overlap_range_str = rangify(@overlap);
	my @overlap_segs;
	if ($overlap_range_str =~ /\,/) { 
		@overlap_segs = split(/\,/, $overlap_range_str);
	}else{
		push (@overlap_segs, $overlap_range_str);
	}

	if ($DEBUG) { 
		printf "DEBUG $sub: %i residue overlap\n", scalar(@overlap);
	}

	#my $split = scalar(@overlap) % 2;
	foreach my $seg (@overlap_segs) { 
		if ($DEBUG) { 
			print "DEBUG $sub: $seg\n";
		}
		if ($seg =~ /\-/) { 
			my $seg_overlap_aref = range_expand($seg);
				
			my $split = int(scalar(@$seg_overlap_aref)/2);
			my (@one, @two);
			for (my $i = 0; $i < scalar(@$seg_overlap_aref); $i++) { 
				if ($i < $split) { 
					push (@one, $$seg_overlap_aref[$i]);
					#range_exclude(\@{$overlap[$i]}, $range_aref2);
				}else{
					push (@two, $$seg_overlap_aref[$i]);
					#range_exclude(\@{$overlap[$i]}, $range_aref1);
				}
			}
			if ($$range_aref1[0] < $$range_aref2[0]) { 
				range_exclude(\@one, $range_aref2);
				range_exclude(\@two, $range_aref1);
			}else{
				range_exclude(\@one, $range_aref1);
				range_exclude(\@two, $range_aref2);
			}
		}else{
			my @singleton;
			push (@singleton, $seg);
			range_exclude(\@singleton, $range_aref2);
		}
			
	}

	if ($DEBUG) { 
		my $rect1 = rangify(@$range_aref1);
		my $rect2 = rangify(@$range_aref2);
		print "DEBUG $sub: rect1 $rect1: rect2 $rect2\n";
	}


	return($range_aref1, $range_aref2);

}



sub best_hit_comp { 
	my $sub = 'best_hit_comp';

	my ($scop_domid1, $scop_domid2, $ver) = @_;

	$scop_domid1 = lc($scop_domid1);
	$scop_domid2 = lc($scop_domid2);

	my $scop_cla_file;
	if ($ver eq 'scop1.69') { 
		$scop_cla_file = '/home/rschaeff/data/scop/v1.69/dir.cla.scop.txt_1.69.txt';
	}elsif($ver eq 'scop1.75') { 
		$scop_cla_file = '/home/rschaeff/data/scop/v1.75/dir.cla.scop.txt_1.75.txt';
	}else{
		die "ERROR! $sub: Bad reference $ver (should be scop1.69/scop1.75)\n";
	}
	if (!-f $scop_cla_file) { 
		die "ERROR! $sub: SCOP cla file $scop_cla_file not found\n";
	}

	open (my $fh, "<", $scop_cla_file) or die "ERROR! $sub: Could not open $scop_cla_file for reading:$!\n";

	my ($dom1_ln, $dom2_ln);
	while (my $ln = <$fh>) { 
		p
		$ln =~ /^#/ and next;
	
		my @F = split(/\s+/, $ln);

		if($F[0] eq $scop_domid1) {  
			$dom1_ln	= $ln;
		}
		if($F[0] eq $scop_domid2) { 
			$dom2_ln	= $ln;
		}
	}
	close $fh;

	if (!$dom1_ln) { 
		print "WARNING! $sub: Domain information for $scop_domid1 not found in $scop_cla_file\n";
	}
	if (!$dom2_ln) { 
		print "WARNING! $sub: Domain information for $scop_domid2 not found in $scop_cla_file\n";
	}

	my @dom1 = split(/\s+/, $dom1_ln);
	my @dom2 = split(/\s+/, $dom2_ln);

	my @cla1;
	$dom1[5] =~ /cl=(\d+),cf=(\d+),sf=(\d+),fa=(\d+),dm=(\d+),sp=(\d+),px=(\d+)/;
	$cla1[0] = $1;
	$cla1[1] = $2;
	$cla1[2] = $3;
	$cla1[3] = $4;
	$cla1[4] = $5;
	$cla1[5] = $6;
	$cla1[6] = $7;
	$cla1[7] = $8;

	my @cla2;
	$dom2[5] =~ /cl=(\d+),cf=(\d+),sf=(\d+),fa=(\d+),dm=(\d+),sp=(\d+),px=(\d+)/;
	$cla2[0] = $1;
	$cla2[1] = $2;
	$cla2[2] = $3;
	$cla2[3] = $4;
	$cla2[4] = $5;
	$cla2[5] = $6;
	$cla2[6] = $7;
	$cla2[7] = $8;

	if ($cla1[2] eq $cla2[2]) { 
		return 1;
	}

	return 0;

}


			
sub domain_mode_process_v2 { 
	my $sub = 'domain_mode_process_v2';

	my ($domain_summary_doc, $reference_doc, $mode, $reference, $run_node, $run_list_dir, $run_list_label, $domain_prefix) = @_;

	my $rl_job_file;
	my $rl_job_file_XPath;

	$rl_job_file_XPath = "run_list_job_xml_file[\@mode=\'$mode\']";
	if (existsMode($run_node, $mode)) { 
		$rl_job_file = $run_node->findvalue($rl_job_file_XPath);
		if (!-f "$run_list_dir/$rl_job_file") { 
			print "WARNING! $sub: $mode job xml file $rl_job_file not found, skipping $run_list_dir...\n";
			return 0;
		}
	}else{
		print "WARNING! $sub: no job xml file found for $run_list_dir, $domain_prefix, $mode, skipping...\n";
		return 0;
	}

	my $rl_week_label;
	my $rl_week_label_XPath = 'week_label';
	if ($run_node->exists($rl_week_label_XPath)) { 
		$rl_week_label = $run_node->findvalue($rl_week_label_XPath);
	}else{
		print "WARNING! $sub: No week label found for $run_list_dir, $domain_prefix, $mode...\n";
	}


	open (my $xml_fh, "<",  "$run_list_dir/$rl_job_file") or die "ERROR! $sub: Could not open $rl_job_file for reading:$!\n";
	my $job_xml_doc = XML::LibXML->load_xml( IO => $xml_fh );
	close $xml_fh;
	print "DEBUG_FOUL: $run_list_dir/$rl_job_file\n";

	my $job_XPath;
	if ($mode eq 'hora') { 
		$job_XPath = '//hora_job_set_top';
	}else{
		$job_XPath = '//job_set_top';
	}

	my $bad_domain_xml_file_count = 0;
	if ($job_xml_doc->exists($job_XPath)) { 
		my $fhx_doc_nodes = $job_xml_doc->findnodes($job_XPath);

		my $fhx_doc_node;
		if ($fhx_doc_nodes->size() > 1){ 
			die "ERROR! $sub: No support for multiple job docs per file in $rl_job_file, aborting...\n";
		}else{
			$fhx_doc_node = $fhx_doc_nodes->get_node(1);
		}


		my $job_list_dir;
		my $job_list_dir_XPath;
		if ($mode eq 'hora') { 
			$job_list_dir_XPath = 'hora_job_list_dir';
		}else{
			$job_list_dir_XPath = 'job_list_dir';
		}
		if ($fhx_doc_node->exists($job_list_dir_XPath)) { 
			$job_list_dir = $fhx_doc_node->findvalue($job_list_dir_XPath);
			if (!-d $job_list_dir) { 
				print "ERROR! $sub: job list dir $job_list_dir not found\n";
				next;
			}
		}else{
			print "ERROR! $sub: Fouled XPath/Broken HorA xml for job list dir on $rl_job_file\n";
			next;
		}

		my $job_dump_dir;
		my $job_dump_dir_XPath;
		if ($mode eq 'hora') { 
			$job_dump_dir_XPath = 'hora_job_dump_dir';
		}else{
			$job_dump_dir_XPath = 'job_dump_dir';
		}
		if ($fhx_doc_node->exists($job_dump_dir_XPath)) { 
			$job_dump_dir = $fhx_doc_node->findvalue($job_dump_dir_XPath);
			if (!-d $job_dump_dir) { 
				print "ERROR! $sub: job dump dir $job_dump_dir not found\n";
				next;
			}
		}else{
			print "ERROR! $sub: Fouled XPath/Broken HorA xml for job dump dir on $rl_job_file\n";
			next;
		}

		my $job_XPath = 'job_list/job';
		my $job_asm_XPath = 'job_asm_list/job_asm';

		my $job_nodes = $fhx_doc_node->findnodes($job_XPath);
		my $job_asm_nodes = $fhx_doc_node->findnodes($job_asm_XPath);

		if ($DEBUG) { 
			printf "DEBUG $sub: Found %i jobs in run list dir $run_list_dir for $run_list_label(mode $mode/ reference $reference) \n", $job_nodes->size();
		}

		#HACKY
		my $reps = 0;
		if ($fhx_doc_node->findvalue('//@rep95')) { 
			$reps = 1;
			print "DEBUG: rep95 exists for $job_list_dir\n";
		}

		foreach my $job_node ($job_nodes->get_nodelist() ) { 

			#if ($reps && $job_node->findvalue('@rep95') ne 'true') { next } 

			my $job_inx		= $job_node->findvalue('@id');

			my $query_pdb		= $job_node->findvalue('query_pdb');
			my $query_chain		= $job_node->findvalue('query_chain');

			my $pdb_chain 		= "${query_pdb}_${query_chain}";

			my $job_mode		= $job_node->findvalue('mode');
			my $reference		= $job_node->findvalue('reference');

			my $obsolete 		= $job_node->findvalue('@pdb_obsolete');
			if ($obsolete eq 'true') { next } 
			
			if ($DEBUG > 2) { 
				print "DEBUG $sub: $query_pdb $query_chain $job_mode $reference\n";
			}
			my $rep95 = 'false';
			if ($job_node->exists('@rep95')) { 
				$rep95 = $job_node->findvalue('@rep95');
			}

			my $peptide_filter = 'false';
			if ($job_node->exists('peptide_filter')) { 
				my $peptide_filter_nodes = $job_node->findnodes('peptide_filter');
				my @vers;
				foreach my $pf_node ($peptide_filter_nodes->get_nodelist() ) { 
					push (@vers, $pf_node->findvalue('@ver'));
				}
				@vers = sort {$b <=> $a} @vers;
				$peptide_filter_nodes = $job_node->findnodes("peptide_filter[\@ver='$vers[0]']");
				my $peptide_filter_node;
				if ($peptide_filter_nodes->size() > 1) { 
					warn "WARNING! multiple same version peptide filter nodes in $job_inx, $run_list_dir, weird...\n";
					$peptide_filter_node = $peptide_filter_nodes->get_node(1);
				}else{
					$peptide_filter_node = $peptide_filter_nodes->get_node(1);
				}
				$peptide_filter = $peptide_filter_node->findvalue('@apply');
			}
			my $coiled_coil_filter = 'false';
			if ($job_node->exists('coiled_coil_filter')) { 
				my $cc_nodes = $job_node->findnodes('coiled_coil_filter');
				my @vers;
				foreach my $cc_node ($cc_nodes->get_nodelist()) { 
					push (@vers, $cc_node->findvalue('@ver'));
				}
				@vers = sort {$b <=> $a} @vers;
				$cc_nodes = $job_node->findnodes("coiled_coil_filter[\@ver='$vers[0]']");
				my $cc_node;
				if ($cc_nodes->size() > 1) { 
					die "ERROR! multiple same version cc filter nodes in $job_inx, $run_list_dir, weird...\n";
				}else{
					$cc_node = $cc_nodes->get_node(1);
				}
				$coiled_coil_filter = $cc_node->findvalue('@apply');
			}

			my $mc_domain_component = 'false';
			if ($job_node->findvalue('@sw_mc_complex') eq 'true') { 
				$mc_domain_component = 'true';
			}


			#Generate some summary nodes

			my $ds_pdb_chain_list_XPath = '//domain_parse_summary/pdb_chain_list';	
			my ($ds_pdb_chain_list_node, $ds_pdb_chain_node, $ds_job_node, $ds_domain_parse_node);
			if ($domain_summary_doc->exists($ds_pdb_chain_list_XPath)) { 
				$ds_pdb_chain_list_node = $domain_summary_doc->findnodes($ds_pdb_chain_list_XPath)->get_node(1);

				my $ds_pdb_chain_XPath = "pdb_chain[\@pdb=\'$query_pdb\'][\@chain=\'$query_chain\']";

				if ($ds_pdb_chain_list_node->exists($ds_pdb_chain_XPath)) { 
					my $ds_pdb_chain_nodes = $ds_pdb_chain_list_node->findnodes($ds_pdb_chain_XPath);
					if ($ds_pdb_chain_nodes->size() > 1) { 
						die "ERROR! $sub: More than one pdb_chain node for $query_pdb, $query_chain... $run_list_dir\n";
					}else{
						$ds_pdb_chain_node = $ds_pdb_chain_nodes->get_node(1);
					}
				}else{
					$ds_pdb_chain_node = $domain_summary_doc->createElement('pdb_chain');
		
					$ds_pdb_chain_node->setAttribute('pdb', $query_pdb);
					$ds_pdb_chain_node->setAttribute('chain', $query_chain);
					$ds_pdb_chain_node->setAttribute('week_label', $rl_week_label);
					$ds_pdb_chain_node->setAttribute('peptide_filter', $peptide_filter);
					$ds_pdb_chain_node->setAttribute('coiled_coil_filter', $coiled_coil_filter);
					$ds_pdb_chain_node->setAttribute('mc_domain_component', $mc_domain_component);
					$ds_pdb_chain_node->setAttribute('rep95', $rep95);

					$ds_pdb_chain_list_node->appendChild($ds_pdb_chain_node);
				}

				my $ds_domain_parse_XPath = qq{domain_parse[\@type='parser']/job[\@mode='$job_mode'][\@reference='$reference\']};
				if ($ds_pdb_chain_node->exists($ds_domain_parse_XPath)) { 
					die "WARNING! $sub: domain parse for $query_pdb,$query_chain,$job_mode,$reference,parser already exists?\n";
				}else{
					$ds_job_node = $domain_summary_doc->createElement('job');
					$ds_domain_parse_node = $domain_summary_doc->createElement('domain_parse');

					$ds_job_node->setAttribute('mode', $job_mode);
					$ds_job_node->setAttribute('reference', $reference);
					#my $type 	= "parser_$mode";
					my $type	= $mode;
					$ds_domain_parse_node->setAttribute('type', $type);

					$ds_pdb_chain_node->appendChild($ds_domain_parse_node);
					$ds_domain_parse_node->appendChild($ds_job_node);
				}

			}
			

			
			my $domain_xml_file;
			#This is not the best way to do this.
			if ($mode eq 'seq_iter') { 
				$domain_xml_file	= "$job_dump_dir/$pdb_chain/$domain_prefix.$pdb_chain.$reference.xml";
			}elsif($mode eq 'struct_search') { 
				$domain_xml_file	= "$job_dump_dir/$pdb_chain/$domain_prefix.$pdb_chain.$reference.struct_search.xml";
			}




			if (!-f $domain_xml_file) { 
				#print "WARNING! $sub: Domain XML file $domain_xml_file not found!\n";
				$bad_domain_xml_file_count++; #Just missing, usually because not rep95 chain
				next;
			}

			if (-s $domain_xml_file == 0) { 
				print "WARNING! DOMAIN XML file $domain_xml_file empty!\n";
			}


			my $domain_xml_doc = xml_open($domain_xml_file);
			
			my $domain_doc_XPath	= "//chain_domains_set_top[\@domain_prefix='$domain_prefix']";
			if ($domain_xml_doc->exists($domain_doc_XPath)) { 
				my $domain_doc_nodes = $domain_xml_doc->findnodes($domain_doc_XPath);

				my $domain_doc_node;
				if ($domain_doc_nodes->size() > 1) { 
					die "ERROR! $sub: No support for multiple doc nodes per domain file in $domain_xml_file, aborting...\n";
				}else{
					$domain_doc_node = $domain_doc_nodes->get_node(1);
				}

				if ($domain_xml_doc->findvalue('@known_reference_peptide') eq 'true') { 
					$ds_pdb_chain_node->setAttribute('peptide_filter', 'true');
				}

				my $pdb_id	= lc($domain_doc_node->findvalue('@pdb_id'));
				my $chain_id	= $domain_doc_node->findvalue('@chain_id');
				my $domain_mode = $domain_doc_node->findvalue('@mode');

				#Sanity
				if (lc($pdb_id) ne lc($query_pdb)) { 
					print "WARNING $sub! PDB I/O mismatch for $pdb_id/$query_pdb \n";
				}

				if ($chain_id ne $query_chain){ 
					print "WARNING $sub! CHAIN I/O mismatch for $chain_id/$query_chain";
				}

				if ($job_mode ne $domain_mode && $domain_mode ne 'blast') { 
					print "WARNING $sub! MODE I/O mismatch for $job_mode/$domain_mode for $job_inx\n";
				}

				#my $coverage	= $domain_doc_node->findvalue('chain_domain_coverage');
				if (!$domain_doc_node->exists('chain_domain_coverage')) { 
					die "ERROR! $sub: No coverage node for $pdb_id $chain_id\n";
				}
				my $coverage_node = $domain_doc_node->findnodes('chain_domain_coverage')->get_node(1);
				if (!$ds_domain_parse_node->exists('chain_domain_coverage')) { 
					$ds_domain_parse_node->appendChild($coverage_node);
				}else{
					die "ERROR! $sub: Chain coverage node already exists for $job_inx, $query_pdb, $query_chain, $run_list_dir?\n";
				}

				if (!$domain_doc_node->exists('optimized_chain_domain_coverage')) { 
					#print "WARNING! $sub: No optimized coverage node for $pdb_id $chain_id\n";
				}else{
					my $optimized_coverage_node = $domain_doc_node->findnodes('optimized_chain_domain_coverage')->get_node(1);
					if (!$ds_domain_parse_node->exists('optimized_chain_domain_coverage')) { 
						$ds_domain_parse_node->appendChild($optimized_coverage_node);
					}else{
						die "ERROR! $sub: Optimized chain coverage node already exists for $job_inx, $query_pdb, $query_chain, $run_list_dir?\n";
					}
				}


				my $domain_list_XPath = 'domain_list/domain';
				my $domain_list_nodes = $domain_doc_node->findnodes($domain_list_XPath);

				my $ds_domain_list_node;
				if (!$ds_job_node->exists('domain_list')) { 
					$ds_domain_list_node = $domain_summary_doc->createElement('domain_list');
					$ds_job_node->appendChild($ds_domain_list_node);
				}else{
					die "ERROR! $sub: Domain list already exists for $job_inx, $query_pdb, $query_chain, $run_list_dir?\n";
				}

				foreach my $domain_node ($domain_list_nodes->get_nodelist() ) { 

					my $domain_id 		= $domain_node->findvalue('@domain_id');
					my $hit_domain_id	= $domain_node->findvalue('hit_domain/@ecod_domain_id');	
					my $method		= $domain_node->findvalue('method');

					
					if ($DEBUG > 2) { 
						print "DEBUG: $domain_id $hit_domain_id $method\n";
					}	

					if ($method !~ /self_comp/ && !$domain_node->exists('hit_domain/domain')) { 
						if ($$reference_doc{$hit_domain_id}) { 
							#my $hit_domain_node	= $domain_node->findnodes('hit_domain')->get_node(1);
							#my $ref_hit_domain_node	= $reference_doc->findnodes(qq{//domain[\@ecod_domain_id='$hit_domain_id']})->get_node(1);
							#$hit_domain_node->appendChild($ref_hit_domain_node);
							#my $new_ref_node = $ref_hit_domain_node->cloneNode(1);
							#my $new_ref_node = $domain_summary_doc->createElement('hit_domain');

							#$hit_domain_node->appendChild($new_ref_node);
						}else{
							print "WARNING! $hit_domain_id not found in reference?  $reference $domain_id ($rl_week_label)\n";
						}
					}
					
					$ds_domain_list_node->appendChild($domain_node);
				}
				#Coils

				my $coil_list_XPath = 'coil_list/coil';
				my $coil_list_nodes = $domain_doc_node->findnodes($coil_list_XPath);

				my $ds_coil_list_node;
				if (!$ds_job_node->exists('coil_list')) { 
					$ds_coil_list_node = $domain_summary_doc->createElement('coil_list');
					$ds_job_node->appendChild($ds_coil_list_node);
				}else{
					die "ERROR! $sub: Coil list already exists for $job_inx, $query_pdb, $query_chain, $run_list_dir?\n";
				}

				foreach my $coil_node ($coil_list_nodes->get_nodelist() )  {
					my $coil_id	= $coil_node->findvalue('@coil_id');
					my $method	= $coil_node->findvalue('method');

					$ds_coil_list_node->appendChild($coil_node);
				}
			}else{
				print "ERROR! $sub: Fouled XPath/broken domain XML for top node on $domain_xml_file\n";
			}
				
		}	
		
		foreach my $job_asm_node ($job_asm_nodes->get_nodelist() ) { 

			my $job_inx		= $job_asm_node->findvalue('@id');

			my $query_pdb		= $job_asm_node->findvalue('query_pdb');
			my $query_chains	= $job_asm_node->findvalue('query_chains');
			$query_chains =~ s/,//g;

			my $pdb_chains 		= "${query_pdb}_${query_chains}";

			my $job_mode		= $job_asm_node->findvalue('mode');
			my $reference		= $job_asm_node->findvalue('reference');

			my $obsolete 		= $job_asm_node->findvalue('@pdb_obsolete');
			if ($obsolete eq 'true') { next } 
			
			if ($DEBUG > 2) { 
				print "DEBUG $sub: $query_pdb $query_chains $job_mode $reference\n";
			}
			my $rep95 = 'false';
			if ($job_asm_node->exists('@rep95')) { 
				$rep95 = $job_asm_node->findvalue('@rep95');
			}

			my $peptide_filter = 'false';
			if ($job_asm_node->exists('peptide_filter')) { 
				my $peptide_filter_nodes = $job_asm_node->findnodes('peptide_filter');
				my @vers;
				foreach my $pf_node ($peptide_filter_nodes->get_nodelist() ) { 
					push (@vers, $pf_node->findvalue('@ver'));
				}
				@vers = sort {$b <=> $a} @vers;
				$peptide_filter_nodes = $job_asm_node->findnodes("peptide_filter[\@ver='$vers[0]']");
				my $peptide_filter_node;
				if ($peptide_filter_nodes->size() > 1) { 
					die "ERROR! multiple same version peptide filter nodes in $job_inx, $run_list_dir, weird...\n";
				}else{
					$peptide_filter_node = $peptide_filter_nodes->get_node(1);
				}
				$peptide_filter = $peptide_filter_node->findvalue('@apply');
			}
			my $coiled_coil_filter = 'false';
			if ($job_asm_node->exists('coiled_coil_filter')) { 
				my $cc_nodes = $job_asm_node->findnodes('coiled_coil_filter');
				my @vers;
				foreach my $cc_node ($cc_nodes->get_nodelist()) { 
					push (@vers, $cc_node->findvalue('@ver'));
				}
				@vers = sort {$b <=> $a} @vers;
				$cc_nodes = $job_asm_node->findnodes("coiled_coil_filter[\@ver='$vers[0]']");
				my $cc_node;
				if ($cc_nodes->size() > 1) { 
					die "ERROR! multiple same version cc filter nodes in $job_inx, $run_list_dir, weird...\n";
				}else{
					$cc_node = $cc_nodes->get_node(1);
				}
				$coiled_coil_filter = $cc_node->findvalue('@apply');
			}

			my $mc_domain_component = 'false';
			if ($job_asm_node->findvalue('@sw_mc_complex') eq 'true') { 
				$mc_domain_component = 'true';
			}


			#Generate some summary nodes

			my $ds_pdb_chain_list_XPath = '//domain_parse_summary/pdb_chain_list';	
			my ($ds_pdb_chain_list_node, $ds_pdb_chain_node, $ds_job_asm_node, $ds_domain_parse_node);
			if ($domain_summary_doc->exists($ds_pdb_chain_list_XPath)) { 
				$ds_pdb_chain_list_node = $domain_summary_doc->findnodes($ds_pdb_chain_list_XPath)->get_node(1);

				my $ds_pdb_chain_XPath = "pdb_chain[\@pdb=\'$query_pdb\'][\@chains=\'$query_chains\']";

				if ($ds_pdb_chain_list_node->exists($ds_pdb_chain_XPath)) { 
					my $ds_pdb_chain_nodes = $ds_pdb_chain_list_node->findnodes($ds_pdb_chain_XPath);
					if ($ds_pdb_chain_nodes->size() > 1) { 
						die "ERROR! $sub: More than one pdb_chain node for $query_pdb, $query_chains... $run_list_dir\n";
					}else{
						$ds_pdb_chain_node = $ds_pdb_chain_nodes->get_node(1);
					}
				}else{
					$ds_pdb_chain_node = $domain_summary_doc->createElement('pdb_chains');
		
					$ds_pdb_chain_node->setAttribute('pdb', $query_pdb);
					$ds_pdb_chain_node->setAttribute('chains', $query_chains);
					$ds_pdb_chain_node->setAttribute('week_label', $rl_week_label);
					$ds_pdb_chain_node->setAttribute('peptide_filter', $peptide_filter);
					$ds_pdb_chain_node->setAttribute('coiled_coil_filter', $coiled_coil_filter);
					$ds_pdb_chain_node->setAttribute('mc_domain_component', $mc_domain_component);
					$ds_pdb_chain_node->setAttribute('rep95', $rep95);

					$ds_pdb_chain_list_node->appendChild($ds_pdb_chain_node);
				}

				my $ds_domain_parse_XPath = qq{domain_parse[\@type='parser_asm']/job[\@mode='$job_mode'][\@reference='$reference\']};
				if ($ds_pdb_chain_node->exists($ds_domain_parse_XPath)) { 
					die "WARNING! $sub: domain parse for $query_pdb,$query_chains,$job_mode,$reference,parser already exists?\n";
				}else{
					$ds_job_asm_node = $domain_summary_doc->createElement('job_asm');
					$ds_domain_parse_node = $domain_summary_doc->createElement('domain_parse');

					$ds_job_asm_node->setAttribute('mode', $job_mode);
					$ds_job_asm_node->setAttribute('reference', $reference);
					#my $type 	= "parser_$mode";
					my $type	= $mode;
					$ds_domain_parse_node->setAttribute('type', $type);

					$ds_pdb_chain_node->appendChild($ds_domain_parse_node);
					$ds_domain_parse_node->appendChild($ds_job_asm_node);
				}

			}
			

			
			my $domain_xml_file;
			#This is not the best way to do this.
			if ($mode eq 'seq_iter') { 
				$domain_xml_file	= "$job_dump_dir/$pdb_chains/$domain_prefix.$pdb_chains.$reference.xml";
			}elsif($mode eq 'struct_search') { 
				$domain_xml_file	= "$job_dump_dir/$pdb_chains/$domain_prefix.$pdb_chains.$reference.struct_search.xml";
			}

			if (!-f $domain_xml_file) { 
				print "WARNING! $sub: Domain XML file $domain_xml_file not found!\n";
				$bad_domain_xml_file_count++;
				next;
			}

			if (-s $domain_xml_file == 0) { 
				print "WARNING! DOMAIN XML file $domain_xml_file empty!\n";
			}

			my $domain_xml_doc = xml_open($domain_xml_file);
			
			my $domain_doc_XPath	= "//chain_domains_set_top[\@domain_prefix='$domain_prefix']";
			if ($domain_xml_doc->exists($domain_doc_XPath)) { 
				my $domain_doc_nodes = $domain_xml_doc->findnodes($domain_doc_XPath);

				my $domain_doc_node;
				if ($domain_doc_nodes->size() > 1) { 
					die "ERROR! $sub: No support for multiple doc nodes per domain file in $domain_xml_file, aborting...\n";
				}else{
					$domain_doc_node = $domain_doc_nodes->get_node(1);
				}

				my $pdb_id	= lc($domain_doc_node->findvalue('@pdb_id'));
				my $chain_str	= $domain_doc_node->findvalue('@chain_str');
				my $domain_mode = $domain_doc_node->findvalue('@mode');

				#Sanity
				if (lc($pdb_id) ne lc($query_pdb)) { 
					print "WARNING $sub! PDB I/O mismatch for $pdb_id/$query_pdb \n";
				}

				if ($chain_str ne $query_chains){ 
					print "WARNING $sub! CHAINS I/O mismatch for $chain_str/$query_chains";
				}

				if ($job_mode ne $domain_mode && $domain_mode ne 'blast') { 
					print "WARNING $sub! MODE I/O mismatch for $job_mode/$domain_mode for $job_inx\n";
				}

				#my $coverage	= $domain_doc_node->findvalue('chain_domain_coverage');
				if (!$domain_doc_node->exists('chain_domain_coverage')) { 
					die "ERROR! $sub: No coverage node for $pdb_id $chain_str\n";
				}
				my $coverage_node = $domain_doc_node->findnodes('chain_domain_coverage')->get_node(1);
				if (!$ds_domain_parse_node->exists('chain_domain_coverage')) { 
					$ds_domain_parse_node->appendChild($coverage_node);
				}else{
					die "ERROR! $sub: Chain coverage node already exists for $job_inx, $query_pdb, $query_chains, $run_list_dir?\n";
				}

				if (!$domain_doc_node->exists('optimized_chain_domain_coverage')) { 
					print "WARNING! $sub: No optimized coverage node for $pdb_id $chain_str\n";
				}else{
					my $optimized_coverage_node = $domain_doc_node->findnodes('optimized_chain_domain_coverage')->get_node(1);
					if (!$ds_domain_parse_node->exists('optimized_chain_domain_coverage')) { 
						$ds_domain_parse_node->appendChild($optimized_coverage_node);
					}else{
						die "ERROR! $sub: Optimized chain coverage node already exists for $job_inx, $query_pdb, $query_chains, $run_list_dir?\n";
					}
				}

				my $domain_list_XPath = 'domain_list/domain';
				my $domain_list_nodes = $domain_doc_node->findnodes($domain_list_XPath);

				my $ds_domain_list_node;
				if (!$ds_job_asm_node->exists('domain_list')) { 
					$ds_domain_list_node = $domain_summary_doc->createElement('domain_list');
					$ds_job_asm_node->appendChild($ds_domain_list_node);
				}else{
					die "ERROR! $sub: Domain list already exists for $job_inx, $query_pdb, $query_chains, $run_list_dir?\n";
				}

				foreach my $domain_node ($domain_list_nodes->get_nodelist() ) { 

					my $domain_id 		= $domain_node->findvalue('@domain_id');
					my $hit_domain_id	= $domain_node->findvalue('hit_domain/@ecod_domain_id');	
					my $method		= $domain_node->findvalue('method');

					
					if ($DEBUG > 2) { 
						print "DEBUG: $domain_id $hit_domain_id $method\n";
					}	

					if ($method !~ /self_comp/ && !$domain_node->exists('hit_domain/domain')) { 
						if ($$reference_doc{$hit_domain_id}) { 
							#my $hit_domain_node	= $domain_node->findnodes('hit_domain')->get_node(1);
							#my $ref_hit_domain_node	= $reference_doc->findnodes(qq{//domain[\@ecod_domain_id='$hit_domain_id']})->get_node(1);
							#$hit_domain_node->appendChild($ref_hit_domain_node);
							#my $new_ref_node = $ref_hit_domain_node->cloneNode(1);
							#$hit_domain_node->appendChild($new_ref_node);
						}else{
							print "WARNING! $hit_domain_id not found in reference?  $reference $domain_id ($rl_week_label)\n";
						}
					}
					
					$ds_domain_list_node->appendChild($domain_node);
				}
			}else{
				print "ERROR! $sub: Fouled XPath/broken domain XML for top node on $domain_xml_file\n";
			}
				
		}
		
	}else{
		print "ERROR! $sub: Fouled XPath/broken job XML for top node on $rl_job_file\n";
		next;
	}
	if ($bad_domain_xml_file_count > 0) { 
		print "WARNING! $rl_job_file has $bad_domain_xml_file_count missing domain files!\n";
	}
	return 0;
}


sub register_repair { 
	my $sub = 'register_repair';

	my ($ecod_master_ref_inx_fn, $ecod_master_repair_inx_fn, $run_list_file) = @_;

	if (!-f $ecod_master_ref_inx_fn) { 
		die "ERROR! $sub: ECOD master reference index not found, $ecod_master_ref_inx_fn\n";
	}
	if (!-f $ecod_master_repair_inx_fn) { 
		die "ERROR! $sub ECOD master repair run index not found, $ecod_master_repair_inx_fn\n";
	}

	my $ecod_master_ref_xml = xml_open($ecod_master_ref_inx_fn);
	my $ecod_master_repair_xml = xml_open($ecod_master_repair_inx_fn);

	my $current_version = $ecod_master_ref_xml->findvalue('//@currentVersion');

	my $run_list_label = "run_list.repair.$current_version";

	my $abs_run_list_file = File::Spec->rel2abs($run_list_file);
	load_partition_vars();
	#my $domain_prefix = 'domains_v12'; #Fix this
	my $domain_prefix = $DOMAIN_PREFIX; #Fix this
	my $reference 	= $current_version;
	my $mode	= 'merge';

	my $run_list_node = $ecod_master_repair_xml->createElement('run_list');
	$run_list_node->appendTextChild('run_list_file', $abs_run_list_file);
	$run_list_node->appendTextChild('run_list_label', $run_list_label);
	$run_list_node->appendTextChild('domain_prefix', $domain_prefix);
	$run_list_node->appendTextChild('reference', $current_version);
	$run_list_node->appendTextChild('mode', 'merge');

	$run_list_node->setAttribute('display', 'true');
	$run_list_node->setAttribute('process', 'true');

	my $run_list_summary_node;
	if ($ecod_master_repair_xml->exists('//run_list_summary')) { 
		$run_list_summary_node = $ecod_master_repair_xml->findnodes('//run_list_summary')->get_node(1);
	}else{
		die "ERROR! $sub: No run list summary node in master repair\n";
	}

	my $i = 1;
	foreach my $run_list ($ecod_master_repair_xml->findnodes('//run_list')) { 
		my $old_run_list_label = $run_list->findvalue('run_list_label');
		if ($old_run_list_label eq $run_list_label) { 
			warn "WARNING! $run_list_label already in master repair lst\n";
			return 0;
		}
		$run_list->setAttribute('process', 'false');
		$i++
	}

	$run_list_node->setAttribute('id', $i++);
	$run_list_summary_node->appendChild($run_list_node);
	
	xml_write($ecod_master_repair_xml, $ecod_master_repair_inx_fn);

}
sub buildali_hhblits { 
	my $sub = 'buildali_hhblits';

	my ($job_xml_fn, $use_reps) = @_;

	print "DEBUG: $sub $job_xml_fn\n";

	my $job_xml_doc = xml_open($job_xml_fn);

	my ($job_dump_dir, $job_list_dir) = get_job_dirs_from_job_xml($job_xml_doc);
	my $job_list_nodes = find_job_nodes($job_xml_doc);

	my $reps = 0;
	if ($job_xml_doc->exists('//@rep95')) { $reps = 1 } 

	my @job_ids;
	foreach my $node ($job_list_nodes->get_nodelist() ) { 

		if ($use_reps && $reps && $node->findvalue('@rep95') ne 'true') { next } 

		#if ($node->findvalue('@poor') ne 'true') { next } 

		my $id 		= $node->findvalue('@id');
		my $mode	= $node->findvalue('mode');
		my $reference	= $node->findvalue('reference');

		my $query_pdb	= $node->findvalue('query_pdb');
		my $query_lc_pdb = lc($query_pdb);
		my $query_chain = $node->findvalue('query_chain');
		my $pdb_chain = $query_pdb . "_" . $query_chain;

		if ($node->findvalue('@type') eq 'recurse') { 
			my $seqid_range = $node->findvalue('seqid_range');
			my $clean_seqid_range = $seqid_range;
			$clean_seqid_range =~ s/\,/_/g;
			$pdb_chain .= "_$clean_seqid_range";
		}
		print "i: $id $pdb_chain\n";
		
		my $fa_file = "$job_dump_dir/$pdb_chain/$pdb_chain.fa";

		if (!-f $fa_file) { 
			print "WARNING! No FA file ($fa_file) for $query_pdb $query_chain\n";
			next;
		}

		my $a3m_file = "$job_dump_dir/$pdb_chain/$pdb_chain.a3m";
		my $new_a3m = $a3m_file;
		$new_a3m =~ s/a3m/remake.a3m/;

		my $new_hhm = $new_a3m;
		$new_hhm =~ s/a3m/hhm/;
		
		my $job_fn = "$job_dump_dir/$pdb_chain/hh_blits.$pdb_chain.job";
		my @jobs;
		

		if (-f $fa_file && (! -f $new_a3m || $FORCE_HH_OVERWRITE)) { 
			my $hh_lib = "export HHLIB=$HH_LIB\n";
			push (@jobs, $hh_lib);
			my $hhblits_command = "$HHBLITS_EXE -i $fa_file -oa3m $new_a3m -ohhm $new_hhm -addss -psipred $PSIPRED_EXE -psipred_data $PSIPRED_DATA -d $NR20_TMP_DB -cpu 8 ";
			print "$hhblits_command\n";
			push (@jobs, $hhblits_command);
		}

#		if (! -f $new_hhm ) { 
#			my $hhmake_command = "$HHMAKE_EXE -i $new_a3m";
#			push (@jobs, $hhmake_command);
#		}

		if (scalar(@jobs) > 0) { 
			job_create($job_fn, \@jobs, 1);
			my $job_id = qsub($job_fn);
			if ($job_id =~ /\d+/) { 
				push (@job_ids, $job_id);
			}else{
				croak "?$job_id\n";
			}

		}
	}
	return \@job_ids;
}
sub bsumm { 
	my $sub = 'bsumm';

	my ($job_xml_fn, $use_reps, $force_overwrite) = @_;

	my $xml_doc = xml_open($job_xml_fn);

	my ($job_dump_dir, $job_list_dir) = get_job_dirs_from_job_xml($xml_doc);
	my $job_nodes = find_job_nodes($xml_doc);

	#my $bsumm_script = '/data/ecod/weekly_updates/weeks/bin/single_domain_summary.pl';

	my $reps = 0;
	if ($xml_doc->exists('//@rep95')) { 
		$reps = 1;
	}

	my @job_ids;
	foreach my $node ($job_nodes->get_nodelist() ) { 

		my $mode = $node->findvalue('mode');
		if ($mode ne 'seq_iter') { next } 

		my $bsumm_script;
		if (!$use_reps || $reps &&  $node->findvalue('@rep95') eq 'true') { 
			$bsumm_script = '/data/ecod/weekly_updates/weeks/bin/single_domain_summary.pl';
		}else{
			$bsumm_script = '/data/ecod/weekly_updates/weeks/bin/single_domain_summary.pl --concise';
		}

		my ($query_pdb, $query_chain) = get_pdb_chain($node);

		my $query_lc_pdb = lc($query_pdb);
		my $pdb_chain = "${query_pdb}_${query_chain}";
		my $reference = $node->findvalue('reference');

		my $job_fn = "$job_dump_dir/$pdb_chain/bsumm.$pdb_chain.$reference.job";
		my $bsumm_cmd = "$bsumm_script $query_pdb $query_chain $reference $job_dump_dir";

		my @jobs;
		push (@jobs, $bsumm_cmd);
		job_create($job_fn, \@jobs);
		my $job_id = qsub($job_fn);
		print "bsumm $job_id $job_fn\n";
		push (@job_ids, $job_id);
	}

	return \@job_ids;

}
sub hh_run { 
	my $sub = 'hh_run';

	my ($job_xml_fn, $reference, $use_reps) = @_;

	my %seqid_range;
	my %ecod_domain_id;
	my %uid_lookup;
	load_references();

	my $ref_xml_doc = xml_open($REF_XML{$reference});	
	
	foreach my $ref_node (find_manual_domain_nodes($ref_xml_doc)) { 

		my ($uid, $ecod_domain_id) 	= get_ids($ref_node);
		my $seqid_range			= get_seqid_range($ref_node);

		$seqid_range{$uid}		= $seqid_range;
		$ecod_domain_id{$uid}		= $ecod_domain_id;
		$uid_lookup{$ecod_domain_id}	= $uid;

	}

	my $job_xml_doc	= xml_open($job_xml_fn);

	my ($job_dump_dir, $job_list_dir) = get_job_dirs_from_job_xml($job_xml_doc);
	my $job_nodes = find_job_nodes($job_xml_doc);

	my $reps = 0;
	if ($job_xml_doc->exists('//@rep95')) { 
		$reps = 1;
	}
	my @job_ids;
	foreach my $node ($job_nodes->get_nodelist() ) { 

		if ($reps && $use_reps && $node->findvalue('@rep95') ne 'true') { next } 

		my $query_pdb	= $node->findvalue('query_pdb');
		my $query_chain	= $node->findvalue('query_chain');
		my $pdb_chain = "${query_pdb}_${query_chain}";

		my $input = 'struct_seqid';
		if ($node->exists('mode/@input')) { 
			$input 	= $node->findvalue('mode/@input');
		}

		my $recurse = 0;
		my $recurse_range;
		if ($node->findvalue('@type') eq 'recurse') { 
			my $seqid_range = $node->findvalue('seqid_range');
			$recurse_range = $seqid_range;
			my $clean_seqid_range = $seqid_range;	
			$clean_seqid_range =~ s/\,/_/g;
			$pdb_chain .= "_$clean_seqid_range";
		}

		my $query_reference = $node->findvalue('reference');
		
		if ($query_reference ne $reference) { 
			print "WARNING! query reference $query_reference does not match input reference $reference\n";
			next;
		}

		my $ecod_dir = "$job_dump_dir/$pdb_chain";

		#my $a3m_fn = "$pdb_chain.a3m";
		my $a3m_fn = "$pdb_chain.remake.a3m";
		if (!-f "$ecod_dir/$a3m_fn") { 
			print "WARNING! a3m file $ecod_dir/$a3m_fn not found\n";
			next;
		}

		my $hhm_fn = $a3m_fn;
		$hhm_fn =~ s/a3m/hhm/;

		my $job_fn = "$job_dump_dir/$pdb_chain/hh.$pdb_chain.job";
		my @jobs;
		my $hh_lib = "export HHLIB=$HH_LIB\n";
		push (@jobs, $hh_lib);
 
#		if (!-f $hhm_fn || $FORCE_HH_OVERWRITE) { 
#			if ($DEBUG) { 
#				print "HHmake $a3m_fn\n";
#			}
#			
#			my $hhmake_command =  "$HHMAKE_EXE -i $a3m_fn";
#			push (@jobs, $hhmake_command);
#		}

		#my $result_fn = "$pdb_chain.$reference.result";	
		my $result_fn = "$ecod_dir/$pdb_chain.$reference.rebuild.result";	

		if (!-f $result_fn || $FORCE_HH_OVERWRITE) { 
			if ($DEBUG) { 
				print "HHsearch $hhm_fn\n";
			}

			my $hh_cmd = "$HHSEARCH_EXE -i $ecod_dir/$hhm_fn -o $result_fn -d $HH_REF{$reference} -cpu 8\n";
			push (@jobs, $hh_cmd);
			job_create($job_fn, \@jobs);
			my $job_id = qsub($job_fn);
			push (@job_ids, $job_id);
		}
	}
	return (\@job_ids);
}
sub hh_parse { 
	my $sub = 'hh_parse';

	my ($job_xml_fn, $reference, $use_reps) = @_;
	
	my %seqid_range;
	my %ecod_domain_id;
	my %uid_lookup;
	
	my $ref_xml_doc = xml_open($REF_XML{$reference});

	foreach my $ref_node (find_manual_domain_nodes($ref_xml_doc)) { 

		my ($uid, $ecod_domain_id)	= get_ids($ref_node);
		my $seqid_range			= get_seqid_range($ref_node);

		$seqid_range{$uid}		= $seqid_range;
		$ecod_domain_id{$uid}		= $ecod_domain_id;
		$uid_lookup{$ecod_domain_id}	= $uid;

	}

	my $job_xml_doc = xml_open($job_xml_fn);


	my ($job_dump_dir, $job_list_dir) = get_job_dirs_from_job_xml($job_xml_doc);

	my $job_nodes = find_job_nodes($job_xml_doc);

	my $reps = 0;
	if ($job_xml_doc->exists('//@rep95')) { 
		$reps = 1;
	}
	my @job_ids;
	foreach my $node ($job_nodes->get_nodelist() ) { 

		if ($reps && $use_reps &&  $node->findvalue('@rep95') ne 'true') { next } 

	
		my ($query_pdb, $query_chain) = get_pdb_chain($node);
		my $pdb_chain = $query_pdb . "_" . $query_chain;

		my $input = 'struct_seqid';
		if ($node->exists('mode/@input')) { 
			$input 	= $node->findvalue('mode/@input');
		}

		my $recurse = 0;
		my $recurse_range;
		if ($node->findvalue('@type') eq 'recurse') { 
			my $seqid_range 	= $node->findvalue('seqid_range');
			$recurse_range 		= $seqid_range;
			my $clean_seqid_range 	= $seqid_range;	
			$clean_seqid_range 	=~ s/\,/_/g;
			$pdb_chain 		.= "_$clean_seqid_range";
		}

		my $query_reference = $node->findvalue('reference');
		
		if ($query_reference ne $reference) { 
			print "WARNING! query reference $query_reference does not match input reference $reference\n";
			next;
		}

		my $ecod_dir = "$job_dump_dir/$pdb_chain";

		my $a3m_fn = "$pdb_chain.remake.a3m";
		if (!-f "$ecod_dir/$a3m_fn") { 
			print "WARNING! a3m file $ecod_dir/$a3m_fn not found\n";
			next;
		}

		my $hhm_fn = $a3m_fn;
		$hhm_fn =~ s/a3m/hhm/;

		my $job_fn = "$job_dump_dir/$pdb_chain/hh_parse.$pdb_chain.job";
		my @jobs;
		
		my $result_fn 	= "$ecod_dir/$pdb_chain.$reference.rebuild.result";	
		my $summ_fn 	= "$ecod_dir/$pdb_chain.$reference.rebuild.hh_summ.xml";	

		my $hh_parse_script = "/data/ecod/weekly_updates/weeks/bin/hh_parse.pl";

		if (!-f $summ_fn || $FORCE_HH_OVERWRITE) { 
			if ($DEBUG) { 
				print "HHparse $hhm_fn\n";
			}
	
			my $input = 'struct_seqid';
			my $hh_cmd = "$hh_parse_script $query_pdb $query_chain $result_fn $reference $summ_fn $input";
			#print "$hh_cmd\n";

			push (@jobs, $hh_cmd);
			job_create($job_fn, \@jobs);
			my $job_id = qsub($job_fn);
			push (@job_ids, $job_id);
			
		}
	}
	return (\@job_ids);
}
sub optimize { 
	my $sub = 'optimize';

	my ($job_xml_fn, $force_overwrite) = @_;

	my $job_xml_doc	= xml_open($job_xml_fn);

	my ($job_dump_dir, $job_list_dir) = get_job_dirs_from_job_xml($job_xml_doc);
	my $job_nodes = find_job_nodes($job_xml_doc);

	my $reps = hasReps_job_list($job_xml_doc);
	
	load_partition_vars();
	my $domain_prefix = $DOMAIN_PREFIX;

	my @optimize_job_ids;
	if ($job_nodes->size() == 0) { 
		print "WARNING! : No jobs nodes found for $job_xml_fn, skipping...\n";
		return 0;
	}else{
		foreach my $job_node (find_job_nodes($job_xml_doc)) { 

			my ($pdb, $chain) = get_pdb_chain($job_node);
			my $mode	= $job_node->findvalue('mode');
			
			my $pdb_chain	= $pdb . "_" . $chain;

			my $reference	= $job_node->findvalue('reference');

			my $recurse = 0;
			my $recurse_range;
			if ($job_node->findvalue('@type') eq 'recurse') { 
				my $seqid_range = $job_node->findvalue('seqid_range');
				$recurse_range = $seqid_range;
				$seqid_range =~ s/\,/_/;
				$pdb_chain .= "_$seqid_range";
				$recurse++;
			}

			my $domain_fn;
			my $domain_dir;
			if ($mode eq 'seq_iter') { 
				$domain_dir	 = "$job_dump_dir/$pdb_chain";
				$domain_fn	 = "$domain_prefix.$pdb_chain.$reference.xml";
			}elsif($mode eq 'struct_search') { 
				$domain_dir	= "$job_dump_dir/$pdb_chain";
				$domain_fn 	= "$domain_prefix.$pdb_chain.$reference.struct_search.xml";
			}else{
				die "ERROR! : Unknown run mode $mode\n";
			}

			if (! -f "$domain_dir/$domain_fn") { 
				print "WARNING! Could not find domain file $domain_dir $domain_fn, skipping...\n";
				next;
			}else{
				#boundary_optimize($pdb_chain, $domain_dir, $domain_fn, $recurse, $recurse_range);
				my $optim_cmd;
				if ($force_overwrite) { 
					print "$SINGLE_BOUNDARY_OPTIMIZE $pdb_chain $domain_dir $domain_fn --force_overwrite \n";
					$optim_cmd =  "$SINGLE_BOUNDARY_OPTIMIZE $pdb_chain $domain_dir $domain_fn --force_overwrite\n";
				}else{
					print "$SINGLE_BOUNDARY_OPTIMIZE $pdb_chain $domain_dir $domain_fn \n";
					$optim_cmd =  "$SINGLE_BOUNDARY_OPTIMIZE $pdb_chain $domain_dir $domain_fn \n";

				}
				my @jobs;
				push (@jobs, $optim_cmd);
				my $job_fn = "$job_dump_dir/$pdb_chain/optim.$pdb_chain.job";

				job_create($job_fn, \@jobs);

				my $job_id = qsub($job_fn);
				push (@optimize_job_ids, $job_id);
			}
		}
	}
	return \@optimize_job_ids;
}

sub peptide_filter { 
	my $sub = 'peptide_filter';
	my ($job_list_xml_fn) = @_;

	my $job_xml_doc 			= xml_open($job_list_xml_fn);
	my ($job_dump_dir, $job_list_dir) 	= get_job_dirs_from_job_xml($job_xml_doc);

	my @peptide_job_ids;
	foreach my $job_node (find_job_nodes($job_xml_doc)) { 
		my ($query_pdb, $query_chain) = split '_', job_node_pdb_chain($job_node);		 
		my $pc = $query_pdb . "_" . $query_chain;
		
		my $cmd = "$SINGLE_PEPTIDE_FILTER $query_pdb $query_chain $job_dump_dir/$pc/\n";
		my @jobs;
		push (@jobs, $cmd);
		my $job_fn = "$job_dump_dir/$pc/peptide.$pc.job";
		job_create($job_fn, \@jobs);

		my $job_id = qsub($job_fn);
		push (@peptide_job_ids, $job_id);

	}

	return \@peptide_job_ids;


}
sub peptide_collate { 
	my $sub = 'peptide_collate';
	my ($job_list_xml_fn) = @_;

	print "peptide collate\n";
	my $job_xml_doc 			= xml_open($job_list_xml_fn);
	my ($job_dump_dir, $job_list_dir) 	= get_job_dirs_from_job_xml($job_xml_doc);
	
	foreach my $job_node (find_job_nodes($job_xml_doc)) { 
		my $pc = job_node_pdb_chain($job_node);		 
		#my ($query_pdb, $query_chain) = job_node_pdb_chain($job_node);		 
		#my $pc = $query_pdb . "_" . $query_chain;
		
		my $path = "$job_dump_dir/$pc/$pc.peptide.v$FILTER_VERSION.xml";
		if (-f $path) { 
			my $peptide_xml_doc = xml_open($path);
			if ($peptide_xml_doc->exists('//peptide_filter') && !$job_node->exists('peptide_filter')) { 
				$job_node->appendChild($peptide_xml_doc->findnodes('//peptide_filter')->get_node(1));
			}
		}
	}
	my $out_fn = "$job_list_xml_fn.peptide";
	xml_write($job_xml_doc, $out_fn);
	if (-f $out_fn) { 
		move($out_fn, $job_list_xml_fn);
	}	
}


sub domains { 
	my $sub = 'domains';
	my ($job_list_xml_fn, $reference, $use_reps, $force_overwrite) = @_;
	#MOVE executable locations
	#my $single_domain_partition = '/data/ecod/weekly_updates/weeks/bin/single_domain_partition_v12.pl';
	my $single_domain_partition = '/data/ecod/weekly_updates/weeks/bin/single_domain_partition_v13.pl'; #Experimental
	if (!-f $single_domain_partition) { 
		die "ERROR! Could not find single domain partition script $single_domain_partition\n";
	}
	load_references();
	load_partition_vars();

	my $job_xml_doc = xml_open($job_list_xml_fn);
	my ($job_dump_dir, $job_list_dir) = get_job_dirs_from_job_xml($job_xml_doc);
	if (!-d $job_dump_dir) { 
		print "Could not find job dump dir $job_dump_dir, aborting...\n";
	}
	my $job_list_nodes = find_job_nodes($job_xml_doc);

	my $global_reference = $reference ? $reference : $LATEST_REFERENCE;
	
	load_partition_vars();
	my $domain_prefix = $DOMAIN_PREFIX;
	#Ref range, ecod_domain_id, $uid, $range
	my ($ref_range_cache, $ref_domain_uid_lookup) = reference_cache_load($global_reference);
	my $ref_chain_domains = reference_chainwise_transform($ref_range_cache);

	my $reps = hasReps_job_list($job_xml_doc);
	
	my @bdom_job_ids;
	foreach my $node (find_job_nodes($job_xml_doc)) { 

		my $concise = 0;

		if ($node->findvalue('@poor') ne 'true' && $POOR_ONLY) { next }  
		if ($reps && $use_reps & $node->findvalue('@rep95') ne 'true') {
			if ($PERFORM_NREP_CONCISE) { 
				$concise = 1;
			}else{ 
				next;
			}
		} 

		my $id 		= $node->findvalue('@id');
		my $mode	= $node->findvalue('mode'); #Jobs have either seq_iter or struct_search mode
		if ($mode ne 'seq_iter') { next } 


		my $input_mode = $node->exists('mode/@input') ? $node->findvalue('mode/@input') : 'struct_seqid';

		#Right now we assume all jobs in a run have the same reference, this may need to be broken
		my $reference	= $node->findvalue('reference');
		if ($reference ne $global_reference) { 
			print "WARNING! $reference != $global_reference\n";
			$ref_range_cache = reference_cache_load($global_reference);
		}


		my ($query_pdb, $query_chain) = get_pdb_chain($node);
		my $query_lc_pdb = lc($query_pdb);
		my $pdb_chain = "${query_pdb}_${query_chain}";
		print "pc: $pdb_chain\n";

#Recursion has been turned off for the moment. But domain partition should take a range argument in the future
#	my $recurse = 0;
#	my $recurse_range;
#	if ($node->findvalue('@type') eq 'recurse') { 
#		my $seqid_range = $node->findvalue('seqid_range');
#		$recurse_range = $seqid_range;
#		$seqid_range =~ s/\,/_/g;
#		$pdb_chain .= "_$seqid_range";
#		$recurse++;
#	}

		my $node_dir = "$job_dump_dir/$pdb_chain";

		if (!-d $node_dir) { 
			print "WARNING! Node dir ($node_dir) not found for $query_pdb, $query_chain! Skipping...\n";
			next;
		}

		my $domain_fn 	= "$node_dir/$domain_prefix.$pdb_chain.$reference.xml";

		my $job_fn 	= "$node_dir/bdom.$pdb_chain.job";
		my @job_lns;

		if (-f $domain_fn && !$force_overwrite) { 
			print "SKIPPING $pdb_chain\n";
			next ; 
		}else{
			print "$domain_fn\n";
		}
		my $job_cmd = "$single_domain_partition $query_pdb $query_chain $node_dir $input_mode $reference $concise single\n";
		push (@job_lns, $job_cmd);
		job_create($job_fn, \@job_lns);
		my $job_id = qsub($job_fn);
		push (@bdom_job_ids, $job_id);
	}
	return \@bdom_job_ids;
}
sub chblastp_job_file { 
	my $sub = 'chblastp_job_file';
	
	my ($job_xml_fn, $reps_only, $force_overwrite) = @_;
	if ($reps_only) { 
		print "DEBUG: $sub $reps_only reps_only\n";
	}
	load_references();

	my $job_xml_doc = xml_open($job_xml_fn);

	my @job_ids; 

	my ($job_dump_dir, $job_list_dir) = get_job_dirs_from_job_xml($job_xml_doc);

	if (!-d $job_dump_dir) { 
		die "Could not find job_dump dir $job_dump_dir, aborting...\n";
	}

	my $job_list_nodes = find_job_nodes($job_xml_doc);

	my $hn = `hostname`;

	my $BLAST_CMD;
	if ($hn =~ /lotta/) { 
		$BLAST_CMD = '/usr1/ncbi-blast-2.2.25+/bin/blastp';
	}else{
		$BLAST_CMD = '/usr1/ncbi-blast-2.2.25+-x64/bin/blastp';
	}


	foreach my $node ($job_list_nodes->get_nodelist() ) { 
		my $id 		= $node->findvalue('@id');
		#my $hora_id	= $node->findvalue('@hora_id');
		
		print "$id \n";
		#if ($node->findvalue('@poor') ne 'true') { next }

		my $mode	= $node->findvalue('mode');
		my $reference	= $node->findvalue('reference');

		if ($mode ne 'seq_iter') { next } 

		my ($pdb, $chain) = get_pdb_chain($node);
		my $pdb_chain = "${pdb}_${chain}";

		my $db;
		if ($reps_only) { 
			if ($CHAIN_REPSONLY_REF{$reference}) { 
				$db = $CHAIN_REPSONLY_REF{$reference};
				#print "DEBUG: chain_repsonly_ref $CHAIN_REPSONLY_REF{$reference}\n";
			}else{
				die "ERROR! No chain repsonly for $reference\n";
			}
		}else{
			$db = $CHAIN_REF{$reference};
		}


		if ($node->findvalue('@type') eq 'recurse') { 
			my $seqid_range = $node->findvalue('seqid_range');
			$seqid_range =~ s/\,/_/g;
			$pdb_chain .= "_$seqid_range";
		}

		if (!-f "$job_dump_dir/$pdb_chain/$pdb_chain.fa") { 
			print "WARNING: No fasta file found for $id $pdb_chain\n";
			next;
		}
		my $blast_fn = "$job_dump_dir/$pdb_chain/$pdb_chain.$reference.chainwise_blast.xml";
		if (-f $blast_fn && !$force_overwrite && -s $blast_fn > 0) { next } 

		#Run Blast
		my $blast_cmd =  "$BLAST_CMD -query $job_dump_dir/$pdb_chain/$pdb_chain.fa  -db $db -outfmt 5 -num_alignments 5000 -evalue 0.002 > $job_dump_dir/$pdb_chain/$pdb_chain.$reference.chainwise_blast.xml\n";
		#print `$BLAST_CMD -query $job_dump_dir/$pdb_chain/$pdb_chain.fa  -db $CHAIN_REF{$reference} -outfmt 5 -num_alignments 5000 -evalue 0.002 > $job_list_dump_dir/$pdb_chain/$pdb_chain.$reference.chainwise_blast.xml`;

		my $job_fn = "$job_dump_dir/$pdb_chain/chblastp.$pdb_chain.$reference.job";

		my $immediate = 0;

		if ($immediate) { 
			print `$blast_cmd`;
		}else{
			my @jobs;
			push (@jobs, $blast_cmd);
			job_create($job_fn, \@jobs);
			my $job_id = qsub ($job_fn);
			push (@job_ids, $job_id);
		}
	}

	return (\@job_ids);

}
sub blastp_job_file { 
	my $sub  = 'blastp_job_file';

	my ($job_xml_fn, $reps_only, $force_overwrite) = @_;
	print "DEBUG $sub: top $job_xml_fn $reps_only\n";

	my @job_ids;
	load_references();

	my $job_xml_doc = xml_open($job_xml_fn);

	my ($job_dump_dir, $job_list_dir) = get_job_dirs_from_job_xml($job_xml_doc);

	if (!-d $job_dump_dir) { 
		die "Could not find job_dump dir $job_dump_dir, aborting...\n";
	}

	my $BLAST_CMD = '/usr1/ncbi-blast-2.2.25+/bin/blastp';
	
	foreach my $node (find_job_nodes($job_xml_doc)) { 
		my $id 		= $node->findvalue('@id');
		
		print "$id \n";

		my $mode	= $node->findvalue('mode');
		my $reference	= $node->findvalue('reference');

		my $db;
		if ($reps_only) { 
			$db = $DOMAIN_REPSONLY_REF{$reference};
			#print "DEBUG: domain_repsonly_ref $DOMAIN_REPSONLY_REF{$reference}\n";
		}else{
			$db = $DOMAIN_REF{$reference};
		}


		if ($mode ne 'seq_iter') { next } 

		my ($pdb, $chain) = get_pdb_chain($node);
		my $pdb_chain = "${pdb}_${chain}";

		if ($node->findvalue('@type') eq 'recurse') { 
			my $seqid_range = $node->findvalue('seqid_range');
			$seqid_range =~ s/\,/_/g;
			$pdb_chain .= "_$seqid_range";
		}
		if (!-f "$job_dump_dir/$pdb_chain/$pdb_chain.fa") { 
			print "WARNING: No fasta file found for $id $pdb_chain\n";
			next;
		}

		#Run Blast
		my $job_fn = "$job_dump_dir/$pdb_chain/blast.$pdb_chain.$reference.job";

		my $immediate = 0;

		if ($force_overwrite || ! -f "$job_dump_dir/$pdb_chain/$pdb_chain.$reference.blast.xml" || -s "$job_dump_dir/$pdb_chain/$pdb_chain.$reference.blast.xml" == 0) { 
			my $blast_cmd =  "$BLAST_CMD -query $job_dump_dir/$pdb_chain/$pdb_chain.fa  -db $db -outfmt 5 -num_alignments 5000 -evalue 0.002 > $job_dump_dir/$pdb_chain/$pdb_chain.$reference.blast.xml";
			#print $blast_cmd;
			if ($immediate) { 
				print `$blast_cmd`;
			}else{
				my @jobs;
				push (@jobs, $blast_cmd);
				job_create($job_fn, \@jobs);
				my $job_id = qsub($job_fn);
				push (@job_ids, $job_id);
			}
		}
	}

	return (\@job_ids);
}
sub job_list_maintain { 
	my ($job_list_xml_fn, $new_week, $reference, $reps_only, $force_overwrite, $no_peptide) = @_;

	#Generate FASTA files.
	print "Generate FASTA...\n";

	my $job_xml_doc = xml_open($job_list_xml_fn);

	my ($job_dump_dir, $job_list_dir) = get_job_dirs_from_job_xml($job_xml_doc);

	load_references();
	load_partition_vars();

	print "#Build FASTA\n";
	foreach my $job_node (find_job_nodes($job_xml_doc)) { 

		my ($query_pdb, $query_chain) = get_pdb_chain($job_node);	
		my $pdb_chain	= $query_pdb . "_" . $query_chain;

		if (!-d "$job_dump_dir/$pdb_chain") { 
			if (!mkdir("$job_dump_dir/$pdb_chain")) { 
				die "ERROR! Could not create $job_dump_dir/$pdb_chain\n";
			}else{
				chown($UID, $GID, "$job_dump_dir/$pdb_chain");
			}
		}
		my $fasta_fn = "$job_dump_dir/$pdb_chain/$pdb_chain.fa";
	
		if (-f $fasta_fn && $force_overwrite < 2) { 
			print "WARNING! $fasta_fn exists, skipping...\n";
			next;
		}else{
			my ($seqid_aref, $struct_seqid_aref, $pdbnum_aref, $asym_id) 
				= pdbml_seq_parse($query_pdb, $query_chain);

			if (!$struct_seqid_aref || $struct_seqid_aref == 0 || scalar(@$struct_seqid_aref) == 0) {
				print "WARNING! No struct_seq for $query_pdb, $query_chain\n";
				$job_node->setAttribute('unstructured', 'true');
				next;
			}

			my $struct_seqid_range    = rangify(@$struct_seqid_aref);
			my $ungapped_struct_range = ungap_range($struct_seqid_range, $GAP_TOL);
			my $ungap_aref            = range_expand($ungapped_struct_range);

			if ($ungapped_struct_range eq 0) { 
				print "WARNING! empty range for $query_pdb, $query_chain\n";
				next;
			}

			my $fasta_string = pdbml_fasta_fetch($query_pdb, $asym_id, $query_chain, $ungap_aref);

			if (!$fasta_string) { 
				print "WARNING! No fasta string for $query_pdb, $query_chain\n";
				next;
			}

			open (my $fh, ">", $fasta_fn) or die "ERROR! Could not open $fasta_fn for writing:$!\n";
			print $fh ">$query_pdb,$query_chain\n$fasta_string\n";
			close $fh; 
		}
	}

	#Generate run list cluster
	if ($DEBUG) { 
		print "DEBUG: cluster: $new_week\n";
	}
	my $cluster_script = '/data/ecod/weekly_updates/weeks/bin/run_list_cluster.pl';
	if (!-f $cluster_script) { 
		die "run list cluster script not found! $cluster_script\n";
	}
	open (my $fh, ">", "cluster.$new_week.job") or die "Could not open cluster.$new_week.job for writing:$!\n";
	print $fh "#!/bin/bash\n";
	print $fh "#\$ -cwd\n";
	print $fh "#\$ -j y \n";
	print $fh "#\$ -S /bin/bash\n";
	print $fh "#\$ -M dustin.schaeffer\@gmail.com\n";
	print $fh "$cluster_script $job_list_xml_fn\n";
	close $fh;

	my $clust_job_id = qsub("cluster.$new_week.job");

	print "$clust_job_id\n";

	while (qstat_wait($clust_job_id)) { 
		print "SLEEPING cluster\n";
		sleep(10);
	}

	#Build profiles
	my $buildali_script = '/data/ecod/weekly_updates/weeks/bin/xml_hhblits_profile_build.pl';
	if (!-f $buildali_script) { 
		die "buildali script not found! $buildali_script\n";
	}

	if ($DEBUG) { 
		print "DEBUG: buildali_hhblits: $new_week\n";
	}
	my $ali_job_ids = buildali_hhblits($job_list_xml_fn);

	while (qstat_wait_list($ali_job_ids)) { 
		print "SLEEPING...\n";
		sleep(30); #5 min job checks;
	}

	unless($no_peptide) { 
		#Generate peptide filter
		if ($DEBUG) { 
			print "DEBUG: peptide_prefilter $new_week\n";
		}
		#my $peptide_filter_script = '/data/ecod/weekly_updates/weeks/bin/xml_peptide_prefilter.pl';
		#if (!-f $peptide_filter_script) { 
		#	die "peptide filter script not found! $peptide_filter_script\n";
		#}
		#open (OUT, ">pep.$new_week.job") or die "Coudl not open pep.$new_week.job for writing:$!\n";
#		print OUT "#!/bin/bash\n";
#		print OUT "#\$ -cwd\n";
#		print OUT "#\$ -j y \n";
#		print OUT "#\$ -S /bin/bash\n";
#		print OUT "#\$ -M dustin.schaeffer\@gmail.com\n";
#		print OUT "$peptide_filter_script $job_list_xml_fn\n";
#		close OUT;
#
		#my $pep_job_id = qsub("pep.$new_week.job");
		my $pep_job_ids = peptide_filter($job_list_xml_fn);
		
		while (qstat_wait_list($pep_job_ids)) { 
			print "SLEEPING peptide!\n";
			sleep(10);
		}
		peptide_collate($job_list_xml_fn);
	}

	#Run self_comp jobs
	if ($DEBUG) { 
		print "DEBUG: self_comp: $new_week\n";
	}
	#my $self_comp_script = '/home/rschaeff/work/domain_partition/ecod_domain_parser/bin/ecod_self_comparison.pl';
	my $self_comp_script = '/data/ecod/weekly_updates/weeks/bin/ecod_self_comparison.pl';
	if (!-f $self_comp_script) { 
		die "self comp script not found! $self_comp_script\n";
	}

	open ($fh, ">", "self_comp.$new_week.job") or die "Could not open self_comp.$new_week.job for writing:$!\n";
	print $fh "#!/bin/bash\n";
	print $fh "#\$ -cwd\n";
	print $fh "#\$ -j y \n";
	print $fh "#\$ -S /bin/bash\n";
	print $fh "#\$ -M dustin.schaeffer\@gmail.com\n";
	print $fh "$self_comp_script $job_list_xml_fn\n";
	close $fh;

	my $self_comp_job_id = qsub("$new_week.self_comp.job");

	#Run Chain-blast jobs
	if ($DEBUG ) { 
		print "DEBUG: chblastp: $new_week\n";
	}
	print "CHBLAST $job_list_xml_fn $reps_only\n";
	my $chblast_job_ids = chblastp_job_file($job_list_xml_fn, $reps_only, $force_overwrite);

	while (qstat_wait_list($chblast_job_ids)) { 
		print "SLEEPING...\n";
		sleep(30); #5 min job checks;
	}
	
	#Run blastp jobs 
	if ($DEBUG) { 
		print "DEBUG blastp: $new_week\n";
	}
	print "BLAST $job_list_xml_fn $reps_only\n";
	my $blast_job_ids = blastp_job_file($job_list_xml_fn, $reps_only, $force_overwrite);

	while(qstat_wait_list($blast_job_ids)) { 
		print "SLEEPING...\n";
		sleep(30);
	}

	#Run hhsearch jobs dependent on ali
	#Split hh search into run and parse components and allow for parallelization

	if ($DEBUG) { 
		print "DEBUG hh: $new_week\n";
	}

	my $hh_run_job_ids = hh_run($job_list_xml_fn, $reference, 1);

	printf "hhrun %i\n", scalar(@$hh_run_job_ids);
	while (qstat_wait_list($hh_run_job_ids)) { 
		print "SLEEPING...\n";
		sleep(30);
	}

	if ($DEBUG) { 
		print "DEBUG hh_parse: $new_week\n";
	}

	my $hh_parse_job_ids = hh_parse($job_list_xml_fn, $reference, 1);

	printf "hhparse %i\n", scalar(@$hh_parse_job_ids);
	while (qstat_wait_list($hh_parse_job_ids)) { 
		print "SLEEPING...\n";
		sleep(30);
	}
	

	#Run blast domain summary
	if ($DEBUG) { 
		print "DEBUG bsumm: $new_week\n";
	}

	my $bsumm_job_ids = bsumm($job_list_xml_fn, 1);
	printf "bsumm %i\n", scalar(@$bsumm_job_ids);
	while (qstat_wait_list($bsumm_job_ids)) { 
		print "SLEEPING...\n";
		sleep(30);
	}

	#Run seq_iter domain partition

	if ($DEBUG) { 
		print "DEBUG domains: $new_week\n";
	}

	my $domains_job_ids = domains($job_list_xml_fn, $reference, 1, $FORCE_OVERWRITE);
	while (qstat_wait_list($domains_job_ids)) { 
		print "SLEEPING...\n";
		sleep(30);
	}
#
	#Run boundary optimize
	my $optimize_job_ids = optimize($job_list_xml_fn);
	while (qstat_wait_list($optimize_job_ids)) { 
		print "SLEEPING...\n";
		sleep(30);
	}

#
#	my $DALI_SEARCH = 1;
#	if ($DALI_SEARCH) { 
#		foreach my $run_node ($run_list_xml_doc->findnodes('//domain_parse_run')) { 
#			if ($run_node->exists('run_list_job_xml_file[@mode="struct_search"]')) { 
#
#				my $run_list_dir 	= $run_node->findvalue('run_list_dir');
#				my $job_list_xml_fn 	= $run_node->findvalue('run_list_job_xml_file[@mode="struct_search"]');
#
#				my $job_list_xml_path = "$run_list_dir/$job_list_xml_fn";
#				print "jlxf: $job_list_xml_fn\n";			
#				print "jlxp: $job_list_xml_path\n";
#				struct_search_dali_query_gen($job_list_xml_path);
#				my $dali_jobs_aref = query_glob($job_list_xml_path);
#				throttled_qsub($dali_jobs_aref, 20000);
#				dali_summary($job_list_xml_path, $LATEST_REFERENCE);
#			}
#		}
#	}
}

#Think about moving this to Domains::Dali
sub dali_summary { 
	my $sub = 'dali_summary';
	my ($job_list_xml_fn,  $reference, $force_overwrite) = @_;

	load_references();
	my $job_xml_doc 	= xml_open($job_list_xml_fn);
	my $ecod_xml_doc 	= xml_open($REF_XML{$reference});

	my $DALI_SUMM_EXE = '/data/ecod/weekly_updates/weeks/bin/single_dali_summary.pl';
	if (!-f $DALI_SUMM_EXE) { die "ERROR! $DALI_SUMM_EXE not found\n"; } 

	my $cache;
	if ($REF_DALI_CACHE{$reference}) { 
		$cache = $REF_DALI_CACHE{$reference};
	}else{
		die "ERROR! $sub: No cache for $reference\n";
	}

	my %rep_stats;

	if (-f $cache) { 
		%rep_stats = %{retrieve($cache)};
		printf "Found %i reps in cache\n", scalar keys %rep_stats;
	}else{
		die "WARNING! $cache file not found\n";
	}
	my $rep_nodes = find_manual_domain_nodes($ecod_xml_doc);
	my @ecod_reps = find_uids($rep_nodes);
	build_rep_stats($rep_nodes, \%rep_stats);
	my $job_nodes = find_job_nodes($job_xml_doc);

	if ($DEBUG) { 
		printf "DEBUG: Found %i struct_search jobs\n", $job_nodes->size();
	}

	my ($job_dump_dir, $job_list_dir) = get_job_dirs_from_job_xml($job_xml_doc);

	my $reps = hasReps_job_list($job_xml_doc);
	my %query;
	my @job_ids;
	foreach my $job_node ($job_nodes->get_nodelist() ) { 

		my $job_id 	= $job_node->findvalue('@id');

		if ($reps && !$job_node->findvalue('@rep95') eq 'true') { next } 
		
		my ($pdb, $chain) = get_pdb_chain($job_node);
		my $pdb_chain	= $pdb . "_" . $chain;

		my $reference	= $job_node->findvalue('reference');
		my $mode 	= $job_node->findvalue('mode');

		if ($mode ne 'struct_search') { next } 

		my $recurse = 0;
		my $recurse_range;
		if ($job_node->findvalue('@type') eq 'recurse') { 
			my $seqid_range = $job_node->findvalue('seqid_range');
			$recurse_range = $seqid_range;
			$seqid_range =~ s/\,/_/g;
			$pdb_chain .= "_$seqid_range";
			$recurse++;
		}

		my $dali_summ_fn = "$job_dump_dir/$pdb_chain/$pdb_chain.$reference.dali_summ.xml";
		
		if (-f $dali_summ_fn && !$force_overwrite) { 
			print "WARNING! $dali_summ_fn already exists, skipping...\n";
			next;
		}

		my $dali_tmp_dir	= "$job_dump_dir/$pdb_chain/dali_tmp";

		my ($seqid_aref, $struct_seqid_aref, $pdbnum_href, $asym_id) = pdbml_seq_parse($pdb, $chain);

		if ($recurse) { 
			$seqid_aref = isect(range_expand($recurse_range), $seqid_aref);
			$struct_seqid_aref = isect(range_expand($recurse_range), $struct_seqid_aref);
		}
		#my $seq_href = pdbml_seq_fetch($pdb, $asym_id, $chain, $seqid_aref);

#		$query{pdb} 	= $pdb;
#		$query{chain} 	= $chain;
#		$query{asym_id} = $asym_id;
#		$query{struct_seqid_aref} = $struct_seqid_aref;
#		$query{seq}	= $seq_href;
#
		if (!-d $dali_tmp_dir) { 
			print "WARNING! Dali tmp dir $dali_tmp_dir not found\n";
			next;
		}

		my @dali_files;
		my %dali_files;
		if (0) { 
			my $missing_files = 0;
			foreach my $ecod_rep (@ecod_reps) { 

				my $dali_file = "$dali_tmp_dir/q.$pdb_chain.$ecod_rep.$reference.dali";
				if (-f $dali_file && !$dali_files{$dali_file}) { 
					push (@dali_files, $dali_file);
					$dali_files{$dali_file}++;
				}else{
					$missing_files++;
				}
			}
			printf "Found %i dali files\n", scalar(@dali_files);
			printf "%i dali files missing\n", $missing_files;
			if (scalar(@dali_files) == 0) { next } 
		}
		
		my @job_lns;
		if (! -f $dali_summ_fn || $force_overwrite) { 
			my $job_ln = "$DALI_SUMM_EXE $pdb $chain $job_dump_dir $reference";
			push (@job_lns, $job_ln);
			#my $job_fn = jobify($job_ln, "dali_summ.$pdb_chain.job");
			my $job_fn = "dali_summ.$pdb_chain.job";
			job_create($job_fn, \@job_lns);
			my $job_id = qsub($job_fn);
			push (@job_ids, $job_id);
		}else{
			print "WARNING! $dali_summ_fn found\n";
		}
	}

	return \@job_ids;
}

sub build_rep_stats { 
	my ($node_list, $rep_stats_href) = @_;

	foreach my $rep_node ($node_list->get_nodelist() ) { 

		my ($uid, $ecod_domain_id) = get_ids($rep_node);
		my $short_uid = substr($uid, 2, 5);

		next if isObsolete($rep_node);

		my $seqid_range = get_seqid_range($rep_node);

		my ($ecod_pdb, $ecod_chain) = get_pdb_chain($rep_node);

		if ($ecod_chain eq '.') { next}  #This is hacky and should be fixed 1/11/2013

		if ($$rep_stats_href{$uid}) { 
			#push (@ecod_reps, $ecod_domain_id);
			$$rep_stats_href{$uid}{range} = $seqid_range;
			next ;
		} else{
			print "$ecod_domain_id not in cache?\n";
		}

		my $pdb_path = "$DOMAIN_DATA_DIR/$short_uid/$uid/$uid.seqres.pdb";
		if (!-f "$DOMAIN_DATA_DIR/$short_uid/$uid/$uid.seqres.pdb") { 
			print "WARNING! No clean domain file for $ecod_domain_id\n";
		}else{
			#push (@ecod_reps, $ecod_domain_id);
			#$ecod_rep_range{$ecod_domain_id} = $seqid_range;
			$$rep_stats_href{$uid}{range} = $seqid_range;
		}

		my ($ecod_seqid_aref, $ecod_struct_seqid_aref, $ecod_pdbnum_href, $ecod_asym_id) = pdbml_seq_parse($ecod_pdb, $ecod_chain);
		if (!$ecod_seqid_aref) { 
			print "WARNING! $ecod_pdb, $ecod_chain obsolete in reps\n";
			next;
		}
		print "$ecod_pdb, $ecod_chain, $ecod_domain_id\n";
		my ($ecod_seq) = pdbml_seq_fetch($ecod_pdb, $ecod_asym_id, $ecod_chain, $ecod_seqid_aref);	
		if ($ecod_seqid_aref) { 
			$$rep_stats_href{$uid}{seqid_aref} = $ecod_seqid_aref;
			$$rep_stats_href{$uid}{struct_seqid_aref} = $ecod_struct_seqid_aref;
			$$rep_stats_href{$uid}{pdbnum_href} = $ecod_pdbnum_href;
			$$rep_stats_href{$uid}{asym_id} 	= $ecod_asym_id;
			$$rep_stats_href{$uid}{seq}	= $ecod_seq;
		}
	}
}
sub run_list_maintain { 
	my $sub = 'run_list_maintain';
	my ($run_list_xml_fn, $run_label, $run_dir, $reference, $start_date_ctime) = @_;;
	my $run_list_xml_doc;
	print "$sub: $run_list_xml_fn\n";
	if (-f $run_list_xml_fn)  { 
		$run_list_xml_doc = xml_open($run_list_xml_fn);
	}else{
		my ($run_list_xml_doc, $run_list_doc_node) = xml_create('domain_parse_run_list');
		$run_list_xml_doc->setDocumentElement($run_list_doc_node);
		$run_list_doc_node->setAttribute('tree_struct', 'pdb_chain');
	}

	my $run_list_dp_XPath = '//domain_parse_run';
	my @run_list_week_dirs;
	my $MAX_RUN_ID = 0;
	foreach my $dp_node (find_domain_parse_nodes($run_list_xml_doc)) { 
		my $week_label 	 = $dp_node->findvalue('week_label');
		my $run_list_dir = $dp_node->findvalue('run_list_dir');
		my $run_id	 = $dp_node->findvalue('@run_id');
		push (@run_list_week_dirs, $week_label);
		if ($run_id > $MAX_RUN_ID) { $MAX_RUN_ID = $run_id }
	}

	opendir(DIR, $status_dir) or die "ERROR! Could not open status directory for reading:$!\n";
	my @week_dirs = grep {$_ =~ /\d{7}/} readdir(DIR);
	close DIR;

	my $obsolete_file = '/usr2/pdb/data/status/obsolete.dat';
	if (!-f $obsolete_file) { 
		die "ERROR! Could not find obsolete list $obsolete_file\n";
	}

	my %obsolete;
	open (my $fh, "<", $obsolete_file) or die "ERROR! Could not open obsolete manifest file $obsolete_file: $!\n";
	while (my $ln = <$fh>) { 
		if (substr($ln, 0, 6) ne 'OBSLTE') { next } 
		my $obs = substr($ln, 20, 4);
		$obsolete{$obs}++;
	}
	close $fh;

	my $run_list_head_node = find_first_node($run_list_xml_doc, "domain_parse_run_list");

	my @new_week_dirs;
	foreach my $week (sort {$a <=> $b }  @week_dirs) { 

		my $week_time = str2time($week);

		if ($week_time <= $start_date_ctime) { next } 
		if (grep {$_ eq $week} @run_list_week_dirs) { next } 
		
		push (@new_week_dirs, $week);

	}

	printf "Found %i new PDB weeks\n", scalar(@new_week_dirs);

#Generate run_list nodes, directories, and job files;

	my $new_run_id = $MAX_RUN_ID+1;
	foreach my $new_week (@new_week_dirs) { 
		
		print "n: $new_week\n";

		if (! -d "$run_dir/$new_week") { 
			print "WARNING! Creating week dir $run_dir/$new_week\n";
			if (! mkdir("$run_dir/$new_week")) { 
				die "ERROR! Could not create $run_dir/$new_week\n";
			}else{
				chown($UID, $GID, "$run_dir/$new_week");
			}
		}

		my @new_pdbs;
		my $warning = 0;
		my $fh;
		if (!open($fh, "<",  "$status_dir/$new_week/added.pdb")) { 
			die "ERROR! Could not open $status_dir/$new_week/added.pdb for reading:$!\n";

		}else{
			while (my $ln = <$fh>) { 
				my $pdb = $ln;
				chomp $pdb;
				my $two = substr($pdb, 1, 2);
				if ($obsolete{$pdb}) { 
					next;
				#}elsif (-f "$pdb_dir/$two/pdb$pdb.ent.gz") { 
				}elsif (-f "$pdbml_dir/$two/$pdb-noatom.xml.gz") { 
					push (@new_pdbs, $pdb);
				}else{
					print "WARNING! $pdb not found in $pdb_dir/$two\n";
					$warning++;
				}
			}
		}
		if ($warning > 10) { 
			print "WARNING! week $new_week not complete (missing_pdb = $warning), skipping\n";
			next;
		}

		my $domain_parse_run_node	= $run_list_xml_doc->createElement('domain_parse_run');
		$domain_parse_run_node->setAttribute('run_id', $new_run_id++);
		$new_run_id++;

		$run_list_head_node->appendChild($domain_parse_run_node);

		my $week_label_node		= $run_list_xml_doc->createElement('week_label');
		$week_label_node->appendTextNode($new_week);
		$domain_parse_run_node->appendChild($week_label_node);

		my $run_list_dir_node	= $run_list_xml_doc->createElement('run_list_dir');
		$run_list_dir_node->appendTextNode("$run_dir/$new_week");
		$domain_parse_run_node->appendChild($run_list_dir_node);

		#generate run list;
		my $job_list_xml_fn = "$run_dir/$new_week/$new_week.$run_label.job.xml";
		my $job_xml_doc;
		my ($job_dump_dir, $job_list_dir);
		if (-f $job_list_xml_fn) { 
			print "WARNING! $job_list_xml_fn exists, set \$FORCE_OVERWRITE to overwrite\n";
			$job_xml_doc = xml_open($job_list_xml_fn);
			($job_dump_dir, $job_list_dir) = get_job_dirs_from_job_xml($job_xml_doc);
		}elsif (!-f $job_list_xml_fn || $FORCE_OVERWRITE) {

			my $domain_parse_job_xml_node 	= $run_list_xml_doc->createElement('run_list_job_xml_file');
			$domain_parse_job_xml_node->appendTextNode("$new_week.$run_label.job.xml");
			$domain_parse_job_xml_node->setAttribute('mode', 'seq_iter');
			$domain_parse_run_node->appendChild($domain_parse_job_xml_node);
			
			$job_xml_doc	= XML::LibXML->createDocument();
			my $run_list_document_node	= $job_xml_doc->createElement("job_set_top");
			$job_xml_doc->setDocumentElement($run_list_document_node);

			$job_list_dir = "$run_dir/$new_week";
			my $job_list_dir_node = $job_xml_doc->createElement('job_list_dir');
			$job_list_dir_node->appendTextNode($job_list_dir);
			$run_list_document_node->appendChild($job_list_dir_node);

			$job_dump_dir = "$run_dir/$new_week/ecod_dump";
			if (!-d $job_dump_dir) { 
				if (!mkdir($job_dump_dir)) { 
					die "ERROR! Could not create $job_dump_dir\n";
				}else{
					chown($UID, $GID, $job_dump_dir);
				}
			}
			my $job_dump_dir_node	= $job_xml_doc->createElement('job_dump_dir');
			$job_dump_dir_node->appendTextNode($job_dump_dir);
			$run_list_document_node->appendChild($job_dump_dir_node);

			my $job_list_node	= $job_xml_doc->createElement('job_list');
			$run_list_document_node->appendChild($job_list_node);
			
			my $job_id = 0;
			for (my $i = 0; $i < scalar(@new_pdbs); $i++) { 
				my $pdb = $new_pdbs[$i];
				my $chain_aref = pdbml_fetch_chains($pdb);
				
				for (my $j = 0; $j < scalar(@$chain_aref); $j++) { 
					my $chain = $$chain_aref[$j];
					#print "$pdb, $chain\n";

					my $query_pdb_node	= $job_xml_doc->createElement('query_pdb');
					$query_pdb_node->appendTextNode($pdb);
					
					my $query_chain_node	= $job_xml_doc->createElement('query_chain');
					$query_chain_node->appendTextNode($chain);

					my $reference_node	= $job_xml_doc->createElement('reference');
					$reference_node->appendTextNode($reference);

					my $mode_node		= $job_xml_doc->createElement('mode');
					$mode_node->appendTextNode("seq_iter");
					$mode_node->setAttribute('input', "struct_seqid");

					my $job_node	= $job_xml_doc->createElement('job');
					$job_node->setAttribute('id', $job_id);

					$job_node->appendChild($query_pdb_node);
					$job_node->appendChild($query_chain_node);
					$job_node->appendChild($reference_node);
					$job_node->appendChild($mode_node);

					$job_list_node->appendChild($job_node);

					$job_id++;
				}
			}
			if ($DEBUG) { 
				print "DEBUG: doc string output $new_week\n";
			}

			xml_write($job_xml_doc, $job_list_xml_fn);
		}

		#job_list_maintain($job_list_xml_fn, $new_week, $reference);
		#Generate FASTA files.
		foreach my $job_node (find_job_nodes($job_xml_doc)) { 

			my ($query_pdb, $query_chain) = get_pdb_chain($job_node);
			my $pdb_chain	= $query_pdb . "_" . $query_chain;

			if (!-d "$job_dump_dir/$pdb_chain") { 
				if (!mkdir("$job_dump_dir/$pdb_chain")) { 
					die "ERROR! Could not create $job_dump_dir/$pdb_chain\n";
				}else{
					chown($UID, $GID, "$job_dump_dir/$pdb_chain");
				}
			}
			my $fasta_fn = "$job_dump_dir/$pdb_chain/$pdb_chain.fa";
			if (-f $fasta_fn && $FORCE_OVERWRITE < 2) { 
				print "WARNING! $fasta_fn exists, skipping...\n";
				next;
			}else{
				my ($seqid_aref, $struct_seqid_aref, $pdbnum_aref, $asym_id) = pdbml_seq_parse($query_pdb, $query_chain);
				if (!$struct_seqid_aref || $struct_seqid_aref == 0 || scalar(@$struct_seqid_aref) == 0) {
					print "WARNING! No struct_seq for $query_pdb, $query_chain\n";
					$job_node->setAttribute('unstructured', 'true');
					next;
				}

				#my $fasta_string = pdbml_fasta_fetch($query_pdb, $asym_id, $query_chain, $seqid_aref);

				my $struct_seqid_range 		= rangify(@$struct_seqid_aref);
				my $ungapped_struct_range 	= ungap_range($struct_seqid_range, $GAP_TOL);
				my $ungap_aref 			= range_expand($ungapped_struct_range);

				if ($ungapped_struct_range eq 0) { 
					print "WARNING! empty range for $query_pdb, $query_chain\n";
					next;
				}
				#my $fasta_string = pdbml_fasta_fetch($query_pdb, $asym_id, $query_chain, $struct_seqid_aref);
				my $fasta_string = pdbml_fasta_fetch($query_pdb, $asym_id, $query_chain, $ungap_aref);
				if (!$fasta_string) { 
					print "WARNING! No fasta string for $query_pdb, $query_chain\n";
					next;
				}

				open (my $fh,">", "$fasta_fn") or die "ERROR! Could not open $fasta_fn for writing:$!\n";
				print $fh ">$query_pdb,$query_chain\n$fasta_string\n";
				close $fh;
			}
		}

		#Generate run list cluster
		if ($DEBUG) { 
			print "DEBUG: cluster: $new_week\n";
		}
		my $cluster_script = '/data/ecod/weekly_updates/weeks/bin/run_list_cluster.pl';
		if (!-f $cluster_script) { 
			die "run list cluster script not found! $cluster_script\n";
		}
		open ($fh, ">", "$run_label.$new_week.cluster.job") or die "Could not open $run_label.$new_week.cluster.job for writing:%!\n";
		print $fh "#!/bin/bash\n";
		print $fh "#\$ -cwd\n";
		print $fh "#\$ -j y \n";
		print $fh "#\$ -S /bin/bash\n";
		print $fh "#\$ -M dustin.schaeffer\@gmail.com\n";
		print $fh "$cluster_script $job_list_xml_fn\n";
		close $fh;

		my $clust_job_id = qsub("$run_label.$new_week.cluster.job");

		print "$clust_job_id\n";

		while (qstat_wait($clust_job_id)) { 
			print "SLEEPING cluster\n";
			sleep(10);
		}

		#Build profiles
		my $buildali_script = '/data/ecod/weekly_updates/weeks/bin/xml_hhblits_profile_build.pl';
		if (!-f $buildali_script) { 
			die "buildali script not found! $buildali_script\n";
		}

		if ($DEBUG) { 
			print "DEBUG: buildali_hhblits: $new_week\n";
		}
		my $ali_job_ids = buildali_hhblits($job_list_xml_fn);

		while (qstat_wait_list($ali_job_ids)) { 
			print "SLEEPING...\n";
			sleep(300); #5 min job checks;
		}

		#Generate peptide filter
		if ($DEBUG) { 
			print "DEBUG: peptide_prefilter $new_week\n";
		}
		#my $peptide_filter_script = '/home/rschaeff/work/domain_partition/ecod_domain_parser/bin/xml_peptide_prefilter.pl';
#		my $peptide_filter_script = '/data/ecod/weekly_updates/weeks/bin/xml_peptide_prefilter.pl';
#		if (!-f $peptide_filter_script) { 
#			die "peptide filter script not found! $peptide_filter_script\n";
#		}
#		open (OUT, ">$run_label.$new_week.pep.job") or die "Could not open $run_label.$new_week.pep.job for writing:$!\n";
#		print OUT "#!/bin/bash\n";
#		print OUT "#\$ -cwd\n";
#		print OUT "#\$ -j y \n";
#		print OUT "#\$ -S /bin/bash\n";
#		print OUT "#\$ -M dustin.schaeffer\@gmail.com\n";
#		print OUT "$peptide_filter_script $job_list_xml_fn\n";
#		close OUT;
#
		my $pep_job_ids = peptide_filter($job_list_xml_fn);

		#my $pep_job_id = qsub("$run_label.$new_week.pep.job");
		while (qstat_wait_list($pep_job_ids)) { 
			print "SLEEPING peptide!\n";
			sleep(10);
		}
		peptide_collate($job_list_xml_fn);

		#Run self_comp jobs
		if ($DEBUG) { 
			print "DEBUG: self_comp: $new_week\n";
		}
		#my $self_comp_script = '/home/rschaeff/work/domain_partition/ecod_domain_parser/bin/ecod_self_comparison.pl';
		my $self_comp_script = '/data/ecod/weekly_updates/weeks/bin/ecod_self_comparison.pl';
		if (!-f $self_comp_script) { 
			die "self comp script not found! $self_comp_script\n";
		}

		open (OUT, ">$run_label.$new_week.self_comp.job") or die "Could not open $run_label.$new_week.self_comp.job for writing:$!\n";
		print OUT "#!/bin/bash\n";
		print OUT "#\$ -cwd\n";
		print OUT "#\$ -j y \n";
		print OUT "#\$ -S /bin/bash\n";
		print OUT "#\$ -M dustin.schaeffer\@gmail.com\n";
		print OUT "$self_comp_script $job_list_xml_fn\n";
		close OUT;

		my $self_comp_job_id = qsub("$run_label.$new_week.self_comp.job");

		#Run Chain-blast jobs
		if ($DEBUG ) { 
			print "DEBUG: chblastp: $new_week\n";
		}
		my $reps_only = 0;
		my $force_overwrite = 0;
		my $chblast_job_ids = chblastp_job_file($job_list_xml_fn, $reps_only, $force_overwrite);

		while (qstat_wait_list($chblast_job_ids)) { 
			print "SLEEPING...\n";
			sleep(30); #5 min job checks;
		}
		
		#Run blastp jobs 
		if ($DEBUG) { 
			print "DEBUG blastp: $new_week\n";
		}
		my $blast_job_ids = blastp_job_file($job_list_xml_fn, $reps_only);

		while(qstat_wait_list($blast_job_ids)) { 
			print "SLEEPING...\n";
			sleep(30);
		}

		#Run hhsearch jobs dependent on ali
		#Split hh search into run and parse components and allow for parallelization

		if ($DEBUG) { 
			print "DEBUG hh: $new_week\n";
		}

		my $hh_run_job_ids = hh_run($job_list_xml_fn, $reference, 1);

		while (qstat_wait_list($hh_run_job_ids)) { 
			print "SLEEPING...\n";
			sleep(30);
		}

		if ($DEBUG) { 
			print "DEBUG hh_parse: $new_week\n";
		}

		my $hh_parse_job_ids = hh_parse($job_list_xml_fn, $reference, 1);

		while (qstat_wait_list($hh_parse_job_ids)) { 
			print "SLEEPING...\n";
			sleep(30);
		}
		
		#Run blast domain summary
		if ($DEBUG) { 
			print "DEBUG bsumm: $new_week\n";
		}

		my $bsumm_job_ids = bsumm($job_list_xml_fn, 1);

		while (qstat_wait_list($bsumm_job_ids)) { 
			print "SLEEPING...\n";
			sleep(30);
		}

		#Run seq_iter domain partition

		if ($DEBUG) { 
			print "DEBUG domains: $new_week\n";
		}

		my $domains_job_ids = domains($job_list_xml_fn, $reference, 1, $FORCE_OVERWRITE);
		while (qstat_wait_list($domains_job_ids)) { 
			print "SLEEPING...\n";
			sleep(30);
		}
#
		#Run boundary optimize
		my $optimize_job_ids = optimize($job_list_xml_fn);
		while (qstat_wait_list($optimize_job_ids)) { 
			print "SLEEPING...\n";
			sleep(30);
		}
	}

	generate_struct_search_jobs($run_list_xml_doc);
	xml_write($run_list_xml_doc, $run_list_xml_fn);

	my $DALI_SEARCH = 1;
	if ($DALI_SEARCH) { 
		foreach my $run_node ($run_list_xml_doc->findnodes('//domain_parse_run')) { 
			if (existsMode($run_node, 'struct_search')) { 
				my $run_list_dir 	= $run_node->findvalue('run_list_dir');
				my $job_list_xml_fn 	= $run_node->findvalue('run_list_job_xml_file[@mode="struct_search"]');

				my $job_list_xml_path = "$run_list_dir/$job_list_xml_fn";
				struct_search_dali_query_gen($job_list_xml_path);
				my $dali_jobs_aref = query_glob($job_list_xml_path);
				throttled_qsub($dali_jobs_aref, 20000);
				my $dali_summ_job_ids = dali_summary($job_list_xml_path, $LATEST_REFERENCE);
				while (qstat_wait_list($dali_summ_job_ids)) { 
					print "SLEEPING...\n";
					sleep(30);
				}
			}
		}
	}
}

sub query_glob { 
	my $sub = 'query_glob';

	my ($job_list_xml_fn) = @_;

	my $job_xml_doc = xml_open($job_list_xml_fn);

	my ($job_dump_dir, $job_list_dir) = get_job_dirs_from_job_xml($job_xml_doc);

	my @glob_files;
	foreach my $job_node ($job_xml_doc->findnodes('//job')) { 
		my ($query_pdb, $query_chain) = get_pdb_chain($job_node);
		my $pc = $query_pdb . "_" . $query_chain;

		my @job_files 	= glob("$job_dump_dir/$pc/dali_tmp/*.job");
		my @dali_files 	= glob("$job_dump_dir/$pc/dali_tmp/*.dali");
		printf "#%i %i\n", scalar(@job_files), scalar(@dali_files);

		my %dali_lookup;
		foreach my $dali_file (@dali_files) { 
			$dali_file =~ s/dali$/job/;
			$dali_lookup{$dali_file}++;	
		}

		foreach my $job_file (@job_files)  { 
			if (!$dali_lookup{$job_file}) { 
				push (@glob_files, $job_file);
			}
		}

	}
	
	printf "$sub: Found %i files\n", scalar(@glob_files);

	return \@glob_files;
}
sub job_list_generate_fasta { 
	my $sub = 'job_list_generate_fasta';
	my ($job_xml_doc, $force_overwrite) = @_;
	#Generate FASTA files.
	my ($job_dump_dir, $job_list_dir) = get_job_dirs_from_job_xml($job_xml_doc);

	foreach my $job_node (find_job_nodes($job_xml_doc)) { 

		my ($query_pdb, $query_chain) = get_pdb_chain($job_node);
		my $pdb_chain	= $query_pdb . "_" . $query_chain;

		if (!-d "$job_dump_dir/$pdb_chain") { 
			if (!mkdir("$job_dump_dir/$pdb_chain")) { 
				die "ERROR! Could not create $job_dump_dir/$pdb_chain\n";
			}else{
				chown($UID, $GID, "$job_dump_dir/$pdb_chain");
			}
		}
		my $fasta_fn = "$job_dump_dir/$pdb_chain/$pdb_chain.fa";
		if (-f $fasta_fn && $force_overwrite) { 
			print "WARNING! $fasta_fn exists, skipping...\n";
			next;
		}else{
			my ($seqid_aref, $struct_seqid_aref, $pdbnum_aref, $asym_id) = pdbml_seq_parse($query_pdb, $query_chain);
			if (!$struct_seqid_aref || $struct_seqid_aref == 0 || scalar(@$struct_seqid_aref) == 0) {
				print "WARNING! No struct_seq for $query_pdb, $query_chain\n";
				$job_node->setAttribute('unstructured', 'true');
				next;
			}

			#my $fasta_string = pdbml_fasta_fetch($query_pdb, $asym_id, $query_chain, $seqid_aref);

			my $struct_seqid_range 		= rangify(@$struct_seqid_aref);
			my $ungapped_struct_range 	= ungap_range($struct_seqid_range, $GAP_TOL);
			my $ungap_aref 			= range_expand($ungapped_struct_range);

			#New query strings are generated on ungapped stuct_seqids, so that small gaps filled in by domain partition, are also
			#filled in the query
			my $fasta_string = pdbml_fasta_fetch($query_pdb, $asym_id, $query_chain, $ungap_aref);
			if (!$fasta_string) { 
				print "WARNING! No fasta string for $query_pdb, $query_chain\n";
				next;
			}

			open (OUT, ">$fasta_fn") or die "ERROR! Could not open $fasta_fn for writing:$!\n";
			print OUT ">$query_pdb,$query_chain\n$fasta_string\n";
			close OUT;
		}
	}
}

sub job_list_cluster_fasta { 
	my ($job_list_xml_fn) = @_;
	my $sub = "job_list_cluster_fasta";
	if (!-f $job_list_xml_fn) { die "ERROR! $sub: File not found $job_list_xml_fn:$!\n"; }

	my $cluster_script = '/data/ecod/weekly_updates/weeks/bin/run_list_cluster.pl';
	if (!-f $cluster_script) { 
		die "run list cluster script not found! $cluster_script\n";
	}

	my $label;
	if ($job_list_xml_fn =~ /(.*)\.job/) { 
		$label = $1;
	}else{
		my @chars = ("A" .. "Z", "a" .. "z");
		$label .= $chars[rand @chars] for 1..8;
	}
	my $job_fn = "$label.cluster.job";
	my @job_lns;
	my $cmd = "$cluster_script $job_list_xml_fn\n";
	push (@job_lns, $cmd);
	job_create($job_fn, \@job_lns);
	my $clust_job_id = qsub($job_fn);

	print "$clust_job_id\n";

	while (qstat_wait($clust_job_id)) { 
		print "SLEEPING cluster\n";
		sleep(10);
	}
}



#sub domain_partition_asm { 
#	my $sub = 'domain_partition_asm';
#
#	my ($query_pdb, $query_chains, $dir, $input_mode, $reference, $domain_prefix, $domain_xml_fn, $ref_range_cache, $concise, $global_opt, $ref_chain_domains, $ref_domain_uid_lookup) = @_;
#
#
#	#Build chains strings
#	my $query_chains_str = join("", @$query_chains);
#	my $pdb_chains = $query_pdb . "_" . $query_chains_str; 
#
#	#Skip domain partition if file exists and force overwrite is not specificied
#	if (-f $domain_xml_fn && !$$global_opt{force_overwrite}) { 
#		print "WARNING! $sub: skipping $query_pdb, $query_chains_str, set FORCE_OVERWRITE bit to override...\n"; return; 
#	} 
#	my ($query_seqid_aref, 
#		$query_struct_seqid_aref,
#		$query_pdbnum_aref,
#		$query_asym_id,
#		$query_chain_aref,
#		$query_struct_chain_aref) = pdbml_mc_seq_parse($query_pdb, $query_chains); 
#
#	my %super_seqid_sort;
#	my %super_seqid_pdbnum_lookup;
#	my %super_seqid_chain_lookup;
#	for (my $i = 0; $i < scalar(@$query_seqid_aref); $i++) { 
#		$super_seqid_sort{$$query_chain_aref[$i]}{$$query_seqid_aref[$i]} = $i;
#		$super_seqid_pdbnum_lookup{$i}	= "$$query_chain_aref[$i]$$query_seqid_aref[$i]";
#		$super_seqid_chain_lookup{$i}	= $$query_chain_aref[$i];
#	}
#
#	my %super_struct_seqid_sort;
#	my %super_struct_seqid_pdbnum_lookup;
#	my %super_struct_seqid_chain_lookup;
#	for (my $i = 0; $i < scalar(@$query_struct_seqid_aref); $i++) { 
#		$super_struct_seqid_sort{$$query_struct_chain_aref[$i]}{$$query_struct_seqid_aref[$i]} = $i;
#
#		$super_struct_seqid_pdbnum_lookup{$i}	= "$$query_struct_chain_aref[$i]$$query_struct_seqid_aref[$i]";
#		$super_struct_seqid_chain_lookup{$i}	= $$query_struct_chain_aref[$i];
#	}
#
#	if (!$query_seqid_aref) { 
#		print "WARNING! $sub: $query_pdb not found! skipping...\n";
#		return 0;
#	}
#
#	#if ($recurse) { 
#	#	die "ERROR! recurse not suppported for assembly jobs\n";
#	#}
#
#	my $query_range_str		= multi_chain_rangify($query_seqid_aref, $query_chain_aref);
#	my $query_struct_range_str	= multi_chain_rangify($query_struct_seqid_aref, $query_struct_chain_aref);
#
#	#Changed struct_seqid to be ungapped (i.e. may contain small unstructured regions);
#	print "?$query_struct_range_str $$global_opt{gap_tol}\n";
#	my $query_ungapped_struct_range_str 			= multi_chain_ungap_range($query_struct_range_str, $$global_opt{gap_tol});
#	($query_struct_seqid_aref, $query_struct_chain_aref) 	= multi_chain_range_expand($query_ungapped_struct_range_str);
#
#	if (!$query_struct_range_str) { 
#		print "WARNING! $sub: $query_pdb $query_chains_str.. No structured range $query_struct_range_str\n";
#	}
#
#	#Define your unused set
#	my @unused_seq;
#	my @unused_chain;
#	my %sort_lookup;
#	my %pdbnum_lookup;
#	my %chain_lookup;
#	if ($input_mode	eq 'struct_seqid') { 
#		@unused_seq	= @$query_struct_seqid_aref;
#		@unused_chain	= @$query_struct_chain_aref;
#		%sort_lookup	= %super_struct_seqid_sort;
#		%pdbnum_lookup 	= %super_seqid_pdbnum_lookup;
#		%chain_lookup 	= %super_seqid_chain_lookup;
#	}elsif($input_mode eq 'seqid') { 
#		@unused_seq	= @$query_seqid_aref;
#		@unused_chain 	= @$query_chain_aref;
#		%sort_lookup 	= %super_seqid_sort;
#		%pdbnum_lookup 	= %super_struct_seqid_pdbnum_lookup;
#		%chain_lookup 	= %super_struct_seqid_chain_lookup;
#	}else{
#		die "ERROR! $sub: Unknown input mode $input_mode\n";
#	}
#
#	#Used arrays
#	my @used_seq;
#	my @used_chain;
#
#	if ($DEBUG) { 
#		print "DEBUG $sub: $query_pdb $query_chains_str $dir $query_range_str $query_struct_range_str $domain_xml_fn\n";
#	}
#
#	#Define the path to the BLAST summary file, check if mode is concise
#	my $blast_summ_fn;
#	if ($concise) {
#		$blast_summ_fn	= "$dir/$pdb_chains.$reference.blast_summ.concise.xml";
#	}else{
#		$blast_summ_fn 	= "$dir/$pdb_chains.$reference.blast_summ.xml";
#	}
#
#
#	#If the file isn't there, check for concise anyways, then skip
#	if (!-f $blast_summ_fn) { 
#		if (-f "$dir/$pdb_chains.$reference.blast_summ.concise.xml") { 
#			$blast_summ_fn = "$dir/$pdb_chains.$reference.blast_summ.concise.xml";
#		}else{
#			print "WARNING! $sub: BLAST summary for $query_pdb, $query_chains_str, $dir ($blast_summ_fn) not found\n";
#			next;
#		}
#	}
#	print "bs: $blast_summ_fn\n";
#
#	#Set up the domain XML document
#	my $domain_count	= 1;
#
#	my $domain_xml_doc	= XML::LibXML->createDocument;
#	my $domain_doc_node	= $domain_xml_doc->createElement('chain_domains_set_top');
#	$domain_xml_doc->setDocumentElement($domain_doc_node);
#
#	$domain_doc_node->setAttribute('pdb_id', $query_pdb);
#	$domain_doc_node->setAttribute('chain_str', $query_chains_str);
#	$domain_doc_node->setAttribute('mode','blast');
#	$domain_doc_node->setAttribute('domain_prefix', $domain_prefix);
#
#	my $lc_query_pdb = lc($query_pdb);
#
#	my $domain_list_node	= $domain_xml_doc->createElement('domain_list');
#	$domain_doc_node->appendChild($domain_list_node);
#
#	my $blast_summ_xml_doc = xml_open($blast_summ_fn);
#
#	my $blast_summ_pdb 	= $blast_summ_xml_doc->findvalue('//blast_summ_doc/blast_summ/@pdb');
#	my $blast_summ_chains	= $blast_summ_xml_doc->findvalue('//blast_summ_doc/blast_summ/@chains');
#
#	if ($query_pdb ne $blast_summ_pdb || $query_chains_str ne $blast_summ_chains) { 
#		print "ERROR! Blast file pdb mismatch $query_pdb/$blast_summ_pdb $query_chains_str/$blast_summ_chains\n";
#		die;
#	}
#	#Self_comparison read:
#	my $self_comp_run_XPath = q{//blast_summ_doc/blast_summ/self_comp_run};
#	my $self_comp_run_nodes = $blast_summ_xml_doc->findnodes($self_comp_run_XPath);
#	my $self_comp_run_nodes_size = $self_comp_run_nodes->size();
#	my @self_comps;
#	if ($self_comp_run_nodes_size < 1) { 
#		print "WARNING! $sub: No self_comp nodes found\n";
#	}else{
#		my $self_comp_hit_nodes = $self_comp_run_nodes->get_node(1)->findnodes('hits/hit');
#		my $ci = 0;
#		foreach my $hit_node ($self_comp_hit_nodes->get_nodelist() ) { 
#
#			my $aligner 	= $hit_node->findvalue('@aligner');
#			my $zscore	= $hit_node->findvalue('@z_score');
#			#my $prob	= $hit_node->findvalue('@prob');
#
#			my $query_reg	= $hit_node->findvalue('query_reg');
#			my $hit_reg	= $hit_node->findvalue('hit_reg');
#			
#			if (! $query_reg =~ /\d+/) { 
#				next;
#			}
#
#			if (!$hit_reg =~ /\d+/) { 
#				next;
#			}
#
#			$self_comps[$ci]{aligner}	= $aligner;
#			if ($aligner eq 'dali') { 
#				$self_comps[$ci]{zscore} = $zscore;
#			}elsif($aligner eq 'hhrepid') { 
#				#$self_comps[$ci]{prob}	= $prob;
#			}else{
#				print "WARNING! Unknown self_comparison aligner $aligner, skipping...\n";
#				next;
#			}
#
#			$self_comps[$ci]{query_reg}	= $query_reg;
#			$self_comps[$ci]{hit_reg}	= $hit_reg;
#			$ci++;
#		}
#	}
#
#	#Chainwise BLAST
#	my $chblastp_run_XPath 	= q{//blast_summ_doc/blast_summ/chain_blast_run[@program='blastp']};
#	my $chblastp_run_nodes	= $blast_summ_xml_doc->findnodes($chblastp_run_XPath);
#	my $chblastp_run_nodes_size	= $chblastp_run_nodes->size();
#	if ($chblastp_run_nodes_size != 1) { 
#		print "ERROR! $sub: Odd number of chblastp run nodes found ($chblastp_run_nodes_size), aborting...\n";
#	}
#	my $chblastp_hit_nodes	= $chblastp_run_nodes->get_node(1)->findnodes('hits/hit');
#	foreach my $hit_node ($chblastp_hit_nodes->get_nodelist() ) { 
#		my $hit_num 	= $hit_node->findvalue('@num');
#		my $hit_pdb	= $hit_node->findvalue('@pdb_id');
#		my $hit_chain	= $hit_node->findvalue('@chain_id');
#		my $hsp_count 	= $hit_node->findvalue('@hsp_count');
#		my $evalues 	= $hit_node->findvalue('@evalues');
#
#		#remove when possible
#		if ($hit_chain eq '.') { 
#			next;
#		}
#
#		my $evalue = 'Unk';
#		if ($hsp_count == 1) { 
#			$evalue = $evalues;
#		}
#
#		my $query_reg		= $hit_node->findvalue('query_reg');
#		#ARG 1 is pos-based, ARG2 is seqid_based, OUTPUT is segid
#		#Construct range and chain arrays from the query region
#		my $hit_query_struct_seqid_aref;
#		my $hit_query_struct_chain_aref;
#		if ($input_mode eq 'seqid') { 
#			#$hit_query_struct_seqid_aref = range_expand($query_reg);
#			($hit_query_struct_seqid_aref, $hit_query_struct_chain_aref) = multi_chain_struct_region(range_expand($query_reg), $query_seqid_aref, $query_chain_aref)
#		}elsif($input_mode eq 'struct_seqid') { 
#			($hit_query_struct_seqid_aref, $hit_query_struct_chain_aref) = multi_chain_struct_region(range_expand($query_reg), $query_struct_seqid_aref, $query_struct_chain_aref);
#		}else{
#			die "ERROR! $sub: input_mode $input_mode not recognized\n";
#		}
#		my $query_struct_seqid_range = multi_chain_rangify($hit_query_struct_seqid_aref, $hit_query_struct_chain_aref);
#
#		my $hit_reg		= $hit_node->findvalue('hit_reg');
#		my $query_seq		= $hit_node->findvalue('query_seq');
#		my $hit_seq		= $hit_node->findvalue('hit_seq');
#
#		if ($DEBUG) { 
#			print "DEBUG $sub: chblastp $hit_num $hit_pdb $hit_chain $query_struct_seqid_range\n";
#		}
#		#Assign if overlaps 
#
#		my $unused_coverage 	= multi_chain_region_coverage($query_struct_seqid_aref, $query_struct_chain_aref, \@unused_seq, \@unused_chain);
#		my $used_coverage 	= multi_chain_region_coverage($query_struct_seqid_aref, $query_struct_chain_aref, \@used_seq, \@used_chain);
#
#		if ($unused_coverage < 0.05 || scalar(@unused_seq) < 10) { 
#			last;
#		}
#
#		my $query_new_coverage 	= multi_chain_region_coverage($hit_query_struct_seqid_aref, $hit_query_struct_chain_aref, \@unused_seq, \@unused_chain);
#		my $query_used_coverage = multi_chain_region_coverage($hit_query_struct_seqid_aref, $hit_query_struct_chain_aref, \@used_seq, \@used_chain);
#
#		if ($DEBUG) { 
#			print "DEBUG $sub: uc $used_coverage uuc $unused_coverage qnc $query_new_coverage qoc $query_used_coverage\n";
#		}
#
#
#		if ($query_new_coverage > $$global_opt{new_coverage_threshold} && $query_used_coverage < $$global_opt{old_coverage_threshold}) { 
#
#			my ($hit_seqid_aref, $hit_struct_seqid_aref, $hit_pdbnum_aref, $hit_asym_id, $hit_chain_aref, $hit_struct_chain_aref);
#				($hit_seqid_aref, $hit_struct_seqid_aref, $hit_pdbnum_aref, $hit_asym_id) = pdbml_seq_parse($hit_pdb, $hit_chain);
#
#			if (!$hit_seqid_aref) { next } 
#			print "DEFINE chblastp: $sub: ($hit_pdb $hit_chain $hit_reg) ->  ($query_struct_seqid_range)\n";
#
#			#Calculate hit->query seqid map;
#			my @hit_segs 	= split(",", $hit_reg);
#			my @query_segs 	= split(",", $query_reg);
#
#			my @hit_seqs 	= split(",", $hit_seq);
#			my @query_seqs 	= split(",", $query_seq);
#
#			my %query_map;
#			my %hit_map;
#			my %hit_chain_map;
#			for (my $i = 0; $i < scalar(@hit_seqs); $i++) { 
#				my $hit_seg	= $hit_segs[$i];
#				#Following line is a bug if library is not structured residues.
#				if ($DEBUG) { 
#					print "DEBUG $sub: chblastp struct_range_expand 1 $hit_seg\n";
#				}
#
#				my ($hit_seg_seqid_aref, $hit_seg_chain_aref);
#				$hit_seg_seqid_aref	= struct_region(range_expand($hit_seg), $hit_struct_seqid_aref);
#
#				my @hit_seq	= split("", $hit_seqs[$i]);
#				my $query_seg	= $query_segs[$i];
#				
#				if ($DEBUG) { 
#					printf "DEBUG $sub: chblastp struct_range_expand 2 $query_seg %s\n", rangify(@$hit_query_struct_seqid_aref);
#				}
#
#				my $query_seg_seqid_aref;
#				my $query_seg_chain_aref;
#				if ($input_mode eq 'seqid') { 
#					#$query_seg_seqid_aref	= range_expand($query_seg);
#					$query_seg_seqid_aref	= isect(range_expand($query_seg), $query_struct_seqid_aref);
#				}elsif($input_mode eq 'struct_seqid') { 
#					$query_seg_seqid_aref	= struct_region(range_expand($query_seg), $query_struct_seqid_aref);
#					$query_seg_chain_aref 	= struct_region(range_expand($query_seg), $query_struct_chain_aref);
#				}
#				if ($DEBUG) { 
#					printf "DEBUG $sub: chblastp query_seg_seqid %s\n", rangify(@$query_seg_seqid_aref);
#				}
#
#				my @query_seq	= split("", $query_seqs[$i]);
#
#				my $hit_start = $$hit_seg_seqid_aref[0];
#				my $query_start	= $$query_seg_seqid_aref[0];
#				
#				my $hit_pos 	= 0;
#				my $query_pos 	= 0;
#
#				for (my $j = 0; $j < scalar(@query_seq); $j++) { 
#					if ($query_seq[$j] =~ /[A-Z]/ && $hit_seq[$j] eq '-') { 
#						$query_pos++;
#					}elsif($query_seq[$j] eq '-' && $hit_seq[$j] =~ /[A-Z]/) { 
#						$hit_pos++;
#					}elsif($query_seq[$j] =~ /[A-Z]/ && $hit_seq[$j] =~ /[A-Z]/) { 
#						if ($$query_seg_seqid_aref[$query_pos] && $$hit_seg_seqid_aref[$hit_pos]) { 
#							$query_map{$$query_seg_seqid_aref[$query_pos] } = $$hit_seg_seqid_aref[$hit_pos];
#							$hit_map{$$hit_seg_seqid_aref[$hit_pos]} = $$query_seg_seqid_aref[$query_pos];
#							$hit_chain_map{$$hit_seg_seqid_aref[$hit_pos]} = $$query_seg_chain_aref[$query_pos];
#						}
#						#print "MAPBUG $hit_pos $$hit_seg_seqid_aref[$hit_pos] $hit_start $query_pos $$query_seg_seqid_aref[$query_pos] $query_start $$query_seg_chain_aref[$query_pos]\n";
#						$query_pos++;
#						$hit_pos++;
#					}else{
#						#Foul
#					}
#				}
#			}
#				
#			#This is sloppy. Build this out of a chainwise index? 
#			#my $hit_reference_domain_structure_nodes = $ref_xml_doc->findnodes(qq{//structure[\@pdb_id='$hit_pdb'][\@chain_id='$hit_chain']});
#			#if ($hit_reference_domain_structure_nodes->size() == 0) { 
#			#	print "WARNING! No nodes found in reference for $hit_pdb $hit_chain\n";
#			#	next;
#			#}
#
#			#foreach my $node ($hit_reference_domain_structure_nodes->get_nodelist()) { 
#
#			foreach my $hit_domain_uid ($$ref_chain_domains{$hit_pdb}{$hit_chain}) { 
#
#				my $ref_domain_id = $$ref_range_cache{$hit_domain_uid}{ecod_domain_id};
#
#				my $hit_domain_id_node	= $domain_xml_doc->createElement('hit_domain');
#				$hit_domain_id_node->setAttribute('ecod_domain_id', $ref_domain_id);
#				$hit_domain_id_node->setAttribute('reference', $reference);
#				$hit_domain_id_node->setAttribute('uid', $hit_domain_uid);
#				
#				my $ref_seqid_range = multi_chain_ungap_range($$ref_range_cache{$hit_domain_uid}{seqid_range}, $GAP_TOL);
#				
#				if ($DEBUG) { 
#					print "DEBUG: $sub: chblastp_ref_domain $ref_domain_id $ref_seqid_range\n";
#				}
#				#This section needs to go first... be replaced.
#
#				#my $new_ref_node	= $ref_hit_domain_node->cloneNode(1);
#				#$hit_domain_id_node->appendChild($new_ref_node);
#				
#				my $domain_node		= $domain_xml_doc->createElement('domain');
#				$domain_node->appendChild($hit_domain_id_node);
#				#my $domain_id		= "d$query_pdb$query_chain$domain_count";
#
#				my $method_node	= $domain_xml_doc->createElement('method');
#				my $method	= 'chblastp';
#				$method_node->appendTextNode($method);
#				$domain_node->appendChild($method_node);
#
#				#REPLACE WITH HIT DOMAIN DERIVED RANGE
#				my @ref_segs = split(/,/, $ref_seqid_range);
#				my @query_struct_seqid_range;
#				my @query_struct_chain_range;
#				my @hit_struct_seqid_range;
#				my $warning = 0;
#				foreach my $ref_seg (@ref_segs) { 
#					#print "REFSEG $ref_seg\n";
#					if ($ref_seg =~ /(\w+):(\-?\d+\-\-?\d+)/) { 
#						my $ref_seg_chain = $1;
#						if ($ref_seg_chain ne $hit_chain) { 
#							print "WARNING! ref seg chain $ref_seg_chain is not the same as expected hit chain $hit_chain, multi chain problems...\n";
#						}
#						my $ref_seg_range = $2;
#						my $ref_seg_range_aref = range_expand($ref_seg_range);
#						foreach my $ref_seqid (@$ref_seg_range_aref) { 
#							if ($hit_map{$ref_seqid}) { 
#								print "MESS $ref_seqid $hit_map{$ref_seqid} $ref_seg_chain \n";
#								push (@query_struct_seqid_range, $hit_map{$ref_seqid});
#								#push (@query_struct_chain_range, $ref_seg_chain);
#								push (@query_struct_chain_range, $hit_chain_map{$ref_seqid});
#							}else{
#								#print "WARNING! No map for $ref_seqid\n";
#								$warning++;
#							}
#							push (@hit_struct_seqid_range, $ref_seqid);
#						}
#					}else{
#						print "WARNING! $sub: Ref seg looks weird $ref_seg\n";
#					}
#				}
#				if ($warning) { 
#					print "WARNING! $warning maps failed\n";
#				}
#				my $struct_seqid_range_node	= $domain_xml_doc->createElement('struct_seqid_range');
#				#$struct_seqid_range_node->appendTextNode($query_struct_seqid_range);	
#				if (scalar(@query_struct_seqid_range) == 0) { 
#					print "WARNING! $sub: No struct seqid for $query_pdb $query_chains chblastp, serious problem...\n";
#					next;
#				}
#				my $query_struct_seqid_range = multi_chain_rangify(\@query_struct_seqid_range, \@query_struct_chain_range);
#				$struct_seqid_range_node->appendTextNode($query_struct_seqid_range);	
#				$domain_node->appendChild($struct_seqid_range_node);
#
#				my $ungapped_seqid_range = multi_chain_ungap_range($query_struct_seqid_range, $$global_opt{gap_tol});
#				my ($ungapped_seqid_range_aref, $ungapped_chain_range_aref) = multi_chain_range_expand($ungapped_seqid_range);
#				my $ungapped_seqid_range_node	= $domain_xml_doc->createElement('ungapped_seqid_range');
#				$ungapped_seqid_range_node->setAttribute('gap_tolerance', $$global_opt{gap_tol});
#				$ungapped_seqid_range_node->appendTextNode($ungapped_seqid_range);
#				$domain_node->appendChild($ungapped_seqid_range_node);
#
#				#REPLACE WITH HIT DOMAIN RANGE (THAT OVERLAPS ALLGNED REGION)
#				my $hit_struct_seqid_range_node	 = $domain_xml_doc->createElement('hit_struct_seqid_range');
#				#$hit_struct_seqid_range_node->appendTextNode($hit_hit_reg);
#				my $hit_struct_seqid_range;
#				$hit_struct_seqid_range = scopify_range(rangify(@hit_struct_seqid_range), $hit_chain);
#
#				$hit_struct_seqid_range_node->appendTextNode($hit_struct_seqid_range);
#				$domain_node->appendChild($hit_struct_seqid_range_node);
#
#				my $pdb_range = multi_chain_pdb_rangify(\@query_struct_seqid_range,  $query_pdbnum_aref,  \@query_struct_chain_range);
#				my $pdb_range_node	= $domain_xml_doc->createElement('pdb_range');
#				$pdb_range_node->appendTextNode($pdb_range);
#				$domain_node->appendChild($pdb_range_node);
#
#				my $ungapped_pdb_range = multi_chain_pdb_rangify($ungapped_seqid_range_aref, $query_pdbnum_aref, $ungapped_chain_range_aref );
#				my $ungapped_pdb_range_node	= $domain_xml_doc->createElement('ungapped_pdb_range');
#				$ungapped_pdb_range_node->setAttribute('gap_tolerance', $$global_opt{gap_tol});
#				$ungapped_pdb_range_node->appendTextNode($ungapped_pdb_range);
#				$domain_node->appendChild($ungapped_pdb_range_node);
#
#				my ($range_str, $chain_str) = scop_range_split($ungapped_pdb_range);
#				my $domain_id_chain_str;
#				if ($chain_str =~ /,/) { 
#					$domain_id_chain_str = ".";
#				}else{
#					$domain_id_chain_str = $chain_str;
#				}
#
#				my $domain_id	= "e". lc( $query_pdb) . "$domain_id_chain_str$domain_count";
#				$domain_node->setAttribute('domain_id', $domain_id);
#
#				my $blast_score_node	= $domain_xml_doc->createElement('blast_eval');
#				$blast_score_node->appendTextNode($evalue);
#				$domain_node->appendChild($blast_score_node);
#
#				$domain_list_node->appendChild($domain_node);
#				#REPLACE WITH DOMAIN BASED RANGE
#
#				my ($query_struct_seqid_aref, $query_struct_chain_aref) = multi_chain_range_expand($query_struct_seqid_range);
#
#				#printf "check1 %i %i\n", scalar(@unused_seq), scalar(@unused_chain);
#				multi_chain_range_exclude($query_struct_seqid_aref, $query_struct_chain_aref, \@unused_seq, \@unused_chain, \%sort_lookup);
#				#printf "check2 %i %i\n", scalar(@unused_seq), scalar(@unused_chain);
#				multi_chain_range_include($query_struct_seqid_aref, $query_struct_chain_aref, \@used_seq, \@used_chain, \%sort_lookup);
#
#				$domain_count++;
#			}
#		}
#	}
#	#BLAST
#	my $blastp_run_XPath = q{//blast_summ_doc/blast_summ/blast_run[@program='blastp']};
#	my $blastp_run_nodes = $blast_summ_xml_doc->findnodes($blastp_run_XPath);
#	my $blastp_run_nodes_size = $blastp_run_nodes->size();
#	if ($blastp_run_nodes_size != 1) { 
#		print "ERROR! $sub: Odd number of blastp run nodes found ($blastp_run_nodes_size), aborting...\n";
#
#	}
#	my $blastp_hit_nodes = $blastp_run_nodes->get_node(1)->findnodes('hits/hit');
#	foreach my $hit_node ($blastp_hit_nodes->get_nodelist() ) { 
#
#		my $hit_num 	= $hit_node->findvalue('@num');
#		my $hit_domain_id = $hit_node->findvalue('@domain_id');
#		my $hit_domain_uid = $hit_node->findvalue('@uid');
#		my $hsp_count 	= $hit_node->findvalue('@hsp_count');
#		my $evalues 	= $hit_node->findvalue('@evalues');
#
#		$hit_domain_id =~ /(e|d)(\w{4})(\w+)/; #Define overall hit domain regexp elsewhere, time to obsolete SCOP domains in hits.
#		
#		my $dtype = $1;
#		my $hit_pdb = $2;
#		my $hit_chain = $3;
#		if ($DEBUG) { 
#		}
#		#if ($hit_chain eq '.') { next } 
#		if ($dtype eq 'd') {
#			$hit_chain = uc($hit_chain);
#		}
#
#		my $evalue = 'Unk';
#		if ($hsp_count == 1) { 
#			$evalue = $evalues;
#		}
#
#		#This SHOUUULD be seqid_reange
#		#my $hit_query_seqid_range = $hit_node->findvalue('text()');
#		#It's not anymore...
#		#my $hit_query_reg	= $hit_node->findvalue('text()');
#		my $hit_query_reg	= $hit_node->findvalue('query_reg');
##		if (!$query_struct_seqid_aref) { die "bad qssa\n";}
##		my ($hit_query_struct_seqid_aref, $hit_query_struct_chain_aref);
##		if ($input_mode eq 'struct_seqid') { 
##			#$hit_query_struct_seqid_aref = struct_region(range_expand($hit_query_reg), $query_struct_seqid_aref);
##			($hit_query_struct_seqid_aref, $hit_query_struct_chain_aref) = multi_chain_struct_region(range_expand($hit_query_reg), $query_struct_seqid_aref, $query_struct_chain_aref);
##		}elsif($input_mode eq 'seqid') { 
##			#$hit_query_struct_seqid_aref = range_expand($hit_query_reg);
##			($hit_query_struct_seqid_aref, $hit_query_struct_chain_aref) = multi_chain_struct_region(range_expand($hit_query_reg), $query_seqid_aref, $query_chain_aref)
##		}else{
##			die "ERROR! input_mode $input_mode not recognized\n";
##		}
#		#my $hit_query_struct_seqid_range = rangify(@$hit_query_struct_seqid_aref);
#
#		my ($hit_query_struct_seqid_aref, $hit_query_struct_chain_aref)  =
#			range_decompose($hit_query_reg, 
#					$query_struct_seqid_aref, 
#					$query_seqid_aref, 
#					$query_struct_chain_aref, 
#					$query_chain_aref, 
#					$input_mode);
#		my $hit_query_struct_seqid_range =
#			multi_chain_rangify($hit_query_struct_seqid_aref, 
#						$hit_query_struct_chain_aref);		
#
#
#		my $hit_hit_reg		= $hit_node->findvalue('hit_reg');
#
#		if ($DEBUG) { 
#			print "DEBUG $sub: blastp $hit_num $hit_domain_id $hit_query_struct_seqid_range\n";
#		}
#		#Assign if overlaps 
#
#		my $unused_coverage 	= multi_chain_region_coverage($query_struct_seqid_aref, $query_struct_chain_aref, \@unused_seq, \@unused_chain);
#		my $used_coverage 	= multi_chain_region_coverage($query_struct_seqid_aref, $query_struct_chain_aref, \@used_seq, \@used_chain);
#		
#		if ($unused_coverage < 0.05 || scalar(@unused_seq) < 10) { 
#			last;
#		}
#
#		my $query_new_coverage 	= multi_chain_region_coverage($hit_query_struct_seqid_aref, $hit_query_struct_chain_aref, \@unused_seq, \@unused_chain);
#		my $query_used_coverage = multi_chain_region_coverage($hit_query_struct_seqid_aref, $hit_query_struct_chain_aref, \@used_seq, \@used_chain);
#		if ($DEBUG) { 
#			print "DEBUG $sub: uc $used_coverage uuc $unused_coverage qnc $query_new_coverage qoc $query_used_coverage\n";
#		}
#
#		if ($query_new_coverage > $$global_opt{new_coverage_threshold} && $query_used_coverage < $$global_opt{old_coverage_threshold}) { 
#			print "DEFINE blastp: $sub: $hit_domain_id $hit_query_struct_seqid_range\n";
#
#			my ($hit_seqid_aref, $hit_struct_seqid_aref, $hit_pdbnum_aref, $hit_asym_id, $hit_chain_aref, $hit_struct_chain_aref);
#			if ($hit_chain eq '.') { 
#				my $domain_range = $$ref_range_cache{$hit_domain_uid}{seqid_range};
#				my ($chain_str, $chain_range) = scop_range_split($domain_range);
#				my @chains = split(",", $chain_str);
#				($hit_seqid_aref, $hit_struct_seqid_aref, $hit_pdbnum_aref, $hit_asym_id, $hit_chain_aref, $hit_struct_chain_aref) = pdbml_mc_seq_parse($hit_pdb, \@chains);
#			}else{
#				($hit_seqid_aref, $hit_struct_seqid_aref, $hit_pdbnum_aref, $hit_asym_id) = pdbml_seq_parse($hit_pdb, $hit_chain);
#			}
#			if (!$hit_seqid_aref || !$hit_struct_seqid_aref) { 
#				print "WARNING! $sub: hit parse failure\n";
#				next;
#			}
#
#			if (!$hit_hit_reg) { die "ERROR $sub: hit hit reg failure\n"}
#			my ($hit_reg_struct_seqid_aref, $hit_reg_struct_chain_aref);
#			if ($hit_chain eq '.') { 
#				($hit_reg_struct_seqid_aref, $hit_reg_struct_chain_aref) = multi_chain_struct_region(range_expand($hit_hit_reg), $hit_struct_seqid_aref, $hit_struct_chain_aref);
#			}else{
#				$hit_reg_struct_seqid_aref = struct_region(range_expand($hit_hit_reg), $hit_struct_seqid_aref);
#			}
#
#			my $domain_node		= $domain_xml_doc->createElement('domain');
#			#my $domain_id		= "d$query_pdb$query_chain$domain_count";
#			#my $domain_id	= "e". lc( $query_pdb) . "$query_chains_str$domain_count";
#			#$domain_node->setAttribute('domain_id', $domain_id);
#
#			my $method_node	= $domain_xml_doc->createElement('method');
#			my $method	= 'blastp';
#			$method_node->appendTextNode($method);
#			$domain_node->appendChild($method_node);
#
#			my $struct_seqid_range_node	= $domain_xml_doc->createElement('struct_seqid_range');
#			$struct_seqid_range_node->appendTextNode($hit_query_struct_seqid_range);	
#			$domain_node->appendChild($struct_seqid_range_node);
#
#			my $ungapped_seqid_range = multi_chain_ungap_range($hit_query_struct_seqid_range, $$global_opt{gap_tol});
#			my ($ungapped_seqid_range_aref, $ungapped_seqid_chain_aref) = multi_chain_range_expand($ungapped_seqid_range);
#			my $ungapped_seqid_range_node	= $domain_xml_doc->createElement('ungapped_seqid_range');
#			$ungapped_seqid_range_node->setAttribute('gap_tolerance', $$global_opt{gap_tol});
#			$ungapped_seqid_range_node->appendTextNode($ungapped_seqid_range);
#			$domain_node->appendChild($ungapped_seqid_range_node);
#
#			my $hit_struct_seqid_range_node	 = $domain_xml_doc->createElement('hit_struct_seqid_range');
#			if ($hit_chain eq '.') { 
#				$hit_struct_seqid_range_node->appendTextNode(multi_chain_rangify($hit_reg_struct_seqid_aref, $hit_reg_struct_chain_aref));
#			}else{
#				$hit_struct_seqid_range_node->appendTextNode(rangify(@$hit_reg_struct_seqid_aref));
#			}
#
#			$domain_node->appendChild($hit_struct_seqid_range_node);
#
#			#my $pdb_range = pdb_rangify($hit_query_struct_seqid_aref, $query_pdbnum_aref);
#			my $pdb_range = multi_chain_pdb_rangify($hit_query_struct_seqid_aref, $query_pdbnum_aref, $hit_query_struct_chain_aref);
#			my $pdb_range_node	= $domain_xml_doc->createElement('pdb_range');
#			$pdb_range_node->appendTextNode($pdb_range);
#			$domain_node->appendChild($pdb_range_node);
#
#			#my $ungapped_pdb_range = pdb_rangify($ungapped_seqid_range_aref, $query_pdbnum_aref);
#			my $ungapped_pdb_range = multi_chain_pdb_rangify($ungapped_seqid_range_aref, $query_pdbnum_aref, $ungapped_seqid_chain_aref);
#			my $ungapped_pdb_range_node	= $domain_xml_doc->createElement('ungapped_pdb_range');
#			$ungapped_pdb_range_node->setAttribute('gap_tolerance', $$global_opt{gap_tol});
#			$ungapped_pdb_range_node->appendTextNode($ungapped_pdb_range);
#			$domain_node->appendChild($ungapped_pdb_range_node);
#
#			my ($range_str, $chain_str) = scop_range_split($ungapped_pdb_range);
#			my $domain_id_chain_str;
#			if ($chain_str =~ /,/) { 
#				$domain_id_chain_str = ".";
#			}else{
#				$domain_id_chain_str = $chain_str;
#			}
#
#			my $domain_id	= "e". lc( $query_pdb) . "$domain_id_chain_str$domain_count";
#			$domain_node->setAttribute('domain_id', $domain_id);
#
#			my $hit_domain_id_node	= $domain_xml_doc->createElement('hit_domain');
#			$hit_domain_id_node->setAttribute('ecod_domain_id', $hit_domain_id);
#			$hit_domain_id_node->setAttribute('reference', $reference);
#			$hit_domain_id_node->setAttribute('uid', $$ref_domain_uid_lookup{$hit_domain_id});
#			$domain_node->appendChild($hit_domain_id_node);
#
#			my $blast_score_node	= $domain_xml_doc->createElement('blast_eval');
#			$blast_score_node->appendTextNode($evalue);
#			$domain_node->appendChild($blast_score_node);
#
#			$domain_list_node->appendChild($domain_node);
#			#range_exclude(range_expand($hit_query_struct_seqid_range), \@unused_seq);
#			my ($hit_query_struct_seqid_aref, $hit_query_struct_chain_aref) = multi_chain_range_expand($hit_query_struct_seqid_range);
#			multi_chain_range_exclude($hit_query_struct_seqid_aref, $hit_query_struct_chain_aref, \@unused_seq, \@unused_chain, \%sort_lookup);
#			#range_include(range_expand($hit_query_struct_seqid_range), \@used_seq);
#			multi_chain_range_include($hit_query_struct_seqid_aref, $hit_query_struct_chain_aref, \@used_seq, \@used_chain, \%sort_lookup);
#
#			$domain_count++;
#			#Self_comparison check
#			for (my $i = 0; $i < scalar(@self_comps); $i++) { 
#
#				my $self_comp_aligner 	= $self_comps[$i]{aligner};
#				my $self_comp_query_reg	= $self_comps[$i]{query_reg};
#				my $self_comp_hit_reg	= $self_comps[$i]{hit_reg};
#				
#				my $sc_query_reg_aref 	= range_expand($self_comp_query_reg);
#				my $sc_hit_reg_aref	= range_expand($self_comp_hit_reg);
#				
#				my $sc_query_coverage_1 = region_coverage($sc_query_reg_aref, $ungapped_seqid_range_aref);
#				my $sc_query_coverage_2 = region_coverage($ungapped_seqid_range_aref, $sc_query_reg_aref,);
#				my $sc_hit_coverage_1 = region_coverage($sc_hit_reg_aref, \@unused_seq);
#				#my $sc_hit_coverage_2 = region_coverage(\@unused_seq, $sc_hit_reg_aref);
#				
#				if ($DEBUG) { 
#					print "DEBUG blastp SC $i $self_comp_aligner qr $self_comp_query_reg hr $self_comp_hit_reg sqc1 $sc_query_coverage_1 sqc2 $sc_query_coverage_2 shc1 $sc_hit_coverage_1\n";
#				}
#
#				if ($sc_query_coverage_1 > 0.9 && $sc_query_coverage_2 > 0.9 && $sc_hit_coverage_1 > 0.9 ) { 
#					my $pre_used = rangify(@used_seq);
#					my $pre_unused = rangify(@unused_seq);
#					print "DEFINE blastp self_comp $sub: $self_comp_aligner qr $self_comp_query_reg hr $self_comp_hit_reg pu $pre_used pun $pre_unused\n";
#
#					my $sc_domain_node 	= $domain_xml_doc->createElement('domain');
#
#					my $domain_id	= 'e' . lc($query_pdb) . "$query_chains$domain_count";
#					$sc_domain_node->setAttribute('domain_id', $domain_id);
#
#					my $method_node	= $domain_xml_doc->createElement('method');
#					my $method 	= 'blastp_self_comp';
#					$method_node->appendTextNode($method);
#					$method_node->setAttribute('self_comp_aligner', $self_comp_aligner);
#					$sc_domain_node->appendChild($method_node);
#
#					my $struct_seqid_range_node	= $domain_xml_doc->createElement('struct_seqid_range');
#					$struct_seqid_range_node->appendTextNode($self_comp_hit_reg);
#					$sc_domain_node->appendChild($struct_seqid_range_node);
#
#					my $ungapped_seqid_range	= ungap_range($self_comp_hit_reg, $$global_opt{gap_tol});
#					my $ungapped_seqid_range_aref	= range_expand($ungapped_seqid_range);
#					my $ungapped_seqid_range_node	= $domain_xml_doc->createElement('ungapped_seqid_range');
#					$ungapped_seqid_range_node->setAttribute('gap_tolerance', $$global_opt{gap_tol});
#					$ungapped_seqid_range_node->appendTextNode($ungapped_seqid_range);
#					$sc_domain_node->appendChild($ungapped_seqid_range_node);
#
#					my $hit_struct_seqid_range_node = $domain_xml_doc->createElement('hit_struct_seqid_range');
#					$hit_struct_seqid_range_node->appendTextNode($self_comp_query_reg);
#					$sc_domain_node->appendChild($hit_struct_seqid_range_node);
#					
#					my $pdb_range	= pdb_rangify($sc_hit_reg_aref, $query_pdbnum_aref);
#					my $pdb_range_node	= $domain_xml_doc->createElement('pdb_range');
#					$pdb_range_node->appendTextNode($pdb_range);
#					$sc_domain_node->appendChild($pdb_range_node);
#
#					my $ungapped_pdb_range	= pdb_rangify($ungapped_seqid_range_aref, $query_pdbnum_aref);
#					my $ungapped_pdb_range_node	= $domain_xml_doc->createElement('ungapped_pdb_range');
#					$ungapped_pdb_range_node->setAttribute('gap_tolerance', $$global_opt{gap_tol});
#					$ungapped_pdb_range_node->appendTextNode($ungapped_pdb_range);
#					$sc_domain_node->appendChild($ungapped_pdb_range_node);
#
#					if ($self_comp_aligner eq 'dali') { 
#						my $dali_score_node	= $domain_xml_doc->createElement('dali_zscore');
#						$dali_score_node->appendTextNode($self_comps[$i]{zscore});
#						$sc_domain_node->appendChild($dali_score_node);
#					}elsif($self_comp_aligner eq 'hhrepid') { 
#						my $hh_prob_node	= $domain_xml_doc->createElement('hhrepid_prob');
#						$hh_prob_node->appendTextNode($self_comps[$i]{prob});
#						$sc_domain_node->appendChild($hh_prob_node);
#					}else{
#						print "WARNING! self comp aligner $self_comp_aligner unknown\n";
#					}
#
#					$domain_list_node->appendChild($sc_domain_node);
#					range_include($ungapped_seqid_range_aref, \@used_seq);
#					range_exclude($ungapped_seqid_range_aref, \@unused_seq);
#
#					my $post_used = rangify(@used_seq);
#					my $post_unused = rangify(@unused_seq);
#					
#					$domain_count++;
#					
#				}
#			}
#		}
#	}
#	#my $final_coverage = region_coverage($query_struct_seqid_aref, \@used_seq);
#	my $final_coverage = multi_chain_region_coverage($query_struct_seqid_aref, $query_struct_chain_aref, \@used_seq, \@used_chain);
#	my $chain_domain_coverage_node = $domain_xml_doc->createElement('chain_domain_coverage');
#	$domain_doc_node->appendChild($chain_domain_coverage_node);
#	$chain_domain_coverage_node->appendTextNode($final_coverage);
#
#	my $domain_doc_string = $domain_xml_doc->toString(1);
#
#	open(OUT, ">$domain_xml_fn") or die "ERROR! $sub: Could not open $domain_xml_fn for writing:$!\n";
#	print OUT $domain_doc_string;
#	close OUT;
#
#	return 1;
#
#
#}


sub domain_partition { 
	my $sub = 'domain_partition';

	my ($query_pdb, $query_chain, $dir, $input_mode, $reference, $domain_xml_fn, $ref_range_cache, $concise, $global_opt, $ref_chain_domains, $ref_domain_uid_lookup) = @_;
	print "dxml: $domain_xml_fn\n";

	if (-f $domain_xml_fn && !$$global_opt{force_overwrite}) { 
		print "WARNING! $sub: skipping $query_pdb, $query_chain, set FORCE_OVERWRITE bit to override...\n";
		return;
	}

	#Return the seqids and struct seqids of the query pdb_chain
	my ($query_seqid_aref, $query_struct_seqid_aref, $query_pdbnum_aref, $query_asym_id) = 
			pdbml_seq_parse($query_pdb, $query_chain);

	#Return annotations for the query pdb_chain
	my ($struct_href, 
		$exptl_href, 
		$entity_href, 
		$target_entity_id, 
		$struct_ref_href, 
		$citation_href, 
		$citation_author_href, 
		$entity_src_nat_href, 
		$entity_src_gen_href, 
		$entity_src_syn_href) = pdbml_annot_parse($query_pdb, $query_chain);

	my $pdb_chain = $query_pdb . "_" . $query_chain;

	my $src_method = $$entity_href{$target_entity_id}{src_method};

	if (!$query_seqid_aref) { 
		print "WARNING! $sub: $query_pdb not found! skipping...\n";
		return 0;
	}

	#if ($recurse) { 
	#	$query_seqid_aref = isect(range_expand($recurse_range), $query_seqid_aref);
	#	$query_struct_seqid_aref = isect(range_expand($recurse_range), $query_struct_seqid_aref);
	#}

	my $query_range_str 		= rangify(@$query_seqid_aref);
	my $query_struct_range_str 	= rangify(@$query_struct_seqid_aref);

	#Changed query sequence to be ungapped (i.e. may contain small unstructured regions, extremely important that GAP_TOL is the same for query generation as domain partition).
	my $query_ungapped_struct_range_str 	= ungap_range($query_struct_range_str, $$global_opt{gap_tol});
	my @gapped_query_struct_seqid 	= @$query_struct_seqid_aref; #Need to retain this for dali_fillin (query pdbs are not ungapped)
	$query_struct_seqid_aref 		= range_expand($query_ungapped_struct_range_str);

	if (!$query_struct_range_str) { 
		print "WARNING! $sub: $query_pdb $query_chain NO structured range $query_struct_range_str\n";
	}

	my @unused_seq = $input_mode eq 'struct_seqid' ? @$query_struct_seqid_aref : @$query_seqid_aref;
	my @used_seq;

	if ($DEBUG) { 
		print "DEBUG $sub: $query_pdb $query_chain $dir $query_range_str $query_struct_range_str $domain_xml_fn\n";
	}
	my $blast_summ_fn;
	#if ($recurse) { 
	#	my $seqid_range = $recurse_range;
	#	$seqid_range =~ s/\,/_/g;
	#	$blast_summ_fn = "$dir/${query_pdb}_${query_chain}_${seqid_range}.$reference.blast_summ.xml";
	#}

	if( $concise ) { 
		$blast_summ_fn = "$dir/$pdb_chain.$reference.blast_summ.concise.xml";
	}else{
		$blast_summ_fn = "$dir/$pdb_chain.$reference.blast_summ.xml";
	}
	if (!-f $blast_summ_fn) { 
		if (-f "$dir/$pdb_chain.$reference.blast_summ.concise.xml") { 
			$blast_summ_fn = "$dir/$pdb_chain.$reference.blast_summ.concise.xml";
		}else{
			print "WARNING! $sub: BLAST summary for $query_pdb, $query_chain, $dir ($blast_summ_fn) not found\n";
			next;
		}
	}

	my $domain_count	= 1; #Used for generating new domain ids

	my $domain_xml_doc	= XML::LibXML->createDocument;
	my $domain_doc_node	= $domain_xml_doc->createElement('chain_domains_set_top');
	$domain_xml_doc->setDocumentElement($domain_doc_node);

	$domain_doc_node->setAttribute('pdb_id', $query_pdb);
	$domain_doc_node->setAttribute('chain_id', $query_chain);
	$domain_doc_node->setAttribute('mode','blast');
	$domain_doc_node->setAttribute('domain_prefix', $DOMAIN_PREFIX);

	my $lc_query_pdb = lc($query_pdb);

	if ($src_method eq 'syn') { 
		print "WARNING! $query_pdb, $query_chain is a known synthetic product\n";
		$domain_doc_node->setAttribute('known_synthetic', 'true');
	}

	#Do some work for coil detection
	my $coils_fn = "$dir/$pdb_chain.coils.xml";
	my @coils;
	if (-f $coils_fn) { 
		print "Found coils definition for $pdb_chain\n";
		my $coils_xml_doc = xml_open($coils_fn);
		foreach my $coil_node ($coils_xml_doc->findnodes('//coil')) { 
			my $coil_seqid_range	= $coil_node->findvalue('coil_seqid_range');
			next unless $coil_seqid_range =~ /\d+/;
			my $method	= $coil_node->findvalue('@method');
			my %coil_entry;
			$coil_entry{method} = $method;	
			my ($range_str, $chain_str) = scop_range_split($coil_seqid_range);
			$coil_entry{seqid_range} = $range_str;
			push (@coils, \%coil_entry);
		}
	}

	#Build the stub of the domains xml file
	my $domain_list_node	= $domain_xml_doc->createElement('domain_list');
	$domain_doc_node->appendChild($domain_list_node);

	#Open the sequence summary search file (i.e. the "blast_summ" file)
	my $blast_summ_xml_doc = xml_open($blast_summ_fn);

	my $blast_summ_pdb 		= $blast_summ_xml_doc->findvalue('//blast_summ_doc/blast_summ/@pdb');
	my $blast_summ_chain	= $blast_summ_xml_doc->findvalue('//blast_summ_doc/blast_summ/@chain');

	if ($query_pdb ne $blast_summ_pdb || $query_chain ne $blast_summ_chain) { 
		#I don't think I've seen this error in over 4 years. 08/2015 rds
		print "ERORR! hora file/blast file pdb mismatch $query_pdb/$blast_summ_pdb $query_chain/$blast_summ_chain\n";
		die;
	}
	#Self_comparison read:
	my $self_comp_aref = read_self_comps($blast_summ_xml_doc);
	
	#Chainwise BLAST
	find_chainwise_blast_domains(
		$blast_summ_xml_doc,
		$domain_xml_doc,
		$ref_chain_domains,
		$ref_range_cache,
		$ref_domain_uid_lookup,
		$query_struct_seqid_aref,
		$query_seqid_aref,
		$query_pdbnum_aref,
		\@used_seq,
		\@unused_seq,
		$input_mode,
		$global_opt,
		$reference,
		\@coils,
	);

	#BLAST
	find_domainwise_blast_domains(
		$blast_summ_xml_doc,
		$domain_xml_doc,
		$ref_chain_domains,
		$ref_domain_uid_lookup,
		$ref_range_cache,
		$query_struct_seqid_aref,
		$query_seqid_aref,
		$query_pdbnum_aref,
		$self_comp_aref,
		\@used_seq,
		\@unused_seq,
		$input_mode,
		$global_opt,
		$reference,
		\@coils,
	);

	#PSIBLAST
	unless ($NO_PSIBLAST) { 
		find_psiblast_domains (
			$blast_summ_xml_doc,
			$domain_xml_doc,
			$ref_chain_domains,
			$ref_domain_uid_lookup,
			$ref_range_cache,
			$query_struct_seqid_aref,
			$query_seqid_aref,
			$query_pdbnum_aref,
			$self_comp_aref,
			\@used_seq,
			\@unused_seq,
			$input_mode,
			$global_opt,
			$reference,
			\@coils,
		);
	}
	
	#HHsearch
	find_hhsearch_domains (
		$blast_summ_xml_doc,
		$domain_xml_doc,
		$ref_chain_domains,
		$ref_domain_uid_lookup,
		$ref_range_cache,
		$query_struct_seqid_aref,
		$query_seqid_aref,
		$query_pdbnum_aref,
		$self_comp_aref,
		\@used_seq,
		\@unused_seq,
		$input_mode,
		$global_opt,
		$reference,
		\@coils,
	);

	#DALI fillin 
	print "dir: $dir\n";
	my $dali_summ_fn = generate_dali_fillin ( 
		$blast_summ_xml_doc,		#1
		$domain_xml_doc,		#2
		$ref_chain_domains,		#3
		$ref_domain_uid_lookup,		#4
		$ref_range_cache,		#5
		#$query_struct_seqid_aref,	#6	
		\@gapped_query_struct_seqid,
		$query_seqid_aref,		#7
		$query_pdbnum_aref,		#8
		$query_asym_id,	
		\@used_seq,			#9
		\@unused_seq,			#10
		$input_mode,			#11
		$global_opt,			#12
		$reference,			#13
		$dir,				#14
	);
	if (-f $dali_summ_fn) { 
		my $dali_summ_xml_doc = xml_open($dali_summ_fn);
		find_dali_fillin_domains (
			$domain_xml_doc,
			$dali_summ_xml_doc,
			$ref_chain_domains,
			$ref_domain_uid_lookup,
			$ref_range_cache,
			#$query_struct_seqid_aref,
			\@gapped_query_struct_seqid,
			$query_pdbnum_aref,
			$query_asym_id,
			\@used_seq,
			\@unused_seq,
			$input_mode,
			$global_opt,
			$reference,
			$dir
		);
	}
	
	if (@coils) { 
		define_coils(
			$domain_xml_doc, 
			\@coils,
			\@unused_seq,
			\@used_seq,
			$global_opt,
			$query_pdb,
			$query_chain,
			$query_pdbnum_aref
			);
	}

	my $final_coverage 		= region_coverage($query_struct_seqid_aref, \@used_seq);
	$domain_doc_node->appendTextChild('chain_domain_coverage', $final_coverage);

	xml_write($domain_xml_doc, $domain_xml_fn);

	return 1;
}

sub struct_domain_partition { 
	...
}

sub read_self_comps { 
	my $sub = 'read_self_comps';

	my ($blast_summ_xml_doc) = @_;
	my $self_comp_run_XPath = q{//blast_summ_doc/blast_summ/self_comp_run};
	my $self_comp_run_nodes = $blast_summ_xml_doc->findnodes($self_comp_run_XPath);
	my $self_comp_run_nodes_size = $self_comp_run_nodes->size();
	my @self_comps;
	if ($self_comp_run_nodes_size < 1) { 
		print "WARNING! $sub: No self_comp nodes found\n";
	}else{
		#my $self_comp_hit_nodes = $self_comp_run_nodes->get_node(1)->findnodes('hits/hit');
		my $self_comp_hit_nodes = $blast_summ_xml_doc->findnodes('//blast_summ_doc/blast_summ/self_comp_run/repeat_set_list/repeat_hits/hit');
		my $ci = 0;
		foreach my $hit_node ($self_comp_hit_nodes->get_nodelist() ) { 

			#my $aligner 	= $hit_node->findvalue('@aligner');
			my $aligner = $hit_node->parentNode->parentNode->parentNode->findvalue('@programs');

			if ($aligner eq 'dali') { 
				my $zscore	= $hit_node->findvalue('@z_score');
				$self_comps[$ci]{zscore} = $zscore;
			}elsif ($aligner eq 'hhrepid') { 
				my $prob	= $hit_node->findvalue('@prob');
				$self_comps[$ci]{prob}	= $prob;
			}else{
				print "WARNING! Unknown self_comparison aligner $aligner, skipping...\n";
				next;
			}

			my $query_reg	= $hit_node->findvalue('query_reg');
			my $hit_reg	= $hit_node->findvalue('hit_reg');
			
			if (! $query_reg =~ /\d+/) { 
				next;
			}

			if (!$hit_reg =~ /\d+/) { 
				next;
			}

			$self_comps[$ci]{aligner}	= $aligner;
			
			$self_comps[$ci]{query_reg}	= $query_reg;
			$self_comps[$ci]{hit_reg}	= $hit_reg;
			$ci++;
		}
		$self_comp_hit_nodes = $blast_summ_xml_doc->findnodes('//blast_summ_doc/blast_summ/self_comp_run/hits/hit');
		foreach my $hit_node ($self_comp_hit_nodes->get_nodelist() ) { 

			#my $aligner 	= $hit_node->findvalue('@aligner');
			my $aligner = $hit_node->parentNode->parentNode->findvalue('@programs');

			if ($aligner eq 'dali') { 
				my $zscore	= $hit_node->findvalue('@z_score');
				$self_comps[$ci]{zscore} = $zscore;
			}elsif ($aligner eq 'hhrepid') { 
				my $prob	= $hit_node->findvalue('@prob');
				$self_comps[$ci]{prob}	= $prob;
			}else{
				print "WARNING! Unknown self_comparison aligner $aligner, skipping...\n";
				next;
			}

			my $query_reg	= $hit_node->findvalue('query_reg');
			my $hit_reg	= $hit_node->findvalue('hit_reg');
			
			if (! $query_reg =~ /\d+/) { 
				next;
			}

			if (!$hit_reg =~ /\d+/) { 
				next;
			}

			$self_comps[$ci]{aligner}	= $aligner;
			
			$self_comps[$ci]{query_reg}	= $query_reg;
			$self_comps[$ci]{hit_reg}	= $hit_reg;
			$ci++;
		}
	}
	return \@self_comps;
}
sub define_coils { 

	my $sub = 'define_coils';
	my ($domain_xml_doc, $coil_aref, $unused_seq_aref, $used_seq_aref, $global_opt, $query_pdb, $query_chain, $query_pdbnum_aref) = @_;

	print "DEBUG: COILS top\n";
	my $coil_list_node = $domain_xml_doc->createElement('coil_list');
	my $domain_doc_node = $domain_xml_doc->findnodes('//chain_domains_set_top')->get_node(1);
	$domain_doc_node->appendChild($coil_list_node);
	my $cnum = 1;
	foreach my $coil (@$coil_aref) { 

		my $coil_seqid_str	= $$coil{seqid_range};
		my $coil_method		= $$coil{method};
		my $coil_seqid_aref 	= range_expand($coil_seqid_str);
		my $query_coverage 	= region_coverage($coil_seqid_aref, $unused_seq_aref);
		my $query_used_coverage	= region_coverage($coil_seqid_aref, $used_seq_aref);

		if ($query_coverage > $$global_opt{new_coverage_threshold} &&
			$query_used_coverage < $$global_opt{old_coverage_threshold}) { 
			my $coil_id		= "c$query_pdb$cnum";
			$cnum++;
			print "DEFINE coil $coil_seqid_str $coil_id\n";
			my $coil_node =	$domain_xml_doc->createElement('coil');

			$coil_node->setAttribute('method', $coil_method);
			$coil_node->setAttribute('coil_id', $coil_id);
			
			my $struct_seqid_range_node = $domain_xml_doc->createElement('struct_seqid_range');
			$struct_seqid_range_node->appendTextNode($coil_seqid_str);
			$coil_node->appendChild($struct_seqid_range_node);

			my $ungapped_seqid_range = ungap_range($coil_seqid_str, $$global_opt{gap_tol});
			my $ungapped_seqid_range_aref = range_expand($ungapped_seqid_range);
			my $ungapped_seqid_range_node = $domain_xml_doc->createElement('ungapped_seqid_range');
			$ungapped_seqid_range_node->setAttribute('gap_tolerance', $$global_opt{gap_tol});
			$ungapped_seqid_range_node->appendTextNode($ungapped_seqid_range);
			$coil_node->appendChild($ungapped_seqid_range_node);

			my $pdb_range	= pdb_rangify($coil_seqid_aref, $query_pdbnum_aref);
			my $pdb_range_node = $domain_xml_doc->createElement('pdb_range');
			$pdb_range_node->appendTextNode($pdb_range);
			$coil_node->appendChild($pdb_range_node);

			my $ungapped_pdb_range = pdb_rangify($ungapped_seqid_range_aref, $query_pdbnum_aref);
			my $ungapped_pdb_range_node = $domain_xml_doc->createElement('ungapped_pdb_range');
			$ungapped_pdb_range_node->setAttribute('gap_tolerance', $$global_opt{gap_tol});
			$ungapped_pdb_range_node->appendTextNode($ungapped_pdb_range);
			$coil_node->appendChild($ungapped_pdb_range_node);

			$coil_list_node->appendChild($coil_node);

			range_include($coil_seqid_aref, $used_seq_aref);
			range_exclude($coil_seqid_aref, $unused_seq_aref);

		}

	}

}
sub find_hhsearch_domains { 
	my $sub = 'find_hhsearch_domains';

	my (
		$blast_summ_xml_doc,
		$domain_xml_doc,
		$ref_chain_domains,
		$ref_domain_uid_lookup,
		$ref_range_cache,
		$query_struct_seqid_aref,
		$query_seqid_aref,
		$query_pdbnum_aref, 
		$self_comp_aref,
		$used_seq_aref,
		$unused_seq_aref,
		$input_mode,
		$global_opt,
		$reference,
		$coil_aref
		) = @_;

	my $domain_count = $domain_xml_doc->exists('//domain') 
		? $domain_xml_doc->findnodes('//domain')->size() + 1
		: 1;
	
	my $domain_list_node = $domain_xml_doc->findnodes('//domain_list')->get_node(1);

	my $query_pdb = $domain_xml_doc->findvalue('//@pdb_id');
	my $query_chain = $domain_xml_doc->findvalue('//@chain_id');

	my $hhsearch_run_XPath 	= q{//blast_summ_doc/blast_summ/hh_run[@program='hhsearch']};
	my $hhsearch_run_nodes 	= $blast_summ_xml_doc->findnodes($hhsearch_run_XPath);
	my $hhsearch_run_nodes_size = $hhsearch_run_nodes->size();
	if ($hhsearch_run_nodes_size > 1) { 
		print "ERROR! $sub: Odd number of hhsearch run nodies found {$hhsearch_run_nodes_size), aborting...\n";
	}

	if ($hhsearch_run_nodes_size == 1) { 
		my $hhsearch_hit_nodes = $hhsearch_run_nodes->get_node(1)->findnodes('hits/hit');
		foreach my $hit_node ($hhsearch_hit_nodes->get_nodelist() ) { 
			if ($hit_node->findvalue('@structure_obsolete') eq 'true' ) { 
				print "obsolete hit, skipping\n",
				next;
			}
			#if ($hit_node->findvalue('@low_coverage_hit') eq 'true') { 
			#	next;
			#}

			my $hit_num 	= $hit_node->findvalue('@num');
			my $hit_domain_id	= $hit_node->findvalue('@domain_id');
			if ($hit_node->findvalue('@structure_obsolete') eq 'true' ) { 
				print "obsolete hit $hit_num $hit_domain_id, skipping\n",
				next;
			}

			if ($hit_node->findvalue('@hit_cover') < $$global_opt{hit_coverage_threshold}) { #what
				print "Skipping $hit_num $hit_domain_id $$global_opt{hit_coverage_threshold}\n";
				next;
			}

			my $hh_prob	= $hit_node->findvalue('@hh_prob');
			my $hh_score	= $hit_node->findvalue('@hh_score');

			my $hit_query_reg 	= $hit_node->findvalue('query_reg');
			my $hit_hit_seqid 	= $hit_node->findvalue('hit_reg');

			#9/20/2012 why would this be translated from pos-inx to seqid? It's already in seqid?
			my $hit_query_struct_seqid_aref;
			if ($input_mode eq 'seqid') { 
				$hit_query_struct_seqid_aref = struct_region(range_expand($hit_query_reg), $query_seqid_aref);
			}elsif($input_mode eq 'struct_seqid') { 
				$hit_query_struct_seqid_aref = struct_region(range_expand($hit_query_reg), $query_struct_seqid_aref);
			}else{
				die "ERROR! Unknown input_mode $input_mode\n";
			}
			#my $hit_query_struct_seqid_aref = range_expand($hit_query_reg); #This is not structured residues only
			#my $hit_query_struct_seqid_aref = isect(range_expand($hit_query_reg), $query_struct_seqid_aref);
			my $query_struct_seqid = rangify(@$query_struct_seqid_aref);
			if (!$query_struct_seqid_aref || !$query_struct_seqid) { 
				printf "WARNING! No structured query range? $hit_domain_id $hit_query_reg $query_struct_seqid $$query_struct_seqid_aref[0] $$query_struct_seqid_aref[1] %i\n", scalar(@$query_struct_seqid_aref);
				next;
			}

			if (!$hit_query_struct_seqid_aref || ! rangify(@$hit_query_struct_seqid_aref)) { 
				print "WARNING! No range for hit $hit_num $hit_domain_id $query_pdb $query_chain\n";
				next;
			}
			
			my $hit_query_struct_seqid = rangify(@$hit_query_struct_seqid_aref);
			if ($DEBUG) { 
				print "DEBUG $sub: hhsearch $hit_num $hit_domain_id socore:$hh_score posx:$hit_query_reg seqid:$hit_query_struct_seqid\n";
			}

			my $unused_coverage = region_coverage($query_seqid_aref, $unused_seq_aref);

			if ($unused_coverage < 0.05 || scalar(@$unused_seq_aref) < 10) { 
				#print "query complete\n";
				last;
			}
			my $query_coverage = region_coverage($hit_query_struct_seqid_aref, $unused_seq_aref);
			my $query_used_coverage = residue_coverage($hit_query_struct_seqid_aref, $used_seq_aref);
			if ($DEBUG) { 
				my $range1 = rangify(@$hit_query_struct_seqid_aref);
				my $range2 = rangify(@$unused_seq_aref);
				my $query_residue_coverage = residue_coverage($hit_query_struct_seqid_aref, $unused_seq_aref);
				my $unused_rescount = scalar(@$unused_seq_aref);
				my $query_rescount = scalar(@$hit_query_struct_seqid_aref);
				if ($DEBUG) { 
					print "DEBUG $sub: unused_coverage: $unused_coverage query_coverage $query_coverage query_used_coverage $query_used_coverage\n";
					print "DEBUG $sub: r1 $range1 r2 $range2 qrc $query_residue_coverage urc $unused_rescount qrescount $query_rescount\n";
				}
			}

			if ($query_coverage > $$global_opt{new_coverage_threshold} && $query_used_coverage < $$global_opt{gap_tol} && $hh_prob > 90 && scalar(@$hit_query_struct_seqid_aref) > $$global_opt{gap_tol}) { 

				print "DEFINE hhsearch: $sub: $hit_domain_id $hit_query_struct_seqid\n";

				my $domain_node	= $domain_xml_doc->createElement('domain');
				my $domain_id	= "e". lc( $query_pdb) . "$query_chain$domain_count";
				$domain_node->setAttribute('domain_id', $domain_id);

				my $method = 'hh_full';
				my $method_node = $domain_xml_doc->createElement('method');
				$method_node->appendTextNode($method);
				$domain_node->appendChild($method_node);

				my $struct_seqid_range_node	= $domain_xml_doc->createElement('struct_seqid_range');
				$struct_seqid_range_node->appendTextNode($hit_query_struct_seqid);
				$domain_node->appendChild($struct_seqid_range_node);

				my $hit_struct_seqid_range_node	 = $domain_xml_doc->createElement('hit_struct_seqid_range');
				$hit_struct_seqid_range_node->appendTextNode($hit_hit_seqid);
				$domain_node->appendChild($hit_struct_seqid_range_node);

				my $ungapped_seqid_range = ungap_range($hit_query_struct_seqid, $$global_opt{gap_tol});
				my $ungapped_seqid_range_aref = range_expand($ungapped_seqid_range);
				my $ungapped_seqid_range_node	= $domain_xml_doc->createElement('ungapped_seqid_range');
				$ungapped_seqid_range_node->setAttribute('gap_tolerance', $$global_opt{gap_tol});
				$ungapped_seqid_range_node->appendTextNode($ungapped_seqid_range);
				$domain_node->appendChild($ungapped_seqid_range_node);

				my $pdb_range = pdb_rangify($hit_query_struct_seqid_aref, $query_pdbnum_aref);
				my $pdb_range_node	= $domain_xml_doc->createElement('pdb_range');
				$pdb_range_node->appendTextNode($pdb_range);
				$domain_node->appendChild($pdb_range_node);

				my $ungapped_pdb_range = pdb_rangify($ungapped_seqid_range_aref, $query_pdbnum_aref);
				my $ungapped_pdb_range_node	= $domain_xml_doc->createElement('ungapped_pdb_range');
				$ungapped_pdb_range_node->setAttribute('gap_tolerance', $$global_opt{gap_tol});
				$ungapped_pdb_range_node->appendTextNode($ungapped_pdb_range);
				$domain_node->appendChild($ungapped_pdb_range_node);

				my $hit_domain_id_node	= $domain_xml_doc->createElement('hit_domain');
				$hit_domain_id_node->setAttribute('ecod_domain_id', $hit_domain_id);
				$hit_domain_id_node->setAttribute('reference', $reference);
				$hit_domain_id_node->setAttribute('uid', $$ref_domain_uid_lookup{$hit_domain_id});

				$domain_node->appendChild($hit_domain_id_node);

				my $hh_score_node	= $domain_xml_doc->createElement('hh_score');
				$hh_score_node->setAttribute('hh_prob', $hh_prob);
				$hh_score_node->setAttribute('hh_score', $hh_score);
				$domain_node->appendChild($hh_score_node);

				$domain_list_node->appendChild($domain_node);

				range_exclude($ungapped_seqid_range_aref, $unused_seq_aref);
				range_include($ungapped_seqid_range_aref, $used_seq_aref);

				$domain_count++;
				#Self_comparison check
				for (my $i = 0; $i < scalar(@$self_comp_aref); $i++) { 

					my $self_comp_aligner 	= $$self_comp_aref[$i]{aligner};
					my $self_comp_query_reg	= $$self_comp_aref[$i]{query_reg};
					my $self_comp_hit_reg	= $$self_comp_aref[$i]{hit_reg};
					
					my $sc_query_reg_aref 	= range_expand($self_comp_query_reg);
					my $sc_hit_reg_aref	= range_expand($self_comp_hit_reg);
					
					my $sc_query_coverage_1 = region_coverage($sc_query_reg_aref, $ungapped_seqid_range_aref);
					my $sc_query_coverage_2 = region_coverage($ungapped_seqid_range_aref, $sc_query_reg_aref,);
					my $sc_hit_coverage_1 = region_coverage($sc_hit_reg_aref, $unused_seq_aref);
					#my $sc_hit_coverage_2 = region_coverage(\@unused_seq, $sc_hit_reg_aref);

					if ($DEBUG) { 
						print "DEBUG hhsearch SC $i $self_comp_aligner qr $self_comp_query_reg hr $self_comp_hit_reg sqc1 $sc_query_coverage_1 sqc2 $sc_query_coverage_2 shc1 $sc_hit_coverage_1\n";
					}

					if ($sc_query_coverage_1 > 0.9 && $sc_query_coverage_2 > 0.9 && $sc_hit_coverage_1 > 0.9 ) { 
						print "DEFINE hhsearch self_comp $sub: $self_comp_aligner\n";

						my $sc_domain_node 	= $domain_xml_doc->createElement('domain');

						my $domain_id	= 'e' . lc($query_pdb) . "$query_chain$domain_count";
						$sc_domain_node->setAttribute('domain_id', $domain_id);

						my $method_node	= $domain_xml_doc->createElement('method');
						my $method 	= 'hhsearch_self_comp';
						$method_node->appendTextNode($method);
						$method_node->setAttribute('self_comp_aligner', $self_comp_aligner);
						$sc_domain_node->appendChild($method_node);

						my $struct_seqid_range_node	= $domain_xml_doc->createElement('struct_seqid_range');
						$struct_seqid_range_node->appendTextNode($self_comp_hit_reg);
						$sc_domain_node->appendChild($struct_seqid_range_node);

						my $ungapped_seqid_range	= ungap_range($self_comp_hit_reg, $$global_opt{gap_tol});
						my $ungapped_seqid_range_aref	= range_expand($ungapped_seqid_range);
						my $ungapped_seqid_range_node	= $domain_xml_doc->createElement('ungapped_seqid_range');
						$ungapped_seqid_range_node->setAttribute('gap_tolerance', $$global_opt{gap_tol});
						$ungapped_seqid_range_node->appendTextNode($ungapped_seqid_range);
						$sc_domain_node->appendChild($ungapped_seqid_range_node);

						my $hit_struct_seqid_range_node = $domain_xml_doc->createElement('hit_struct_seqid_range');
						$hit_struct_seqid_range_node->appendTextNode($self_comp_query_reg);
						$sc_domain_node->appendChild($hit_struct_seqid_range_node);
						
						my $pdb_range	= pdb_rangify($sc_hit_reg_aref, $query_pdbnum_aref);
						my $pdb_range_node	= $domain_xml_doc->createElement('pdb_range');
						$pdb_range_node->appendTextNode($pdb_range);
						$sc_domain_node->appendChild($pdb_range_node);

						my $hit_domain_id_node	= $domain_xml_doc->createElement('hit_domain');
						$hit_domain_id_node->setAttribute('ecod_domain_id', $hit_domain_id);
						$hit_domain_id_node->setAttribute('reference', $reference);
						$hit_domain_id_node->setAttribute('uid', $$ref_domain_uid_lookup{$hit_domain_id});
						$sc_domain_node->appendChild($hit_domain_id_node);

						my $ungapped_pdb_range	= pdb_rangify($ungapped_seqid_range_aref, $query_pdbnum_aref);
						my $ungapped_pdb_range_node	= $domain_xml_doc->createElement('ungapped_pdb_range');
						$ungapped_pdb_range_node->setAttribute('gap_tolerance', $$global_opt{gap_tol});
						$ungapped_pdb_range_node->appendTextNode($ungapped_pdb_range);
						$sc_domain_node->appendChild($ungapped_pdb_range_node);


						if ($self_comp_aligner eq 'dali') { 
							my $dali_score_node	= $domain_xml_doc->createElement('dali_zscore');
							$dali_score_node->appendTextNode($$self_comp_aref[$i]{zscore});
							$sc_domain_node->appendChild($dali_score_node);
						}elsif($self_comp_aligner eq 'hhrepid') { 
							my $hh_prob_node	= $domain_xml_doc->createElement('hhrepid_prob');
							$hh_prob_node->appendTextNode($$self_comp_aref[$i]{prob});
							$sc_domain_node->appendChild($hh_prob_node);
						}else{
							print "WARNING! self comp aligner $self_comp_aligner unknown\n";
						}

						$domain_list_node->appendChild($sc_domain_node);
						range_include($ungapped_seqid_range_aref, $used_seq_aref);
						range_exclude($ungapped_seqid_range_aref, $unused_seq_aref);
						

						$domain_count++;
					}

				}
			}else{
				if ($DEBUG) { 
					print "query_coverage $query_coverage query_used_coverage $query_used_coverage hh_prob $hh_prob scalar " .  scalar(@$hit_query_struct_seqid_aref). "\n";
				}

			}
		}
	}
}
sub find_hh_dali_domains { 
	my $sub = 'find_hh_dali_domains';
	
	my ($blast_summ_xml_doc,
		$dali_summ_xml_doc,
		$domain_xml_doc,
		$ref_chain_domains,
		$ref_range_cache,
		$ref_domain_uid_lookup,
		$query_struct_seqid_aref,
		$query_seqid_aref,
		$query_pdbnum_aref,
		$self_comp_aref,
		$used_seq_aref,
		$unused_seq_aref,
		$input_mode,
		$global_opt,
		$reference,
		$coil_aref) = @_;

	my $domain_count = $domain_xml_doc->exists('//domain')
		? $domain_xml_doc->findnodes('//domain')->size() + 1
		: 1;

	my $domain_list_node = $domain_xml_doc->findnodes('//domain_list')->get_node(1);		

	my $query_pdb 	= $domain_xml_doc->findvalue('//@pdb_id');
	my $query_chain = $domain_xml_doc->findvalue('//@chain_id');

	my $dali_i = 0;
	my %dali_hits;
	foreach my $dali_hit_node (find_dali_hit_nodes($dali_summ_xml_doc)) { 
		my ($uid, $hit_ecod_domain_id) 	= get_ids($dali_hit_node);
		my ($z, $rmsd, $id) 			= get_dali_scores($dali_hit_node);	
		my ($hit_reg, $query_reg, $coverage) = get_dali_regions($dali_hit_node);

		my %dali_hit;
		$dali_hit{ecod_domain_id}	= $hit_ecod_domain_id;
		$dali_hit{query_range}		= $query_reg;
		$dali_hit{hit_range}		= $hit_reg;
		$dali_hit{zscore}			= $z;
		push (@{$dali_hits{$hit_ecod_domain_id}}, \%dali_hit);
		$dali_i++;

	}

	my $hh_i = 0;
	my %hh_hits;
	foreach my $hh_hit_node (find_hhsearch_hit_nodes($blast_summ_xml_doc) ) { 
		my $hit_num = $hh_hit_node->findvalue('@num');

		my $hit_ecod_domain_id = $hh_hit_node->findvalue('@domain_id');

		my $hh_prob		= $hh_hit_node->findvalue('@hh_prob');
		my $query_reg	= $hh_hit_node->findvalue('query_reg');
		my $hit_reg 	= $hh_hit_node->findvalue('hit_reg');
		my $hit_query_struct_seqid_aref; # i.e. Query residues covered by hit
		if ($input_mode  eq 'seqid') { 
			$hit_query_struct_seqid_aref = struct_region(range_expand($query_reg), $query_seqid_aref);
		}elsif ($input_mode  eq 'struct_seqid') { 
			$hit_query_struct_seqid_aref = struct_region(range_expand($query_reg), $query_struct_seqid_aref);
		}
		next unless ref $hit_query_struct_seqid_aref eq 'ARRAY';
		my $hit_query_struct_seqid = rangify(@$hit_query_struct_seqid_aref);

		my %hh_hit;
		$hh_hit{ecod_domain_id} = $hit_ecod_domain_id;
		$hh_hit{query_range}	= $hit_query_struct_seqid;
		$hh_hit{hit_range}		= $hit_reg;
		$hh_hit{hh_prob}		= $hh_prob;
		push (@{$hh_hits{$hit_ecod_domain_id}}, \%hh_hit);
		$hh_i++;
	}

	my %mixed;
	foreach my $key (keys %hh_hits, keys %dali_hits) { 
		$mixed{$key}++;
	}

	foreach my $key (keys %mixed) { 
		my $hit_ecod_domain_id = $key;
		foreach my $dali_align (@{$dali_hits{$key}}) { 
			my $dali_range_aref = range_expand($dali_align->{query_range});
			foreach my $hh_align (@{$hh_hits{$key}}) { 
				if ($hh_align->{query_range} eq 'NA') { next } 
				my $hh_range_aref = range_expand($hh_align->{query_range});
				if (bidirectional_coverage(range_expand($dali_align->{query_range}), range_expand($hh_align->{query_range}), 0.5)) { 
					my $mixed_range_aref = union($dali_range_aref, $hh_range_aref); #Not isect?
					my $mixed_range = rangify(@$mixed_range_aref);

					my $unused_coverage = region_coverage($query_seqid_aref, $unused_seq_aref);

					if ($unused_coverage < 0.05 || scalar(@$unused_seq_aref) < 10) { 
						print "query complete\n";
						last;
					}
					my $query_coverage		 	= region_coverage($mixed_range_aref, $unused_seq_aref);
					my $query_used_coverage		= residue_coverage($mixed_range_aref, $used_seq_aref);
					#print "qc:$query_coverage quc: $query_used_coverage nct:$$global_opt{new_coverage_threshold}\n";


					if ($query_coverage > $$global_opt{new_coverage_threshold} && 
						$query_used_coverage < $$global_opt{gap_tol} && 
						$dali_align->{zscore} > 4 && 
						$hh_align->{hh_prob} > 50 &&
						scalar(@$mixed_range_aref) > $$global_opt{gap_tol}) { 
						
						print "DEFINE hh_dali: $sub: $hit_ecod_domain_id\n";

						my $domain_node	= $domain_xml_doc->createElement('domain');
						my $domain_id = 'e' . lc($query_pdb) . $query_chain . $domain_count;
						$domain_node->setAttribute('domain_id', $domain_id);

						$domain_node->appendTextChild('method', 'hh_dali');

						$domain_node->appendTextChild('struct_seqid_range', $mixed_range);

						$domain_node->appendTextChild('hit_dali_range', $dali_align->{hit_range});
						$domain_node->appendTextChild('hit_hh_range', 	$hh_align->{hit_range});

						my $ungapped_seqid_range = ungap_range($mixed_range, $$global_opt{gap_tol});
						my $ungapped_seqid_range_aref = range_expand($ungapped_seqid_range);
						my $usr_node = $domain_xml_doc->createElement('ungapped_seqid_range');
						$usr_node->appendTextNode($ungapped_seqid_range);
						$usr_node->setAttribute('gap_tolerance', $$global_opt{gap_tol});
						$domain_node->appendChild($usr_node);

						$domain_node->appendTextChild('pdb_range', pdb_rangify($mixed_range_aref, $query_pdbnum_aref));
						
						my $ungapped_pdb_range = pdb_rangify($ungapped_seqid_range_aref, $query_pdbnum_aref);
						my $ugp_node = $domain_xml_doc->createElement('ungapped_pdb_range');
						$ugp_node->appendTextNode($ungapped_pdb_range);
						$ugp_node->setAttribute('gap_tolerance', $$global_opt{gap_tol});
						$domain_node->appendChild($ugp_node);

						my $hdi_node = $domain_xml_doc->createElement('hit_domain');
						$hdi_node->setAttribute('ecod_domain_id', $hit_ecod_domain_id);
						$hdi_node->setAttribute('reference', $reference);
						$hdi_node->setAttribute('uid', $$ref_domain_uid_lookup{$hit_ecod_domain_id});
						$domain_node->appendChild($hdi_node);

						my $hh_score_node		= $domain_xml_doc->createElement('hh_score');
						$hh_score_node->setAttribute('hh_prob', $hh_align->{hh_prob});
						$domain_node->appendChild($hh_score_node);

						my $dali_score_node 	= $domain_xml_doc->createElement('dali_score');
						$dali_score_node->setAttribute('z-score', $dali_align->{zscore});
						$domain_node->appendChild($dali_score_node);

						$domain_list_node->appendChild($domain_node);

						range_exclude($ungapped_seqid_range_aref, $unused_seq_aref);
						range_include($ungapped_seqid_range_aref, $used_seq_aref);

						$domain_count++;
	
					}
				}
			}
		}
	}
}
sub find_psiblast_domains { 
	my $sub = 'find_psiblast_domains';

	my (
		$blast_summ_xml_doc, 
		$domain_xml_doc, 
		$ref_chain_domains, 
		$ref_range_cache,
		$ref_domain_uid_lookup,
		$query_struct_seqid_aref,
		$query_seqid_aref, 
		$query_pdbnum_aref,
		$self_comp_aref,
		$used_seq_aref,
		$unused_seq_aref,
		$input_mode,
		$global_opt,
		$reference,
		$coil_aref) = @_;

	my $domain_count = $domain_xml_doc->exists('//domain') 
		? $domain_xml_doc->findnodes('//domain')->size() + 1
		: 1;
	
	my $domain_list_node = $domain_xml_doc->findnodes('//domain_list')->get_node(1);


	my $query_pdb 	= $domain_xml_doc->findvalue('//@pdb_id');
	my $query_chain = $domain_xml_doc->findvalue('//@chain_id');

	my $psiblast_run_XPath = q{//blast_summ_doc/blast_summ/blast_run[@program='psiblast']};
	my $psiblast_run_nodes = $blast_summ_xml_doc->findnodes($psiblast_run_XPath);
	my $psiblast_run_nodes_size = $psiblast_run_nodes->size();
	if ($psiblast_run_nodes_size > 1){ 
		print "ERROR! $sub: Odd number of psiblast run nodes found ($psiblast_run_nodes_size), aborting...\n";
	}

	if ($psiblast_run_nodes_size == 1 && !$$global_opt{no_psiblast}) { 
		my $psiblast_hit_nodes = $psiblast_run_nodes->get_node(1)->findnodes('hits/hit');
		foreach my $hit_node ($psiblast_hit_nodes->get_nodelist() ) { 

			my $hit_num	= $hit_node->findvalue('@num');
			my $hit_domain_id = $hit_node->findvalue('@domain_id');
			my $hsp_count 	= $hit_node->findvalue('@hsp_count');
			my $evalues 	= $hit_node->findvalue('@evalues');

			my $evalue = 'Unk';
			if ($hsp_count == 1) { 
				$evalue = $evalues;
			}

			#This SHOUUULD be seqid_reange
			#my $hit_query_seqid_range = $hit_node->findvalue('text()');
			#my $hit_query_reg = $hit_node->findvalue('text()');
			my $hit_query_reg = $hit_node->findvalue('query_reg');
			my $hit_query_struct_seqid_aref;
			if ($input_mode eq 'struct_seqid') { 
				$hit_query_struct_seqid_aref = struct_region(range_expand($hit_query_reg), $query_struct_seqid_aref);
			}elsif ($input_mode eq 'seqid') { 
				$hit_query_struct_seqid_aref = range_expand($hit_query_reg);
			}else{
				die "ERROR! input_mode $input_mode not recognized\n";
			}
			my $hit_query_struct_seqid_range = rangify(@$hit_query_struct_seqid_aref);

			my $hit_hit_reg = $hit_node->findvalue('hit_reg');


			if ($DEBUG) { 
				print "DEBUG $sub: psiblast $hit_num $hit_domain_id $hit_query_struct_seqid_range\n";
			}
			#Assign if overlaps 

			my $unused_coverage = region_coverage($query_seqid_aref, $unused_seq_aref);
			
			if ($unused_coverage < 0.05 || scalar(@$unused_seq_aref) < 10) { 
				last;
			}
			my $query_coverage = region_coverage(range_expand($hit_query_struct_seqid_range), $unused_seq_aref);


			if ($DEBUG) { 
				print "DEBUG $sub: uc $unused_coverage qc $query_coverage\n";
			}

			if ($query_coverage > $$global_opt{new_coverage_threshold}) { 
				print "DEFINE psiblast: $sub: $hit_domain_id $hit_query_struct_seqid_range\n";

				my $domain_node		= $domain_xml_doc->createElement('domain');
				#my $domain_id		= "d$query_pdb$query_chain$domain_count";
				my $domain_id	= "e". lc( $query_pdb) . "$query_chain$domain_count";
				$domain_node->setAttribute('domain_id', $domain_id);

				my $method_node	= $domain_xml_doc->createElement('method');
				my $method	= 'psiblast';
				$method_node->appendTextNode($method);
				$domain_node->appendChild($method_node);

				my $struct_seqid_range_node	= $domain_xml_doc->createElement('struct_seqid_range');
				$struct_seqid_range_node->appendTextNode($hit_query_struct_seqid_range);	
				$domain_node->appendChild($struct_seqid_range_node);

				my $hit_struct_seqid_range_node	= $domain_xml_doc->createElement('hit_struct_seqid_range');
				$hit_struct_seqid_range_node->appendTextNode($hit_hit_reg);
				$domain_node->appendChild($hit_struct_seqid_range_node);

				my $ungapped_seqid_range = ungap_range($hit_query_struct_seqid_range, $$global_opt{gap_tol});
				my $ungapped_seqid_range_aref = range_expand($ungapped_seqid_range);
				my $ungapped_seqid_range_node	= $domain_xml_doc->createElement('ungapped_seqid_range');
				$ungapped_seqid_range_node->setAttribute('gap_tolerance', $$global_opt{gap_tol});
				$ungapped_seqid_range_node->appendTextNode($ungapped_seqid_range);
				$domain_node->appendChild($ungapped_seqid_range_node);

				my $pdb_range = pdb_rangify($hit_query_struct_seqid_aref, $query_pdbnum_aref);
				my $pdb_range_node	= $domain_xml_doc->createElement('pdb_range');
				$pdb_range_node->appendTextNode($pdb_range);
				$domain_node->appendChild($pdb_range_node);

				my $ungapped_pdb_range = pdb_rangify($ungapped_seqid_range_aref, $query_pdbnum_aref);
				my $ungapped_pdb_range_node	= $domain_xml_doc->createElement('ungapped_pdb_range');
				$ungapped_pdb_range_node->setAttribute('gap_tolerance', $$global_opt{gap_tol});
				$ungapped_pdb_range_node->appendTextNode($ungapped_pdb_range);
				$domain_node->appendChild($ungapped_pdb_range_node);

				my $hit_domain_id_node	= $domain_xml_doc->createElement('hit_domain');
				$hit_domain_id_node->setAttribute('ecod_domain_id', $hit_domain_id);
				$hit_domain_id_node->setAttribute('reference', $reference);
				$hit_domain_id_node->setAttribute('uid', $$ref_domain_uid_lookup{$hit_domain_id});
				$domain_node->appendChild($hit_domain_id_node);

				my $blast_score_node	= $domain_xml_doc->createElement('blast_eval');
				$blast_score_node->appendTextNode($evalue);
				$domain_node->appendChild($blast_score_node);

				$domain_list_node->appendChild($domain_node);
			
				range_exclude(range_expand($hit_query_struct_seqid_range), $unused_seq_aref);
				range_include(range_expand($hit_query_struct_seqid_range), $used_seq_aref);

				$domain_count++;
				#Self_comparison check
				for (my $i = 0; $i < scalar(@$self_comp_aref); $i++) { 

					my $self_comp_aligner 	= $$self_comp_aref[$i]{aligner};
					my $self_comp_query_reg	= $$self_comp_aref[$i]{query_reg};
					my $self_comp_hit_reg	= $$self_comp_aref[$i]{hit_reg};
					
					my $sc_query_reg_aref 	= range_expand($self_comp_query_reg);
					my $sc_hit_reg_aref	= range_expand($self_comp_hit_reg);
					
					my $sc_query_coverage_1 = region_coverage($sc_query_reg_aref, $ungapped_seqid_range_aref);
					my $sc_query_coverage_2 = region_coverage($ungapped_seqid_range_aref, $sc_query_reg_aref,);
					my $sc_hit_coverage_1 = region_coverage($sc_hit_reg_aref, $unused_seq_aref);
					#my $sc_hit_coverage_2 = region_coverage(\@unused_seq, $sc_hit_reg_aref);

					if ($DEBUG) { 
						print "DEBUG psiblast SC $i $self_comp_aligner qr $self_comp_query_reg hr $self_comp_hit_reg sqc1 $sc_query_coverage_1 sqc2 $sc_query_coverage_2 shc1 $sc_hit_coverage_1\n";
					}

					if ($sc_query_coverage_1 > 0.9 && $sc_query_coverage_2 > 0.9 && $sc_hit_coverage_1 > 0.9 ) { 
						print "DEFINE psiblast self_comp $sub: $self_comp_aligner\n";

						my $sc_domain_node 	= $domain_xml_doc->createElement('domain');

						my $domain_id	= 'e' . lc($query_pdb) . "$query_chain$domain_count";
						$sc_domain_node->setAttribute('domain_id', $domain_id);

						my $method_node	= $domain_xml_doc->createElement('method');
						my $method 	= 'psiblast_self_comp';
						$method_node->appendTextNode($method);
						$method_node->setAttribute('self_comp_aligner', $self_comp_aligner);
						$sc_domain_node->appendChild($method_node);

						my $struct_seqid_range_node	= $domain_xml_doc->createElement('struct_seqid_range');
						$struct_seqid_range_node->appendTextNode($self_comp_hit_reg);
						$sc_domain_node->appendChild($struct_seqid_range_node);

						my $ungapped_seqid_range	= ungap_range($self_comp_hit_reg, $$global_opt{gap_tol});
						my $ungapped_seqid_range_aref	= range_expand($ungapped_seqid_range);
						my $ungapped_seqid_range_node	= $domain_xml_doc->createElement('ungapped_seqid_range');
						$ungapped_seqid_range_node->setAttribute('gap_tolerance', $$global_opt{gap_tol});
						$ungapped_seqid_range_node->appendTextNode($ungapped_seqid_range);
						$sc_domain_node->appendChild($ungapped_seqid_range_node);

						my $hit_struct_seqid_range_node = $domain_xml_doc->createElement('hit_struct_seqid_range');
						$hit_struct_seqid_range_node->appendTextNode($self_comp_query_reg);
						$sc_domain_node->appendChild($hit_struct_seqid_range_node);
						
						my $pdb_range	= pdb_rangify($sc_hit_reg_aref, $query_pdbnum_aref);
						my $pdb_range_node	= $domain_xml_doc->createElement('pdb_range');
						$pdb_range_node->appendTextNode($pdb_range);
						$sc_domain_node->appendChild($pdb_range_node);

						my $ungapped_pdb_range	= pdb_rangify($ungapped_seqid_range_aref, $query_pdbnum_aref);
						my $ungapped_pdb_range_node	= $domain_xml_doc->createElement('ungapped_pdb_range');
						$ungapped_pdb_range_node->setAttribute('gap_tolerance', $$global_opt{gap_tol});
						$ungapped_pdb_range_node->appendTextNode($ungapped_pdb_range);
						$sc_domain_node->appendChild($ungapped_pdb_range_node);


						if ($self_comp_aligner eq 'dali') { 
							my $dali_score_node	= $domain_xml_doc->createElement('dali_zscore');
							$dali_score_node->appendTextNode($$self_comp_aref[$i]{zscore});
							$sc_domain_node->appendChild($dali_score_node);
						}elsif($self_comp_aligner eq 'hhrepid') { 
							my $hh_prob_node	= $domain_xml_doc->createElement('hhrepid_prob');
							$hh_prob_node->appendTextNode($$self_comp_aref[$i]{prob});
							$sc_domain_node->appendChild($hh_prob_node);
						}else{
							print "WARNING! self comp aligner $self_comp_aligner unknown\n";
						}
						$domain_list_node->appendChild($sc_domain_node);
						range_include($ungapped_seqid_range_aref, $used_seq_aref);
						range_exclude($ungapped_seqid_range_aref, $unused_seq_aref);
						

						$domain_count++;
					}
				}
			}
		}	
	}
}

sub find_domainwise_blast_domains { 
	my $sub = 'find_domainwise_blast_domains';

	my (
		$blast_summ_xml_doc,
		$domain_xml_doc,
		$ref_chain_domains,
		$ref_domain_uid_lookup,
		$ref_range_cache,
		$query_struct_seqid_aref,
		$query_seqid_aref,
		$query_pdbnum_aref,
		$self_comp_aref,
		$used_seq_aref,
		$unused_seq_aref,
		$input_mode,
		$global_opt,
		$reference,
		$coil_aref
		) = @_;

	my $domain_count = $domain_xml_doc->exists('//domain') 
		? $domain_xml_doc->findnodes('//domain')->size() + 1
		: 1;
	
	my $domain_list_node = $domain_xml_doc->findnodes('//domain_list')->get_node(1);

	my $query_pdb 	= $domain_xml_doc->findvalue('//@pdb_id');
	my $query_chain = $domain_xml_doc->findvalue('//@chain_id');

	my $blastp_run_XPath = q{//blast_summ_doc/blast_summ/blast_run[@program='blastp']};
	my $blastp_run_nodes = $blast_summ_xml_doc->findnodes($blastp_run_XPath);
	my $blastp_run_nodes_size = $blastp_run_nodes->size();
	if ($blastp_run_nodes_size != 1) { 
		print "ERROR! $sub: Odd number of blastp run nodes found ($blastp_run_nodes_size), aborting...\n";

	}
	my $blastp_hit_nodes = $blastp_run_nodes->get_node(1)->findnodes('hits/hit');
	foreach my $hit_node ($blastp_hit_nodes->get_nodelist() ) { 

		my $hit_num 	= $hit_node->findvalue('@num');
		my $hit_domain_id = $hit_node->findvalue('@domain_id');
		my $hsp_count 	= $hit_node->findvalue('@hsp_count');
		my $evalues 	= $hit_node->findvalue('@evalues');

		$hit_domain_id =~ /(e|d)(\w{4})(.)/;
		my $dtype = $1;
		my $hit_pdb = $2;
		my $hit_chain = $3;
		if ($DEBUG) { 
			print "hp $hit_pdb hc: $hit_chain\n";
		}
		if ($dtype eq 'd') {
			$hit_chain = uc($hit_chain);
		}

		my $evalue = 'Unk';
		if ($hsp_count == 1) { 
			$evalue = $evalues;
		}

		my $hit_query_reg	= $hit_node->findvalue('query_reg');
		if (!$query_struct_seqid_aref) { die "bad qssa\n";}
		my $hit_query_struct_seqid_aref;
		if ($input_mode eq 'struct_seqid') { 
			$hit_query_struct_seqid_aref = struct_region(range_expand($hit_query_reg), $query_struct_seqid_aref);
		}elsif($input_mode eq 'seqid') { 
			$hit_query_struct_seqid_aref = range_expand($hit_query_reg);
		}else{
			die "ERROR! input_mode $input_mode not recognized\n";
		}
		my $hit_query_struct_seqid_range = rangify(@$hit_query_struct_seqid_aref);
		
		my $hit_hit_reg		= $hit_node->findvalue('hit_reg');
		

		if ($DEBUG) { 
			print "DEBUG $sub: blastp $hit_num $hit_domain_id $hit_query_struct_seqid_range\n";
		}
		#Assign if overlaps 

		my $unused_coverage = region_coverage($query_struct_seqid_aref, $unused_seq_aref);
		my $used_coverage = region_coverage($query_struct_seqid_aref, $used_seq_aref);
		
		if ($unused_coverage < 0.05 || scalar(@$unused_seq_aref) < 10) { 
			last;
		}
		#my $query_coverage = region_coverage(range_expand($hit_query_seqid_range), \@unused_seq);
		my $query_new_coverage = region_coverage($hit_query_struct_seqid_aref, $unused_seq_aref);
		my $query_used_coverage = region_coverage($hit_query_struct_seqid_aref, $used_seq_aref);

		if ($DEBUG) { 
			print "DEBUG $sub: uc $unused_coverage qnc $query_new_coverage qoc $query_used_coverage\n";
		}

		if ($query_new_coverage > $$global_opt{new_coverage_threshold} && $query_used_coverage < $$global_opt{old_coverage_threshold}) { 
			print "DEFINE blastp: $sub: $hit_domain_id $hit_query_struct_seqid_range\n";

			my ($hit_seqid_aref, $hit_struct_seqid_aref, $hit_pdbnum_aref, $hit_asym_id) = pdbml_seq_parse($hit_pdb, $hit_chain);
			if (!$hit_seqid_aref || !$hit_struct_seqid_aref) { 
				print "WARNING! $sub: hit parse failure\n";
				next;
			}

			if (!$hit_hit_reg) { die "ERROR $sub: hit hit reg failure\n"}
			my $hit_reg_struct_seqid_aref = struct_region(range_expand($hit_hit_reg), $hit_struct_seqid_aref);
			my $domain_node		= $domain_xml_doc->createElement('domain');
			#my $domain_id		= "d$query_pdb$query_chain$domain_count";
			my $domain_id	= "e". lc( $query_pdb) . "$query_chain$domain_count";
			$domain_node->setAttribute('domain_id', $domain_id);

			my $method_node	= $domain_xml_doc->createElement('method');
			my $method	= 'blastp';
			$method_node->appendTextNode($method);
			$domain_node->appendChild($method_node);

			my $struct_seqid_range_node	= $domain_xml_doc->createElement('struct_seqid_range');
			$struct_seqid_range_node->appendTextNode($hit_query_struct_seqid_range);	
			$domain_node->appendChild($struct_seqid_range_node);

			my $ungapped_seqid_range = ungap_range($hit_query_struct_seqid_range, $$global_opt{gap_tol});
			my $ungapped_seqid_range_aref = range_expand($ungapped_seqid_range);
			my $ungapped_seqid_range_node	= $domain_xml_doc->createElement('ungapped_seqid_range');
			$ungapped_seqid_range_node->setAttribute('gap_tolerance', $$global_opt{gap_tol});
			$ungapped_seqid_range_node->appendTextNode($ungapped_seqid_range);
			$domain_node->appendChild($ungapped_seqid_range_node);

			my $hit_struct_seqid_range_node	 = $domain_xml_doc->createElement('hit_struct_seqid_range');
			$hit_struct_seqid_range_node->appendTextNode(rangify(@$hit_reg_struct_seqid_aref));
			$domain_node->appendChild($hit_struct_seqid_range_node);

			my $pdb_range = pdb_rangify($hit_query_struct_seqid_aref, $query_pdbnum_aref);
			my $pdb_range_node	= $domain_xml_doc->createElement('pdb_range');
			$pdb_range_node->appendTextNode($pdb_range);
			$domain_node->appendChild($pdb_range_node);

			my $ungapped_pdb_range = pdb_rangify($ungapped_seqid_range_aref, $query_pdbnum_aref);
			my $ungapped_pdb_range_node	= $domain_xml_doc->createElement('ungapped_pdb_range');
			$ungapped_pdb_range_node->setAttribute('gap_tolerance', $$global_opt{gap_tol});
			$ungapped_pdb_range_node->appendTextNode($ungapped_pdb_range);
			$domain_node->appendChild($ungapped_pdb_range_node);

			my $hit_domain_id_node	= $domain_xml_doc->createElement('hit_domain');
			$hit_domain_id_node->setAttribute('ecod_domain_id', $hit_domain_id);
			$hit_domain_id_node->setAttribute('reference', $reference);
			$hit_domain_id_node->setAttribute('uid', $$ref_domain_uid_lookup{$hit_domain_id});
			$domain_node->appendChild($hit_domain_id_node);

			my $blast_score_node	= $domain_xml_doc->createElement('blast_eval');
			$blast_score_node->appendTextNode($evalue);
			$domain_node->appendChild($blast_score_node);

			$domain_list_node->appendChild($domain_node);
			
			range_exclude(range_expand($hit_query_struct_seqid_range), $unused_seq_aref);
			range_include(range_expand($hit_query_struct_seqid_range), $used_seq_aref);

			$domain_count++;
			#Self_comparison check
			for (my $i = 0; $i < scalar(@$self_comp_aref); $i++) { 

				my $self_comp_aligner 	= $$self_comp_aref[$i]{aligner};
				my $self_comp_query_reg	= $$self_comp_aref[$i]{query_reg};
				my $self_comp_hit_reg	= $$self_comp_aref[$i]{hit_reg};
				
				my $sc_query_reg_aref 	= range_expand($self_comp_query_reg);
				my $sc_hit_reg_aref	= range_expand($self_comp_hit_reg);
				
				my $sc_query_coverage_1 = region_coverage($sc_query_reg_aref, $ungapped_seqid_range_aref);
				my $sc_query_coverage_2 = region_coverage($ungapped_seqid_range_aref, $sc_query_reg_aref,);
				my $sc_hit_coverage_1 = region_coverage($sc_hit_reg_aref, $unused_seq_aref);
				
				if ($DEBUG) { 
					print "DEBUG blastp SC $i $self_comp_aligner qr $self_comp_query_reg hr $self_comp_hit_reg sqc1 $sc_query_coverage_1 sqc2 $sc_query_coverage_2 shc1 $sc_hit_coverage_1\n";
				}

				if ($sc_query_coverage_1 > 0.9 && $sc_query_coverage_2 > 0.9 && $sc_hit_coverage_1 > 0.9 ) { 
					my $pre_used = rangify(@$used_seq_aref);
					my $pre_unused = rangify(@$unused_seq_aref);
					print "DEFINE blastp self_comp $sub: $self_comp_aligner qr $self_comp_query_reg hr $self_comp_hit_reg pu $pre_used pun $pre_unused\n";

					my $sc_domain_node 	= $domain_xml_doc->createElement('domain');

					my $domain_id	= 'e' . lc($query_pdb) . "$query_chain$domain_count";
					$sc_domain_node->setAttribute('domain_id', $domain_id);

					my $method_node	= $domain_xml_doc->createElement('method');
					my $method 	= 'blastp_self_comp';
					$method_node->appendTextNode($method);
					$method_node->setAttribute('self_comp_aligner', $self_comp_aligner);
					$sc_domain_node->appendChild($method_node);

					my $struct_seqid_range_node	= $domain_xml_doc->createElement('struct_seqid_range');
					$struct_seqid_range_node->appendTextNode($self_comp_hit_reg);
					$sc_domain_node->appendChild($struct_seqid_range_node);

					my $ungapped_seqid_range	= ungap_range($self_comp_hit_reg, $$global_opt{gap_tol});
					my $ungapped_seqid_range_aref	= range_expand($ungapped_seqid_range);
					my $ungapped_seqid_range_node	= $domain_xml_doc->createElement('ungapped_seqid_range');
					$ungapped_seqid_range_node->setAttribute('gap_tolerance', $$global_opt{gap_tol});
					$ungapped_seqid_range_node->appendTextNode($ungapped_seqid_range);
					$sc_domain_node->appendChild($ungapped_seqid_range_node);

					my $hit_struct_seqid_range_node = $domain_xml_doc->createElement('hit_struct_seqid_range');
					$hit_struct_seqid_range_node->appendTextNode($self_comp_query_reg);
					$sc_domain_node->appendChild($hit_struct_seqid_range_node);
					
					my $hit_domain_id_node	= $domain_xml_doc->createElement('hit_domain');
					$hit_domain_id_node->setAttribute('ecod_domain_id', $hit_domain_id);
					$hit_domain_id_node->setAttribute('reference', $reference);
					$hit_domain_id_node->setAttribute('uid', $$ref_domain_uid_lookup{$hit_domain_id});
					$sc_domain_node->appendChild($hit_domain_id_node);

					my $pdb_range	= pdb_rangify($sc_hit_reg_aref, $query_pdbnum_aref);
					my $pdb_range_node	= $domain_xml_doc->createElement('pdb_range');
					$pdb_range_node->appendTextNode($pdb_range);
					$sc_domain_node->appendChild($pdb_range_node);

					my $ungapped_pdb_range	= pdb_rangify($ungapped_seqid_range_aref, $query_pdbnum_aref);
					my $ungapped_pdb_range_node	= $domain_xml_doc->createElement('ungapped_pdb_range');
					$ungapped_pdb_range_node->setAttribute('gap_tolerance', $$global_opt{gap_tol});
					$ungapped_pdb_range_node->appendTextNode($ungapped_pdb_range);
					$sc_domain_node->appendChild($ungapped_pdb_range_node);

					if ($self_comp_aligner eq 'dali') { 
						my $dali_score_node	= $domain_xml_doc->createElement('dali_zscore');
						$dali_score_node->appendTextNode($$self_comp_aref[$i]{zscore});
						$sc_domain_node->appendChild($dali_score_node);
					}elsif($self_comp_aligner eq 'hhrepid') { 
						my $hh_prob_node	= $domain_xml_doc->createElement('hhrepid_prob');
						$hh_prob_node->appendTextNode($$self_comp_aref[$i]{prob});
						$sc_domain_node->appendChild($hh_prob_node);
					}else{
						print "WARNING! self comp aligner $self_comp_aligner unknown\n";
					}

					$domain_list_node->appendChild($sc_domain_node);
					range_include($ungapped_seqid_range_aref, $used_seq_aref);
					range_exclude($ungapped_seqid_range_aref, $unused_seq_aref);

					my $post_used = rangify(@$used_seq_aref);
					my $post_unused = rangify(@$unused_seq_aref);
					
					$domain_count++;
					
				}
			}
		}
	}
}

sub find_chainwise_blast_domains  { 
	my $sub = 'find_chainwise_blast_domains';

	my (	$blast_summ_xml_doc,
		$domain_xml_doc,
		$ref_chain_domains,
		$ref_range_cache,
		$ref_domain_uid_lookup,
		$query_struct_seqid_aref,
		$query_seqid_aref,
		$query_pdbnum_aref,
		$used_seq_aref,
		$unused_seq_aref,
		$input_mode,
		$global_opt,
		$reference, 
		$coil_aref
		) = @_;

	my $GAP_TOL = $$global_opt{gap_tol};
	my $domain_count = $domain_xml_doc->exists('//domain') 
		? $domain_xml_doc->findnodes('//domain')->size() + 1
		: 1;
	
	my $domain_list_node = $domain_xml_doc->findnodes('//domain_list')->get_node(1);

	my $query_pdb = $domain_xml_doc->findvalue('//@pdb_id');
	my $query_chain = $domain_xml_doc->findvalue('//@chain_id');

	my $chblastp_run_XPath 	= q{//blast_summ_doc/blast_summ/chain_blast_run[@program='blastp']};
	my $chblastp_run_nodes	= $blast_summ_xml_doc->findnodes($chblastp_run_XPath);
	my $chblastp_run_nodes_size	= $chblastp_run_nodes->size();
	if ($chblastp_run_nodes_size != 1) { 
		print "ERROR! $sub: Odd number of chblastp run nodes found ($chblastp_run_nodes_size), aborting...\n";
	}
	my $chblastp_hit_nodes	= $chblastp_run_nodes->get_node(1)->findnodes('hits/hit');
	foreach my $hit_node ($chblastp_hit_nodes->get_nodelist() ) {

		my $hit_num 	= $hit_node->findvalue('@num');
		my $hit_pdb	= $hit_node->findvalue('@pdb_id');
		my $hit_chain	= $hit_node->findvalue('@chain_id');
		my $hsp_count 	= $hit_node->findvalue('@hsp_count');
		my $evalues 	= $hit_node->findvalue('@evalues');

		if ($hit_chain eq '.') { 
			die "FATAL ERROR! chain blast not correctly configured to handle multi-chain domain hits\n";
		}

		my $evalue = 'Unk';
		if ($hsp_count == 1) { 
			$evalue = $evalues;
		}
		my $query_reg		= $hit_node->findvalue('query_reg');
		#ARG 1 is pos-based, ARG2 is seqid_based, OUTPUT is segid
		my $hit_query_struct_seqid_aref;
		if ($input_mode eq 'seqid') { 
			$hit_query_struct_seqid_aref = range_expand($query_reg);
		}elsif($input_mode eq 'struct_seqid') { 
			$hit_query_struct_seqid_aref = struct_region(range_expand($query_reg), $query_struct_seqid_aref);
		}else{
			die "ERROR! $sub: input_mode $input_mode not recognized\n";
		}
		my $query_struct_seqid_range 	= rangify(@$hit_query_struct_seqid_aref);
		my $hit_reg			= ungap_range($hit_node->findvalue('hit_reg'), $$global_opt{gap_tol});
		my $query_seq		= $hit_node->findvalue('query_seq');
		my $hit_seq		= $hit_node->findvalue('hit_seq');

		if ($DEBUG) { 
			print "DEBUG $sub: chblastp $hit_num $hit_pdb $hit_chain $query_struct_seqid_range\n";
		}
		#Assign if overlaps 
		
		my $unused_coverage 	= region_coverage($query_struct_seqid_aref, $unused_seq_aref);
		if ($unused_coverage < 0.05 || scalar(@$unused_seq_aref) < 10) { 
			last;
		}
		my $query_new_coverage 	= region_coverage($hit_query_struct_seqid_aref, $unused_seq_aref);
		my $query_used_coverage = region_coverage($hit_query_struct_seqid_aref, $used_seq_aref);

		if ($DEBUG) { 
			print "DEBUG $sub: uc $unused_coverage qnc $query_new_coverage qoc $query_used_coverage\n";
		}

		if ($query_new_coverage > $$global_opt{new_coverage_threshold} && $query_used_coverage < $$global_opt{old_coverage_threshold}) { 
			
			print "??$hit_pdb, $hit_chain\n";
			my ($hit_seqid_aref, $hit_struct_seqid_aref, $hit_pdbnum_aref, $hit_asym_id) = pdbml_seq_parse($hit_pdb, $hit_chain);

			if (!$hit_seqid_aref) { next } 
			print "GT: $GAP_TOL gt: $$global_opt{gap_tol}\n";
			#ungap_range_aref(@$hit_struct_seqid_aref, $GAP_TOL); #This wasn't working!?
			my $test = rangify(@$hit_struct_seqid_aref);
			$test = ungap_range($test, $$global_opt{gap_tol});
			$hit_struct_seqid_aref = range_expand($test);


			print "DEFINE chblastp: $sub: $hit_pdb $hit_chain $query_struct_seqid_range $test\n";

			#Calculate hit->query seqid map;
			my @hit_segs = split(",", $hit_reg);
			my @query_segs = split(",", $query_reg);

			my @hit_seqs = split(",", $hit_seq);
			my @query_seqs = split(",", $query_seq);

			my %query_map;
			my %hit_map;
			my %hit_chain_map;
			for (my $i = 0; $i < scalar(@hit_segs); $i++) { 
				printf "%i :  %s\n", scalar(@hit_segs) , $hit_segs[$i];
				my $hit_seg	= $hit_segs[$i];
				#Following line is a bug if library is not structured residues.
				if ($DEBUG) { 
					print "DEBUG $sub: chblastp struct_range_expand 1 $hit_seg\n";
				}
				my $hit_seg_seqid_aref	= struct_region(range_expand($hit_seg), $hit_struct_seqid_aref);
				my @hit_seq	= split("", $hit_seqs[$i]);

				my $query_seg	= $query_segs[$i];
				if ($DEBUG) { 
					printf "DEBUG $sub: chblastp struct_range_expand 2 $query_seg %s\n", rangify(@$hit_query_struct_seqid_aref);
				}
				#Why is this used?
				#my $query_seg_seqid_aref	= struct_region(range_expand($query_seg), $hit_query_struct_seqid_aref);
				my $query_seg_seqid_aref;
				if ($input_mode eq 'seqid') { 
					$query_seg_seqid_aref	= isect(range_expand($query_seg), $query_struct_seqid_aref);
				}elsif($input_mode eq 'struct_seqid') { 
					$query_seg_seqid_aref	= struct_region(range_expand($query_seg), $query_struct_seqid_aref);
				}
				if ($DEBUG) { 
					printf "DEBUG $sub: chblastp query_seg_seqid %s\n", rangify(@$query_seg_seqid_aref);
				}
				my @query_seq	= split("",$query_seqs[$i]);

				my $hit_start = $$hit_seg_seqid_aref[0];
				my $query_start	= $$query_seg_seqid_aref[0];
				
				my $hit_pos 	= 0;
				my $query_pos 	= 0;

				for (my $j = 0; $j < scalar(@query_seq); $j++) { 
					if ($query_seq[$j] =~ /[A-Z]/ && $hit_seq[$j] eq '-') { 
						$query_pos++;
					}elsif($query_seq[$j] eq '-' && $hit_seq[$j] =~ /[A-Z]/) { 
						$hit_pos++;
					}elsif($query_seq[$j] =~ /[A-Z]/ && $hit_seq[$j] =~ /[A-Z]/) { 
						if ($$query_seg_seqid_aref[$query_pos] && $$hit_seg_seqid_aref[$hit_pos]) { 
							$query_map{$$query_seg_seqid_aref[$query_pos] } = $$hit_seg_seqid_aref[$hit_pos];
							$hit_map{$$hit_seg_seqid_aref[$hit_pos]} = $$query_seg_seqid_aref[$query_pos];
						}
						$query_pos++;
						$hit_pos++;
					}else{
						#Foul
					}
				}
			}
				
			HIT_DOMAIN:
			foreach my $hit_domain_uid (@{$$ref_chain_domains{$hit_pdb}{$hit_chain}}) { 
				
				my $ref_domain_id = $$ref_range_cache{$hit_domain_uid}{ecod_domain_id};
				my $hit_domain_id_node	= $domain_xml_doc->createElement('hit_domain');
				$hit_domain_id_node->setAttribute('ecod_domain_id', $ref_domain_id);
				$hit_domain_id_node->setAttribute('uid', $hit_domain_uid);
				$hit_domain_id_node->setAttribute('reference', $reference);

				my $ref_seqid_range = multi_chain_ungap_range($$ref_range_cache{$hit_domain_uid}{seqid_range}, $$global_opt{gap_tol});

				#print "$hit_domain_uid $ref_domain_id $ref_seqid_range\n";

				if ($DEBUG) { 
					print "DEBUG: $sub: chblastp_domain $ref_domain_id $ref_seqid_range\n";
				}
				my $domain_node		= $domain_xml_doc->createElement('domain');
				$domain_node->appendChild($hit_domain_id_node);
				my $domain_id	= "e". lc( $query_pdb) . "$query_chain$domain_count";
				$domain_node->setAttribute('domain_id', $domain_id);

#				my $method_node	= $domain_xml_doc->createElement('method');
#				my $method	= 'chblastp';
#				$method_node->appendTextNode($method);
#				$domain_node->appendChild($method_node);
#
				$domain_node->appendTextChild('method', 'chblastp');

				#REPLACE WITH HIT DOMAIN DERIVED RANGE
				my @ref_segs = split(/,/, $ref_seqid_range);
				my @query_struct_seqid_range;
				my @query_struct_chain_range;
				my @hit_struct_seqid_range;
				my $warning = 0;
				foreach my $ref_seg (@ref_segs) { 
					#print "REFSEG $ref_seg\n";
					if ($ref_seg =~ /(.+):(\-?\d+\-\-?\d+)/) { 
						my $ref_seg_chain = $1;
						if ($ref_seg_chain ne $hit_chain) { 
							print "WARNING! ref seg chain $ref_seg_chain is not the same as expected hit chain $hit_chain, multi chain problems...\n";
						}
						my $ref_seg_range = $2;
						my $ref_seg_range_aref = range_expand($ref_seg_range);
						foreach my $ref_seqid (@$ref_seg_range_aref) { 
							if ($hit_map{$ref_seqid}) { 
								#print "MESS $ref_seqid $hit_map{$ref_seqid}\n";
								push (@query_struct_seqid_range, $hit_map{$ref_seqid});
								push (@query_struct_chain_range, $hit_chain_map{$ref_seqid});
							}else{
								print "WARNING! No map for $ref_seqid\n";
								$warning++;
							}
							push (@hit_struct_seqid_range, $ref_seqid);
						}
					}else{
						print "WARNING! $sub: Ref seg looks weird $ref_seg\n";
					}
				}
				if ($warning) { 
					print "WARNING! $warning maps failed\n";
				}
				#Prevent very short chblastp domains (does not solve fragment problem, just too short domain problem)
				if (scalar @query_struct_seqid_range < $$global_opt{gap_tol}) { 
					print "WARNING! Short chblastp domain, skipping...\n";
					next HIT_DOMAIN;
				}
				print "DEBUG: scalar qssr: " . scalar(@query_struct_seqid_range) . "\n";

				my $struct_seqid_range_node	= $domain_xml_doc->createElement('struct_seqid_range');
				#$struct_seqid_range_node->appendTextNode($query_struct_seqid_range);	
				if (scalar(@query_struct_seqid_range) == 0) { 
					print "WARNING! $sub: No struct seqid for $query_pdb $query_chain chblastp, serious problem...\n";
					next;
				}
				my $query_struct_seqid_range =rangify(@query_struct_seqid_range);
				$struct_seqid_range_node->appendTextNode($query_struct_seqid_range);	
				$domain_node->appendChild($struct_seqid_range_node);

				my $ungapped_seqid_range = ungap_range(rangify(@query_struct_seqid_range), $$global_opt{gap_tol});
				my $ungapped_seqid_range_aref = range_expand($ungapped_seqid_range);
				my $ungapped_seqid_range_node	= $domain_xml_doc->createElement('ungapped_seqid_range');
				$ungapped_seqid_range_node->setAttribute('gap_tolerance', $$global_opt{gap_tol});
				$ungapped_seqid_range_node->appendTextNode($ungapped_seqid_range);
				$domain_node->appendChild($ungapped_seqid_range_node);

				#REPLACE WITH HIT DOMAIN RANGE (THAT OVERLAPS ALLGNED REGION)
				my $hit_struct_seqid_range_node	 = $domain_xml_doc->createElement('hit_struct_seqid_range');
				my $hit_struct_seqid_range = scopify_range(rangify(@hit_struct_seqid_range), $hit_chain);
				$hit_struct_seqid_range_node->appendTextNode($hit_struct_seqid_range);
				$domain_node->appendChild($hit_struct_seqid_range_node);

				my $pdb_range = pdb_rangify(\@query_struct_seqid_range, $query_pdbnum_aref);
				if ($pdb_range eq 0) {  next HIT_DOMAIN } 
				my $pdb_range_node	= $domain_xml_doc->createElement('pdb_range');
				$pdb_range_node->appendTextNode($pdb_range);
				$domain_node->appendChild($pdb_range_node);

				my $ungapped_pdb_range = pdb_rangify($ungapped_seqid_range_aref, $query_pdbnum_aref);
				my $ungapped_pdb_range_node	= $domain_xml_doc->createElement('ungapped_pdb_range');
				$ungapped_pdb_range_node->setAttribute('gap_tolerance', $$global_opt{gap_tol});
				$ungapped_pdb_range_node->appendTextNode($ungapped_pdb_range);
				$domain_node->appendChild($ungapped_pdb_range_node);

				my $blast_score_node	= $domain_xml_doc->createElement('blast_eval');
				$blast_score_node->appendTextNode($evalue);
				$domain_node->appendChild($blast_score_node);

				$domain_list_node->appendChild($domain_node);
				
				#REPLACE WITH DOMAIN BASED RANGE
				range_exclude(range_expand($query_struct_seqid_range), $unused_seq_aref);
				range_include(range_expand($query_struct_seqid_range), $used_seq_aref);
				
				$domain_count++;
			}
		}
	}

}

sub generate_dali_fillin { 
	my ( 
		$blast_summ_xml_doc,		#1
		$domain_xml_doc, 		#2
		$ref_chain_domains,		#3
		$ref_domain_uid_lookup,		#4
		$ref_range_cache,		#5
		$query_struct_seqid_aref,	#6	
		$query_seqid_aref,		#7
		$query_pdbnum_aref,		#8
		$query_asym_id,
		$used_seq_aref,			#9
		$unused_seq_aref,		#10	
		$input_mode, 			#11
		$global_opt,			#12
		$reference,			#13
		$pdb_chain_dump_dir)		#14
				= @_;

	#Do domains exist to  hit chains with multiple domains (included some unused).

	my (%used_hit_uids, %possible_hit_uids);
	my (%hit_stats, %hit_chains);

	my $query_pdb		= $domain_xml_doc->findvalue('//@pdb_id');
	my $query_chain		= $domain_xml_doc->findvalue('//@chain_id');
	my $query_pdb_chain = $query_pdb . "_" . $query_chain;	
	#Generate dali jobs for query chain against possible hit uids not used 
	my $query_pdb_chain_fn = "$pdb_chain_dump_dir/$query_pdb_chain.pdb"; 
	if (!-f $query_pdb_chain_fn || $FORCE_OVERWRITE) { #generate query pdb 
		print "gen: $query_pdb_chain_fn\n";
		generate_query_pdb($query_pdb_chain_fn, $query_pdb, $query_chain);
	}else{
		if ($DEBUG) { 
			print "DEBUG: $query_pdb_chain_fn exists!\n";
		}
	}

	foreach my $domain (find_domain_nodes($domain_xml_doc)) { 
		my $hit_domain_id 	= $domain->findvalue('hit_domain/@ecod_domain_id');
		my $domain_id 		= $domain->findvalue('@domain_id');

		my $hit_uid = $$ref_domain_uid_lookup{$hit_domain_id};

		$used_hit_uids{$hit_uid}++;

		my $hit_pdb 	= $$ref_range_cache{$hit_uid}{pdb};
		my $hit_chain 	= $$ref_range_cache{$hit_uid}{chain};

		print "D: $domain_id $query_pdb $query_chain => $hit_domain_id/$hit_uid $hit_pdb $hit_chain\n";

		if (scalar (@{$$ref_chain_domains{$hit_pdb}{$hit_chain}}) > 1) { 
			foreach my $hit_uid (@{$$ref_chain_domains{$hit_pdb}{$hit_chain}}) { 
				if (!$used_hit_uids{$hit_uid}) { 
					$possible_hit_uids{$hit_uid}++;
				}
			}
		}else{
			#print "No fillin for $domain_id\n";
		}
	}
	my @files;
	foreach my $hit_uid (keys %possible_hit_uids) { 
		my $hit_pdb 	= $$ref_range_cache{$hit_uid}{pdb};
		my $hit_chain 	= $$ref_range_cache{$hit_uid}{chain};
		$hit_stats{$hit_uid}{pdb} = $hit_pdb;
		$hit_stats{$hit_uid}{chain} = $hit_chain;
		if (!$hit_chains{$hit_pdb}{$hit_chain}) { 
			($hit_chains{$hit_pdb}{$hit_chain}{seqid_aref},
			$hit_chains{$hit_pdb}{$hit_chain}{struct_seqid_aref},
			$hit_chains{$hit_pdb}{$hit_chain}{pdbnum_aref},
			$hit_chains{$hit_pdb}{$hit_chain}{asym_id}) 	= pdbml_seq_parse($hit_pdb, $hit_chain);
			$hit_chains{$hit_pdb}{$hit_chain}{seq} 		= pdbml_seq_fetch($hit_pdb, $hit_chains{$hit_pdb}{$hit_chain}{asym_id}, $hit_chain, $hit_chains{$hit_pdb}{$hit_chain}{seqid_aref});
		}
		$hit_stats{$hit_uid}{seqid_aref} 	= $hit_chains{$hit_pdb}{$hit_chain}{seqid_aref};
		$hit_stats{$hit_uid}{struct_seqid_aref} = $hit_chains{$hit_pdb}{$hit_chain}{struct_seqid_aref};
		$hit_stats{$hit_uid}{pdbnum_aref}	= $hit_chains{$hit_pdb}{$hit_chain}{pdbnum_aref};
		$hit_stats{$hit_uid}{asym_id}		= $hit_chains{$hit_pdb}{$hit_chain}{asym_id};
		$hit_stats{$hit_uid}{seq}		= $hit_chains{$hit_pdb}{$hit_chain}{seq};
		$hit_stats{$hit_uid}{range}		= $$ref_range_cache{$hit_uid}{seqid_range};
		$hit_stats{$hit_uid}{ecod_domain_id}	= $$ref_range_cache{$hit_uid}{ecod_domain_id};

		my $DaliLite_exe = $DALI_EXE;
		if (!-f $DaliLite_exe) { 
			die "ERROR! Could not find DaliLite executable $DaliLite_exe.\n";
		}

		my $dali_tmp_dir = "$pdb_chain_dump_dir/dali_tmp";
		if (! -d $dali_tmp_dir) { 
			mkdir($dali_tmp_dir);
		}
		my @job_ids;

		my $short_uid = substr($hit_uid, 2, 5);
		#my $hit_domain_id = $reverse{$uid};
		my $dali_output_fn = "q.$query_pdb_chain.$hit_uid.$reference.dali";
		push (@files, "$dali_tmp_dir/$dali_output_fn");
		if (-f "$dali_tmp_dir/$dali_output_fn") { 
			next;
		}
		#Can't do that you ninny
		print "immediately!\n";
		immediate_dali($dali_tmp_dir, "../../$query_pdb_chain.pdb", "$DOMAIN_DATA_DIR/$short_uid/$hit_uid/$hit_uid.seqres.pdb", $dali_output_fn);
		#my $job_fn = jobify_dali($dali_tmp_dir, "../../$query_pdb_chain.pdb", "$DOMAIN_DATA_DIR/$short_uid/$hit_uid/$hit_uid.seqres.pdb", $dali_output_fn);
		#my $job_id = qsub($job_fn);
		#push (@job_ids, $job_id);

		#while (qstat_wait_list(\@job_ids)) { 
		#	print "SLEEPING...\n";
		#	sleep(5);
		#}
	}

	my $dali_fillin_summ_fn = "$pdb_chain_dump_dir/$query_pdb_chain.$reference.fillin.dali_summ.xml";
	if (0 && -f $dali_fillin_summ_fn) { 
		print "WARNING! Dali fillin summary file found ($dali_fillin_summ_fn), skipping...\n";
		return;
	}
	#dali_summ

	my $query_seq_href = pdbml_seq_fetch($query_pdb, $query_asym_id, $query_chain, $query_seqid_aref);

	my %query;
	$query{pdb}	= $query_pdb;
	$query{chain}	= $query_chain;
	$query{asym_id} = $query_asym_id;
	$query{struct_seqid_aref} = $query_struct_seqid_aref;
	$query{pdbnum_aref} = $query_pdbnum_aref;
	$query{seq}	= $query_seq_href;

	dali_summ(\@files, $dali_fillin_summ_fn, \%query, \%hit_stats);
	
	return $dali_fillin_summ_fn;
}

sub find_dali_fillin_domains { 
	my $sub = 'find_dali_fillin_domains';
	my ( 
		$domain_xml_doc,
		$dali_summ_xml_doc,
		$ref_chain_domains,
		$ref_domain_uid_lookup,
		$ref_range_cache,
		$query_struct_seqid_aref,
		$query_pdbnum_aref,
		$query_asym_id,
		$used_seq_aref,
		$unused_seq_aref,
		$input_mode, 
		$global_opt,
		$reference,
		$dir,
		) = @_;

	if (ref $query_pdbnum_aref) { 
		print "pdbnum $query_pdbnum_aref\n";
	}else{
		die;
	}

	my $domain_count = $domain_xml_doc->exists('//domain') 
		? $domain_xml_doc->findnodes('//domain')->size() + 1
		: 1;
	
	my $domain_list_node = $domain_xml_doc->findnodes('//domain_list')->get_node(1);

	my $query_pdb = $domain_xml_doc->findvalue('//@pdb_id');
	my $query_chain = $domain_xml_doc->findvalue('//@chain_id');


	my @hits;
	my $i = 0;
	foreach my $hit_node (find_dali_hit_nodes($dali_summ_xml_doc)) { 

		my ($z, $rmsd, $id) 		= get_dali_scores($hit_node);
		my $dali_fn			= get_dali_fn($hit_node);
		my ($hit_reg, $query_reg, $coverage) = get_dali_regions($hit_node);
		my ($uid, $ecod_domain_id) 	= get_ids($hit_node);

		$hits[$i]{z}	= $z;
		$hits[$i]{rmsd} = $rmsd;
		$hits[$i]{id}	= $id;

		$hits[$i]{dali_fn} = $dali_fn;

		$hits[$i]{hit_reg} 	= $hit_reg;
		$hits[$i]{query_reg} 	= $query_reg;
		$hits[$i]{coverage} 	= $coverage;

		$hits[$i]{uid}		= $uid;
		$hits[$i]{ecod_domain_id} = $ecod_domain_id;
		$i++;

	}
	printf "$sub: hit_nodes %i\n", scalar(@hits);

	@hits = sort { $b->{z} <=> $a->{z} } @hits;

	for (my $i = 0; $i < scalar(@hits); $i++ ) {
		my $dali_fn 	= $hits[$i]{dali_fn};

		my $hit_ecod_domain_id = $hits[$i]{ecod_domain_id};
		my $hit_uid		= $hits[$i]{uid};
		my $hit_pdb		= $$ref_range_cache{$hit_uid}{pdb};
		my $hit_chain	= $$ref_range_cache{$hit_uid}{chain};

		my $hit_range	= $$ref_range_cache{$hit_uid}{seqid_range};
		
		my $mc_hit = 0;
		my ($hit_range_str,$hit_chain_str) = scop_range_split($hit_range);
		if ($hit_chain_str =~ /,/) { 
			$hit_chain = $hit_chain_str;
			$mc_hit = 1;
		}

		my $hit_reg	= $hits[$i]{hit_reg};
		my $query_reg 	= $hits[$i]{query_reg};

		my $dali_z	= $hits[$i]{z};
		my $dali_rmsd	= $hits[$i]{rmsd};
		my $dali_id	= $hits[$i]{id};

		my $hit_coverage = $hits[$i]{coverage};

		my $query_domain_seqid_aref 	= struct_region(range_expand($query_reg), $query_struct_seqid_aref);
		my $query_domain_seqid_range 	= ungap_range(rangify(@$query_domain_seqid_aref), $$global_opt{gap_tol});
		my $query_domain_pdb_range 		= ungap_range(pdb_rangify($query_domain_seqid_aref, $query_pdbnum_aref), $$global_opt{gap_tol});

		my $query_new_coverage 	= region_coverage($query_domain_seqid_aref, $unused_seq_aref);
		my $query_used_coverage	= region_coverage($query_domain_seqid_aref, $used_seq_aref);

		printf "$sub: %.2f %.2f %.2f\n", $query_new_coverage, $query_used_coverage, $hit_coverage;

		if (($query_new_coverage > $$global_opt{new_coverage_threshold} &&
			$query_used_coverage < $$global_opt{old_coverage_threshold} && 
			$hit_coverage > $$global_opt{hit_coverage_threshold}) 
			||
			($query_new_coverage > $$global_opt{new_coverage_threshold} && 
			$query_used_coverage < $$global_opt{old_coverage_threshold} && 
			$dali_z > $$global_opt{frag_z_cutoff}) )  { 

			print "define! $hit_reg $query_reg\n";

			my $fragment_defined = $hit_coverage <= $$global_opt{hit_coverage_threshold} ? 1 : 0 ;
	
			my $hit_struct_seqid_range;
			if ($mc_hit) { 
				my ($hit_seqid_aref, $hit_struct_seqid_aref, $hit_pdbnum_href, $hit_asym_id, $hit_chain_aref) = pdbml_mc_seq_parse($hit_pdb, $hit_chain);
				if (!$hit_seqid_aref) { 
					print "WARNING! mc hit parse failure, $hit_pdb, $hit_chain\n";
					next;
				}
				my ($hit_reg_struct_seqid_aref, $hit_reg_struct_chain_aref) = multi_chain_struct_region(range_expand($hit_reg), $hit_struct_seqid_aref, $hit_chain_aref);
				$hit_struct_seqid_range = multi_chain_rangify($hit_reg_struct_seqid_aref, $hit_reg_struct_chain_aref);
			}else{
				my ($hit_seqid_aref, $hit_struct_seqid_aref, $hit_pdbnum_aref, $hit_asym_id) = pdbml_seq_parse($hit_pdb, $hit_chain);
				if (!$hit_seqid_aref) { 
					print "WARNING! seq parse fail, $hit_pdb, $hit_chain\n";
					next;
				}
				my $hit_reg_struct_seqid_aref = struct_region(range_expand($hit_reg), $hit_struct_seqid_aref);
				my $hit_reg_struct_seqid_range = rangify(@$hit_reg_struct_seqid_aref);
				$hit_struct_seqid_range = scopify_range($hit_chain, $hit_reg_struct_seqid_range);
			}

			my $domain_node	= $domain_xml_doc->createElement('domain');
			my $domain_id	= "e" . lc($query_pdb) . "$query_chain$domain_count";
			$domain_node->setAttribute('domain_id', $domain_id);

			my $method_node	= $domain_xml_doc->createElement('method');
			my $method 	= 'dali_fillin';
			$method_node->appendTextNode($method);
			if ($fragment_defined) { 
				$method_node->setAttribute('fragment', 'true');
			}
			$domain_node->appendChild($method_node);

			my $struct_seqid_range_node	= $domain_xml_doc->createElement('struct_seqid_range');
			$struct_seqid_range_node->appendTextNode(rangify(@$query_domain_seqid_aref));
			$domain_node->appendChild($struct_seqid_range_node);

			my $ungapped_seqid_range = ungap_range(rangify(@$query_domain_seqid_aref), $$global_opt{gap_tol});
			my $ungapped_seqid_range_aref = range_expand($ungapped_seqid_range);
			my $ungapped_seqid_range_node	= $domain_xml_doc->createElement('ungapped_seqid_range');
			$ungapped_seqid_range_node->setAttribute('gap_tolerance', $$global_opt{gap_tol});
			$ungapped_seqid_range_node->appendTextNode($ungapped_seqid_range);
			$domain_node->appendChild($ungapped_seqid_range_node);

			my $hit_struct_seqid_range_node	= $domain_xml_doc->createElement('hit_struct_seqid_range');
			$hit_struct_seqid_range_node->appendTextNode($hit_struct_seqid_range);
			$domain_node->appendChild($hit_struct_seqid_range_node);

			my $pdb_range	= pdb_rangify($query_domain_seqid_aref, $query_pdbnum_aref);
			my $pdb_range_node	= $domain_xml_doc->createElement('pdb_range');
			$pdb_range_node->appendTextNode($pdb_range);
			$domain_node->appendChild($pdb_range_node);

			my $ungapped_pdb_range = pdb_rangify($ungapped_seqid_range_aref, $query_pdbnum_aref);
			my $ungapped_pdb_range_node	= $domain_xml_doc->createElement('ungapped_pdb_range');
			$ungapped_pdb_range_node->setAttribute('gap_toleraance', $GAP_TOL);
			$ungapped_pdb_range_node->appendTextNode($ungapped_pdb_range);
			$domain_node->appendChild($ungapped_pdb_range_node);

			my $hit_domain_node	= $domain_xml_doc->createElement('hit_domain');
			$hit_domain_node->setAttribute('ecod_domain_id', $hit_ecod_domain_id);
			$hit_domain_node->setAttribute('uid', $hit_uid);
			$hit_domain_node->setAttribute('reference', $reference);
			$domain_node->appendChild($hit_domain_node);

			my $dali_score_node	= $domain_xml_doc->createElement('dali_score');
			$dali_score_node->setAttribute('z_score', $dali_z);
			$dali_score_node->setAttribute('rmsd', $dali_rmsd);
			$dali_score_node->setAttribute('identity', $dali_id);
			$dali_score_node->setAttribute('coverage', $hit_coverage);
			$domain_node->appendChild($dali_score_node);

			$domain_list_node->appendChild($domain_node);

			range_exclude($query_domain_seqid_aref, $unused_seq_aref);
			range_include($query_domain_seqid_aref, $used_seq_aref);

			$domain_count++;
	
		}
	}
}

sub range_decompose { 
	#hqr => hit_query_reg
	#qssa => query_struct_seqid_aref
	#qsa => query_seqid_aref
	#qsca => query_struct_chain_aref
	#qca => query_chain_aref;
	my ($hqr, $qssa, $qsa, $qsca, $qca, $input_mode) = @_;

	my ($hqssa, $hqsca) = $input_mode eq 'struct_seqid' ?  
				multi_chain_struct_region(range_expand($hqr), $qssa, $qsca) :
				multi_chain_struct_region(range_expand($hqr), $qsa, $qca);

	return ($hqssa, $hqsca);
}

sub bidirectional_coverage { 
	my ($range_aref1, $range_aref2, $c) = @_;
	if ($range_aref1 eq 'NA' or $range_aref2 eq 'NA') { return 0 } 
	my $c1 = region_coverage($range_aref1, $range_aref2);
	my $c2 = region_coverage($range_aref2, $range_aref1);
	if ($c1 > $c && $c2 > $c) { 
		return 1;
	}
	return 0;
}
sub isect { 
	my $sub = 'isect';
	my ($aref1, $aref2) = @_;

	my %seen;
	my @i;
	foreach my $a1 (@$aref1) { 
		$seen{$a1}++;
	}
	foreach my $a2 (@$aref2) { 
		if ($seen{$a2}) { 
			push (@i, $a2);
		}
	}
	@i = sort { $a <=> $b } @i;
	return \@i;
}

sub union { 
	my $sub = 'unioin';
	my ($aref1, $aref2) = @_;

	my %s;
	foreach my $a1 (@$aref1) { 
		$s{$a1}++;
	}
	foreach my $a2 (@$aref2) { 
		$s{$a2}++;
	}
	my @u = sort {$a <=> $b} keys %s;
	return \@u;
}
1;

