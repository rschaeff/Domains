package Domains::PDBml;
require Exporter;

use warnings;
use strict;

use XML::LibXML;
use XML::LibXML::Reader;

use Domains::Range;


our @ISA = qw(Exporter);
our @EXPORT = 
	(
	"&get_ppss_seq_ids",
	"&pdbml_annot_parse", 
	"&pdbml_date_method",
	 "&pdbml_entity_poly",
	 "&pdbml_doc_seq_parse",
	 "&pdbml_fetch_chains",
	 "&pdbml_quick_fetch_chains",
	"&pdbml_seq_parse",
	"&pdbml_obs_seq_parse",
	"&pdbml_mc_seq_parse",
	"&pdbml_coord_parse",
	"&pdbml_coord_ligand_parse_pull",
	"&pdbml_coord_parse_pull",
	"&pdbml_mc_coord_parse_pull",
	"&pdbml_load",
	"&pdbml_asym_annotate_parse",
	"&pdbml_fasta_fetch",
	"&pdbml_obs_fasta_fetch",
	"&pdbml_seq_fetch",
	"&pdbml2pdb_cartn" );

my $DEBUG = 0;

my $MAX_FILE_ATTEMPTS = 5;

my $pdb_noatom_top_dir = '/usr2/pdb/data/structures/divided/XML-noatom';
(-d $pdb_noatom_top_dir) or die "Could not find pdb top dir: $pdb_noatom_top_dir\n";
my $pdb_obs_top_dir = '/usr2/pdb/data/structures/obsolete/XML/';
my $pdb_atom_top_dir = '/usr2/pdb/data/structures/divided/XML';
(-d $pdb_atom_top_dir) or die "Could not find pdb top dir: $pdb_atom_top_dir\n";
#my $pdb_cache_top_dir = '/home/rschaeff_1/pdb_junk_bin';
#(-d $pdb_cache_top_dir) or die "Coudl not find pdb cache dir: $pdb_cache_top_dir\n";
sub get_ppss_seq_ids { 
	($_[0]->findvalue('@seq_id'),
		$_[0]->findvalue('PDBx:auth_seq_num'),
		$_[0]->findvalue('PDBx:pdb_seq_num'));
}
sub pdbml_fetch_chains { 
        my $sub = 'pdbml_fetch_chains';
        my $pdb = shift (@_);

        if ($DEBUG) { 
                print "DEBUG $sub: $pdb\n";
        }
        chomp $pdb;

        my $pdbml = pdbml_load($pdb);

        my $entity_poly_XPath = '//PDBx:entity_polyCategory/PDBx:entity_poly';

        my @chains;
		my %types;
        foreach my $ep_node ($pdbml->findnodes($entity_poly_XPath)->get_nodelist()) { 

			my $type = $ep_node->findvalue('PDBx:type');
			my $entity_id = $ep_node->findvalue('@entity_id');

			if ($type ne 'polypeptide(L)') { next } 

			my $strand_string = $ep_node->findvalue('PDBx:pdbx_strand_id');
			#print "s:$strand_string\n";
			my @strands = split(/,\s?/, $strand_string);
			push (@chains, @strands);
			foreach my $strand (@strands) { 
				$types{$strand} = $type;
			}
        }

	my $pdbx_poly_seq_schemeCategoryXPath = "//PDBx:pdbx_poly_seq_schemeCategory/PDBx:pdbx_poly_seq_scheme";
	my $pdbx_poly_seq_scheme_nodes = $pdbml->findnodes($pdbx_poly_seq_schemeCategoryXPath);

	my @seq_ids;
	my %seq_ids;
	my %seen_seq_ids;
	my %struct_seq_ids;
	my %pdb_seq_nums;

	foreach my $node ($pdbx_poly_seq_scheme_nodes->get_nodelist() ) { 

		my $seq_id	= $node->findvalue('@seq_id');
		my $asym_id	= $node->findvalue('@asym_id');

		my $auth_seq_num	= $node->findvalue('PDBx:auth_seq_num');
		my $pdb_seq_num		= $node->findvalue('PDBx:pdb_seq_num');
		
		my $hetero		= $node->findvalue('PDBx:hetero');

		my $strand_id 	= $node->findvalue('PDBx:pdb_strand_id');

		if ($seen_seq_ids{$strand_id}{$seq_id}) { 
			print "WARNING! $sub: duplicate seq id $seq_id in asym_id $asym_id\n";
			next;
		}else{
			$seen_seq_ids{$strand_id}{$seq_id}++;
		}
		my $pdb_ins_code;
		if ($node->findvalue('PDBx:pdb_ins_code/@xsi:nil') ne 'true') { 
			$pdb_ins_code = $node->findvalue('PDBx:pdb_ins_code');
		}
		if (defined $pdb_seq_num) { 
			if ($pdb_ins_code) { 
				$pdb_seq_nums{$strand_id}[$seq_id] = $pdb_seq_num . $pdb_ins_code;
			}else{
				$pdb_seq_nums{$strand_id}[$seq_id] = $pdb_seq_num;
			}
		}
		if ($auth_seq_num =~ /\d+/) { 
			push (@{$struct_seq_ids{$strand_id}}, $seq_id);
		}
		push (@{$seq_ids{$strand_id}}, $seq_id);

	}

        return (\@chains, \%seq_ids, \%struct_seq_ids, \%pdb_seq_nums, \%types);
}
sub pdbml_quick_fetch_chains { 
	my $sub = 'fetch_chains';
	my $pdb = shift (@_);

	if ($DEBUG) { 
		print "DEBUG $sub: $pdb\n";
	}
	chomp $pdb;

	my $pdbml = pdbml_load($pdb);

	my $entity_poly_XPath = '//PDBx:entity_polyCategory/PDBx:entity_poly';

	my @chains;
	foreach my $ep_node ($pdbml->findnodes($entity_poly_XPath)->get_nodelist()) { 

		my $type = $ep_node->findvalue('PDBx:type');
		if ($type ne 'polypeptide(L)') { next } 

		my $strand_string = $ep_node->findvalue('PDBx:pdbx_strand_id');
		#print "s:$strand_string\n";
		my @strands = split(/,\s?/, $strand_string);
		push (@chains, @strands);
	}
	return \@chains;
}


sub pdbml_load { 
	my $sub = 'pdbml_load';
	my ($pdb4) = @_;
	my $xml_fh;
	$pdb4 =~ /\w(\w{2})\w/;

	my $two = $1;
	if (-f "$pdb_noatom_top_dir/$two/$pdb4-noatom.xml") { 
		if(!open ($xml_fh, "$pdb_noatom_top_dir/$two/$pdb4-noatom.xml")) { 
			my $i = 0;
			while (!open($xml_fh, "$pdb_noatom_top_dir/$two/$pdb4-noatom.xml") && $i < $MAX_FILE_ATTEMPTS) { 
				sleep(1);
				$i++
			}
			if ($i == $MAX_FILE_ATTEMPTS) { 
				die "ERROR $sub: Could not open $pdb_noatom_top_dir/$two/$pdb4-noatom.xml for reading:$!\n";
			}
		}
	}elsif(-f  "$pdb_noatom_top_dir/$two/$pdb4-noatom.xml.gz") { 
		if (!open ($xml_fh, "gzip -dc $pdb_noatom_top_dir/$two/$pdb4-noatom.xml.gz |") ) { 
			my $i = 0;
			while (!open ($xml_fh, "gzip -dc $pdb_noatom_top_dir/$two/$pdb4-noatom.xml.gz |") && $i < $MAX_FILE_ATTEMPTS ) { 
				sleep(1);
				$i++;
			}
			if ($i == $MAX_FILE_ATTEMPTS) {
				die "ERROR $sub: Could not open $pdb_noatom_top_dir/$two/$pdb4-noatom.xml.gz for reading:$!\n";
			}
		}
	}elsif(-f  "$pdb_obs_top_dir/$two/$pdb4.xml.gz") { 

			if (!open ($xml_fh, "gzip -dc $pdb_obs_top_dir/$two/$pdb4.xml.gz |") ) { 
				my $i = 0;
				while (!open ($xml_fh, "gzip -dc $pdb_obs_top_dir/$two/$pdb4.xml.gz |") && $i < $MAX_FILE_ATTEMPTS ) { 
					sleep(1);
					$i++;
				}
				if ($i == $MAX_FILE_ATTEMPTS) {
					die "ERROR $sub: Could not open $pdb_obs_top_dir/$two/$pdb4.xml.gz for reading:$!\n";
				}
			}
	}else{
		print "WARNING! $sub: PDBml file not found for $pdb4 in $pdb_noatom_top_dir\n";
		return 0;
	}

	my $pdbml = XML::LibXML->load_xml (
		IO => $xml_fh);
	if ($DEBUG) { 
		printf "DEBUG $sub: %s\n", $pdbml->findvalue('//PDBx:entryCategory/PDBx:entry/@id')
	}
	close $xml_fh;
	return $pdbml;
}


sub pdbml_entity_poly { 
	my $sub = 'pdbml_entity_poly';
	my ($pdb4) = @_;
	
	$pdb4 = lc($pdb4);
	my $two;
	if ($pdb4 =~ /\w(\w{2})\w/) { 
		$two = $1;
	}else{
		print "ERROR $sub: $pdb4 does not regexp\n";
		return 0;
	}

	if ($DEBUG) { 
		print "DEBUG $sub: $pdb4\n";
	}
	
	my $pdbml = pdbml_load($pdb4);

	if (!$pdbml) { return 0;  }

	my $entity_poly_XPath = '//PDBx:entity_polyCategory/PDBx:entity_poly';

	my $entity_poly_nodes	= $pdbml->findnodes($entity_poly_XPath);
	
	if ($DEBUG) { 
		printf "DEBUG $sub: %i\n", $entity_poly_nodes->size();
	}

	my @entity_poly;
	foreach my $ep_node ($entity_poly_nodes->get_nodelist() ) { 

		my $entity_id		= $ep_node->findvalue('@entity_id');
		my $nstd_linkage	= $ep_node->findvalue('PDBx:nstd_linkage');
		my $nstd_monomer	= $ep_node->findvalue('PDBx:nstd_monomer');
		my $pdbx_strand_id	= $ep_node->findvalue('PDBx:pdbx_strand_id');
		my $type		= $ep_node->findvalue('PDBx:type');

		$entity_poly[$entity_id]{nstd_linkage}	= $nstd_linkage;
		$entity_poly[$entity_id]{nstd_monomer}	= $nstd_monomer;
		$entity_poly[$entity_id]{pdbx_strand_id}= $pdbx_strand_id;	
		$entity_poly[$entity_id]{type}		= $type;

	}
	if ($DEBUG) { 
		printf "DEBUG $sub: scalar %i\n", scalar(@entity_poly);
	}

	return \@entity_poly;
}

sub pdbml_doc_seq_parse { 
	my $sub = 'pdbml_doc_seq_parse';

	my ($pdbml_xml_doc, $pdb4, $chain) = @_;
	my $pdbx_poly_seq_schemeCategoryXPath = "//PDBx:pdbx_poly_seq_schemeCategory/PDBx:pdbx_poly_seq_scheme/PDBx:pdb_strand_id[text()=\'$chain\']";
	my $pdbx_poly_seq_scheme_nodes = $pdbml_xml_doc->findnodes($pdbx_poly_seq_schemeCategoryXPath);

	my @seq_ids;
	my %seq_ids;
	my @struct_seq_ids;
	my @pdb_seq_nums;

	if ($DEBUG) { 
		printf "DEBUG: $pdb4 $chain seqid nodes %i\n", scalar($pdbx_poly_seq_scheme_nodes->get_nodelist());
	}

	my %target_asym_id;
	my $target_asym_id = 'NA';
	foreach my $node ($pdbx_poly_seq_scheme_nodes->get_nodelist()) { 
		my $parent_node	= $node->parentNode;
		my $auth_seq_num	= $parent_node->findvalue("PDBx:auth_seq_num");
		my $seq_id		= $parent_node->findvalue('@seq_id');
		if ($seq_ids{$seq_id}) { 
			print "WARNING! $sub: duplicate seq id $seq_id, heterogeneity? skipping gt first occurence\n";
			next;
		}else{
			$seq_ids{$seq_id}++;
		}

		my $pdb_seq_num 	= $parent_node->findvalue("PDBx:pdb_seq_num");
		my $asym_id		= $parent_node->findvalue('@asym_id');
		my $pdb_ins_code;
		if ($parent_node->findvalue('PDBx:pdb_ins_code/@xsi:nil') ne 'true') {
			$pdb_ins_code = $parent_node->findvalue('PDBx:pdb_ins_code');
		}
		if ($DEBUG > 1) { 
			print "DEBUG: $seq_id $pdb_seq_num $auth_seq_num\n";
		}
		if (defined($pdb_seq_num)) { 
			if ($pdb_ins_code) { 
				$pdb_seq_nums[$seq_id] = $pdb_seq_num . $pdb_ins_code;
			}else{
				$pdb_seq_nums[$seq_id] 	= $pdb_seq_num;
			}
		}

		if ($auth_seq_num =~ /\d+/) { 
			push (@struct_seq_ids, $seq_id);
		}
		push (@seq_ids, $seq_id);
		$target_asym_id{$asym_id}++;
	}
	if (scalar(keys %target_asym_id) > 1) { 
		die "$sub: More than one asym_id for chain $chain, aborting...\n" 
	}else{
		$target_asym_id = join(' ', keys %target_asym_id);
	}

	return (\@seq_ids, \@struct_seq_ids, \@pdb_seq_nums, $target_asym_id);
}
sub pdbml2pdb_cartn { 
	my $sub = 'pdbml2pdb_cartn';
	my ($fn, $cartn_aref, $mode, $model) = @_;

	if(!$model) { $model = 1 } 
	if ($DEBUG) { 
		print "DEBUG $sub: $fn $cartn_aref $mode\n";
	}

	if ($mode !~ /pdb|seq|inx/) { 
		die "ERROR! $sub: $mode not recognized (try pdb or seq)\n";
	}

	if (scalar(@$cartn_aref) == 0) { 
		print "WARNING $sub: passed empty coord ref for $fn\n";
		return 0
	}

	$fn =~ /\w+/ or die "Weird filename $fn, aborting...\n";
	
	open (PDB, ">$fn") or die "Could not open $fn for writing\n";

	my $atom_counter = 1;
	my $residue_counter = 1;
	my $last_seqid;
	my ($record, $atom_number, $alt_id, $atom, $residue, $chain, $residue_number, $ins_code, $cartn_x, $cartn_y, $cartn_z, $occupancy, $bfactor);

format PDB = 
@<<<<<@>>>> @<<<@@|| @@>>>@   @###.###@###.###@###.###@##.##@##.##
$record, $atom_number, $atom, $alt_id, $residue, $chain, $residue_number, $ins_code, $cartn_x, $cartn_y, $cartn_z, $occupancy, $bfactor
.	

		for (my $i = 0; $i < scalar(@$cartn_aref); $i++) { 
		 $record 	= 'ATOM';
		 $atom_number 	= $atom_counter;
		 $atom	= $$cartn_aref[$i]{label_atom_id};
		 if (length($atom) < 4) { 
			 $atom = ' ' . $atom;
		 }
		 $residue	= $$cartn_aref[$i]{label_comp_id};

		 $chain		= 'A';
		 $alt_id = $$cartn_aref[$i]{alt_id};
		 #if ($alt_id =~ /\w/) { print "DEBUG alt id $alt_id\n"}
		 if ($mode eq 'seq') { 
			 $residue_number = $$cartn_aref[$i]{label_seq_id};
			 $ins_code = ' ';
		 }elsif($mode eq 'pdb') { 
			 $residue_number = $$cartn_aref[$i]{auth_seq_id};
			 if ($$cartn_aref[$i]{ins_code}) { 
				 $ins_code = $$cartn_aref[$i]{ins_code};
			 }else{
				 $ins_code = ' '; 
			 }
		 }elsif($mode eq 'inx') { 
			 if (!$last_seqid) { 
				 $residue_number = $residue_counter;
				 $last_seqid = $$cartn_aref[$i]{label_seq_id};
				 $residue_counter++;
			 }elsif ($last_seqid != $$cartn_aref[$i]{label_seq_id}) {
				 $last_seqid = $$cartn_aref[$i]{label_seq_id};
				 $residue_number = $residue_counter;
				 $residue_counter++;
			 }
			$ins_code = ' ';
		 }
			 

		 $cartn_x	= $$cartn_aref[$i]{Cartn_x};
		 $cartn_y	= $$cartn_aref[$i]{Cartn_y};
		 $cartn_z	= $$cartn_aref[$i]{Cartn_z};

		 $bfactor	= $$cartn_aref[$i]{bfactor};
		 $occupancy	= $$cartn_aref[$i]{occupancy};
		write PDB;
		$atom_counter++;
	}
	print PDB "TER\n";
	close PDB;
	return 1;
	
}
sub pdbml_coord_ligand_parse_pull { 
	my $sub = (caller(0))[3];

	my ($pdb4, $non_poly_chem_comp_ids_href) = @_;

	$pdb4 = lc($pdb4);
	my $two;
	if ($pdb4 =~ /\w(\w{2})\w/) { 
		$two = $1;
	}else{
		die "ERROR! $sub: REGEXP fail on $pdb4\n";
	}

	my $xml_fh;
	my $pdb_fn = "$pdb_atom_top_dir/$two/$pdb4.xml.gz";
	if (-f $pdb_fn) { 
		if (!open ($xml_fh, "gzip -dc $pdb_fn |") ) { 
			warn "ERROR! $sub: Could not open $pdb_fn for reading:$!\n";
			return 0;
		}
	}else{
		warn "ERROR! $sub: Could not find $pdb_fn\n";
		return 0;
	}

	my $pdbml_pull	= new XML::LibXML::Reader(IO => $xml_fh) or
		die "ERROR $sub: Cannot read $pdb_atom_top_dir/$two/$pdb4.xml\n";

	my @atom_site;
	my $default_model_num = 1;

	my $i = 0;
	my $start = 0;
	$pdbml_pull->nextElement("PDBx:atom_site");	
	
	while (($i == 0 && $start == 0) || $pdbml_pull->nextSibling()) {
		$start++;

		my $depth 	= $pdbml_pull->depth;
		my $nodeType	= $pdbml_pull->nodeType;
		my $name	= $pdbml_pull->name;
		my $isempty	= $pdbml_pull->isEmptyElement;
		if ($nodeType == 14) { #NODE_TYPE 14 SIGNIFICANT WHITESPACE
			$start = 1;
			next;
		}
		my $text 	= 'NA';
		if ($pdbml_pull->hasValue) { 
			$text = $pdbml_pull->value;
		}

		if ($DEBUG > 2) { 
			print "DEBUG $sub: $depth $nodeType $name $isempty $text\n";
		}

		my $node = $pdbml_pull->copyCurrentNode(1);

		my $atom_site_id	= $node->findvalue('@id');

		if ($DEBUG > 2) { 
			print "DEBUG $sub: atom_site_id $atom_site_id\n";
		}

		my $label_comp_id	= $node->findvalue('PDBx:label_comp_id');
		if (!$non_poly_chem_comp_ids_href || $$non_poly_chem_comp_ids_href{$label_comp_id}) { 
			my $B_iso_or_equiv	= $node->findvalue('PDBx:B_iso_or_equiv');

			#Coords
			my $Cartn_x	= $node->findvalue('PDBx:Cartn_x');
			my $Cartn_y	= $node->findvalue('PDBx:Cartn_y');
			my $Cartn_z	= $node->findvalue('PDBx:Cartn_z');

			#Residue indices
			my $label_asym_id	= $node->findvalue('PDBx:label_asym_id');
			my $label_atom_id	= $node->findvalue('PDBx:label_atom_id');
			my $label_seq_id	= $node->findvalue('PDBx:label_seq_id');

			my $auth_seq_id		= $node->findvalue('PDBx:auth_seq_id'); # Not sure this is necessarily the same as pdb seq num :(
			my $auth_asym_id	= $node->findvalue('PDBx:auth_asym_id');

			#print "??$label_comp_id $label_asym_id $label_seq_id $auth_asym_id $auth_seq_id\n";

			#occupancy
			my $occupancy	= $node->findvalue('PDBx:occupancy');

			#B-Factor
			my $bfactor	= $node->findvalue('PDBx:B_iso_or_equiv');

			#insertion code
			my $pdbx_PDB_ins_code;
			if ($node->exists('PDBx:pdbx_PDB_ins_code')) { 
				$pdbx_PDB_ins_code = $node->findvalue('PDBx:pdbx_PDB_ins_code');
			}
			#alt_id 
			my $label_alt_id;
			if ($node->exists('PDBx:label_alt_id')) { 
				$label_alt_id	= $node->findvalue('PDBx:label_alt_id');
			}

			if ($DEBUG > 2) { 
				print "DEBUG $sub: label_asym_id $label_asym_id\n";
			}
			#model num
			my $pdbx_PDB_model_num = $node->findvalue('PDBx:pdbx_PDB_model_num');
			if ($pdbx_PDB_model_num && $pdbx_PDB_model_num != $default_model_num) { $start = 1; next } 

			$atom_site[$i]{atom_site_id}	= $atom_site_id;

			$atom_site[$i]{Cartn_x}	= $Cartn_x;
			$atom_site[$i]{Cartn_y}	= $Cartn_y;
			$atom_site[$i]{Cartn_z}	= $Cartn_z;

			$atom_site[$i]{occupancy}	= $occupancy;
			$atom_site[$i]{bfactor}		= $bfactor;

			$atom_site[$i]{label_asym_id}	= $label_asym_id;
			$atom_site[$i]{label_atom_id}	= $label_atom_id;
			$atom_site[$i]{label_comp_id}	= $label_comp_id;
			$atom_site[$i]{label_seq_id}	= $label_seq_id;

			$atom_site[$i]{auth_seq_id}	= $auth_seq_id;
			$atom_site[$i]{auth_asym_id}	= $auth_asym_id;

			$atom_site[$i]{ins_code}	= $pdbx_PDB_ins_code;
			$atom_site[$i]{alt_id}		= $label_alt_id;
			$atom_site[$i]{model_num}	= $pdbx_PDB_model_num;
			$i++;
		}
	}
	return \@atom_site;
}

sub pdbml_coord_parse_pull { 
	my $sub = 'pdbml_coord_parse_pull';

	my ($pdb4, $seq_id_aref, $target_asym_id) = @_;

	if ($DEBUG) { 
		print "DEBUG $sub: pdb4=>$pdb4 asym=>$target_asym_id\n";
	}


	$pdb4 = lc($pdb4);
	my $two;
	if ($pdb4 =~ /\w(\w{2})\w/) { 
		$two = $1;
	}else{
		die "ERROR! $sub: REGEXP fail on $pdb4\n";
	}

	my $xml_fh;
	if (-f "$pdb_atom_top_dir/$two/$pdb4.xml") { 
		if(!open ($xml_fh, "$pdb_atom_top_dir/$two/$pdb4.xml")) { 
			my $i = 0;
			while (!open($xml_fh, "$pdb_atom_top_dir/$two/pdb4.xml") && $i < $MAX_FILE_ATTEMPTS) { 
				sleep(1);
				$i++
			}
			if ($i == $MAX_FILE_ATTEMPTS) { 
				die "ERROR $sub: Could not open $pdb_atom_top_dir/$two/$pdb4.xml for reading:$!\n";
			}
		}
	}elsif(-f  "$pdb_atom_top_dir/$two/$pdb4.xml.gz") { 
		if (!open ($xml_fh, "gzip -dc $pdb_atom_top_dir/$two/$pdb4.xml.gz |") ) { 
			my $i = 0;
			while (!open ($xml_fh, "gzip -dc $pdb_atom_top_dir/$two/$pdb4.xml.gz |") && $i < $MAX_FILE_ATTEMPTS ) { 
				sleep(1);
				$i++;
			}
			if ($i == $MAX_FILE_ATTEMPTS) {
				die "ERROR $sub: Could not open $pdb_atom_top_dir/$two/$pdb4.xml.gz for reading:$!\n";
			}
		}
	}elsif(-f  "$pdb_obs_top_dir/$two/$pdb4.xml.gz") { 
		if (!open ($xml_fh, "gzip -dc $pdb_obs_top_dir/$two/$pdb4.xml.gz |") ) { 
			my $i = 0;
			while (!open ($xml_fh, "gzip -dc $pdb_obs_top_dir/$two/$pdb4.xml.gz |") && $i < $MAX_FILE_ATTEMPTS ) { 
				sleep(1);
				$i++;
			}
			if ($i == $MAX_FILE_ATTEMPTS) {
				die "ERROR $sub: Could not open $pdb_obs_top_dir/$two/$pdb4.xml.gz for reading:$!\n";
			}
		}
	}elsif(-f "$pdb_obs_top_dir/$two/$pdb4.xml.gz") { 
		if (!open ($xml_fh, "gzip-dc $pdb_obs_top_dir/$two/$pdb4.xml.gz|")) { 
			my $i = 0;
			while (!open($xml_fh, "gzip -dc $pdb_obs_top_dir/$two/$pdb4.xml.gz|") && $i < $MAX_FILE_ATTEMPTS) { 
				sleep(1);
				$i++;
			}
			if ($i == $MAX_FILE_ATTEMPTS) { 
				die "ERROR $sub: Could not open $pdb_obs_top_dir/$two/$pdb4.xml.gz for reading:$!\n";
			}
		}
	}else{
		print "ERROR $sub: PDBml file not found for $pdb4 in $pdb_atom_top_dir, aborting...\n";
		return 0;
	}

	my $pdbml_pull	= new XML::LibXML::Reader(IO => $xml_fh) or
		die "ERROR $sub: Cannot read $pdb_atom_top_dir/$two/$pdb4.xml\n";

	my @atom_site;
	my $default_model_num = 1;

	my $i = 0;
	my $start = 0;
	$pdbml_pull->nextElement("PDBx:atom_site");	
	
	while (($i == 0 && $start == 0) || $pdbml_pull->nextSibling()) {

		my $depth 	= $pdbml_pull->depth;
		my $nodeType	= $pdbml_pull->nodeType;
		my $name	= $pdbml_pull->name;
		my $isempty	= $pdbml_pull->isEmptyElement;
		if ($nodeType == 14) { #NODE_TYPE 14 SIGNIFICANT WHITESPACE
			$start = 1;
			next;
		}
		my $text 	= 'NA';
		if ($pdbml_pull->hasValue) { 
			$text = $pdbml_pull->value;
		}

		if ($DEBUG > 3) { 
			print "DEBUG $sub: $depth $nodeType $name $isempty $text\n";
		}

		my $node = $pdbml_pull->copyCurrentNode(1);

		my $atom_site_id	= $node->findvalue('@id');

		if ($DEBUG > 2) { 
			print "DEBUG $sub: atom_site_id $atom_site_id\n";
		}

		my $B_iso_or_equiv	= $node->findvalue('PDBx:B_iso_or_equiv');

		#Coords
		my $Cartn_x	= $node->findvalue('PDBx:Cartn_x');
		my $Cartn_y	= $node->findvalue('PDBx:Cartn_y');
		my $Cartn_z	= $node->findvalue('PDBx:Cartn_z');

		#Residue indices
		my $label_asym_id	= $node->findvalue('PDBx:label_asym_id');
		if ($label_asym_id ne $target_asym_id) { $start = 1;  next } 
		if ($DEBUG) { 
			print "DEBUG $sub: atom_site id $atom_site_id label_asym_id $label_asym_id\n";
		}
		my $label_atom_id	= $node->findvalue('PDBx:label_atom_id');
		my $label_comp_id	= $node->findvalue('PDBx:label_comp_id');
		my $label_seq_id	= $node->findvalue('PDBx:label_seq_id');

		my $auth_seq_id		= $node->findvalue('PDBx:auth_seq_id'); # Not sure this is necessarily the same as pdb seq num :(
		my $auth_asym_id	= $node->findvalue('PDBx:auth_asym_id');

		#occupancy
		my $occupancy	= $node->findvalue('PDBx:occupancy');

		#B-Factor
		my $bfactor	= $node->findvalue('PDBx:B_iso_or_equiv');

		#insertion code
		my $pdbx_PDB_ins_code;
		if ($node->exists('PDBx:pdbx_PDB_ins_code')) { 
			$pdbx_PDB_ins_code = $node->findvalue('PDBx:pdbx_PDB_ins_code');
		}
		#alt_id 
		my $label_alt_id;
		if ($node->exists('PDBx:label_alt_id')) { 
			$label_alt_id	= $node->findvalue('PDBx:label_alt_id');
		}

		if ($DEBUG > 2) { 
			print "DEBUG $sub: label_asym_id $label_asym_id\n";
		}
		#model num
		my $pdbx_PDB_model_num = 'Na';
		if ($node->exists('PDBx:pdbx_PDB_model_num')) { 
			$pdbx_PDB_model_num = $node->findvalue('PDBx:pdbx_PDB_model_num');
			if ($pdbx_PDB_model_num && $pdbx_PDB_model_num != $default_model_num) { $start = 1; next } 
		}

		$atom_site[$i]{atom_site_id}	= $atom_site_id;

		$atom_site[$i]{Cartn_x}	= $Cartn_x;
		$atom_site[$i]{Cartn_y}	= $Cartn_y;
		$atom_site[$i]{Cartn_z}	= $Cartn_z;

		$atom_site[$i]{occupancy}	= $occupancy;
		$atom_site[$i]{bfactor}		= $bfactor;

		$atom_site[$i]{label_asym_id}	= $label_asym_id;
		$atom_site[$i]{label_atom_id}	= $label_atom_id;
		$atom_site[$i]{label_comp_id}	= $label_comp_id;
		$atom_site[$i]{label_seq_id}	= $label_seq_id;

		$atom_site[$i]{auth_seq_id}	= $auth_seq_id;
		$atom_site[$i]{auth_asym_id}	= $auth_asym_id;

		$atom_site[$i]{ins_code}	= $pdbx_PDB_ins_code;
		$atom_site[$i]{alt_id}		= $label_alt_id;
		$atom_site[$i]{model_num}	= $pdbx_PDB_model_num;
		$i++;
	}
	
	if ($DEBUG) { 
		print "DEBUG $sub: " . scalar(@atom_site) . "\n";
	}

	my @output_coords;
	my $j = 0;

	for (my $i = 0; $i < scalar(@atom_site); $i++) { 
		my $atom_site_id	= $atom_site[$i]{atom_site_id};
		if ($atom_site[$i]{label_asym_id} eq $target_asym_id) { 
			if (grep {$_ eq $atom_site[$i]{label_seq_id}} @$seq_id_aref) { 
				if ($DEBUG > 1) { 
					print "DEBUG $sub: label_seq_id $atom_site[$i]{label_seq_id} label_asym_id $atom_site[$i]{label_asym_id} target_asym_id $target_asym_id ins_code $atom_site[$i]{ins_code}\n";
				}

				$output_coords[$j]{atom_site_id}	= $atom_site_id;
				$output_coords[$j]{label_atom_id}	= $atom_site[$i]{label_atom_id};
				$output_coords[$j]{label_comp_id}	= $atom_site[$i]{label_comp_id};
				$output_coords[$j]{label_asym_id}	= $atom_site[$i]{label_asym_id};
				$output_coords[$j]{label_seq_id}	= $atom_site[$i]{label_seq_id};
				$output_coords[$j]{auth_seq_id}		= $atom_site[$i]{auth_seq_id};
				$output_coords[$j]{ins_code}		= $atom_site[$i]{ins_code};
				$output_coords[$j]{alt_id}		= $atom_site[$i]{alt_id};

				$output_coords[$j]{bfactor}	= $atom_site[$i]{bfactor};
				$output_coords[$j]{occupancy}	= $atom_site[$i]{occupancy};
				$output_coords[$j]{Cartn_x}	= $atom_site[$i]{Cartn_x};
				$output_coords[$j]{Cartn_y}	= $atom_site[$i]{Cartn_y};
				$output_coords[$j]{Cartn_z}	= $atom_site[$i]{Cartn_z};

				$j++;
			}
		}
	}
	return \@output_coords;
}

sub pdbml_mc_coord_parse_pull { 
	my $sub = 'pdbml_mc_coord_parse_pull';

	my ($pdb4, $seq_id_aref, $chain_aref, $chain2asym) = @_;

	if ($DEBUG) { 
		print "DEBUG $sub: pdb4=>$pdb4 \n";
	}

	my %asym_lookup;
	if(scalar(@$seq_id_aref) != scalar(@$chain_aref)) { die "$sub: chain seqid array mismatch\n";}
	for (my $i = 0; $i < scalar(@$seq_id_aref); $i++) { 
		my $chain	= $$chain_aref[$i];
		my $seqid	= $$seq_id_aref[$i];
		my $asym_id	= $$chain2asym{$chain}{$seqid};
		$asym_lookup{$asym_id}{$seqid}++;
	}

	$pdb4 = lc($pdb4);
	my $two;
	if ($pdb4 =~ /\w(\w{2})\w/) { 
		$two = $1;
	}else{
		die "ERROR! $sub: REGEXP fail on $pdb4\n";
	}

	my $xml_fh;
	if (-f "$pdb_atom_top_dir/$two/$pdb4.xml") { 
		if(!open ($xml_fh, "$pdb_atom_top_dir/$two/$pdb4.xml")) { 
			my $i = 0;
			while (!open($xml_fh, "$pdb_atom_top_dir/$two/pdb4.xml") && $i < $MAX_FILE_ATTEMPTS) { 
				sleep(1);
				$i++
			}
			if ($i == $MAX_FILE_ATTEMPTS) { 
				die "ERROR $sub: Could not open $pdb_atom_top_dir/$two/$pdb4.xml for reading:$!\n";
			}
		}
	}elsif(-f  "$pdb_atom_top_dir/$two/$pdb4.xml.gz") { 
		if (!open ($xml_fh, "gzip -dc $pdb_atom_top_dir/$two/$pdb4.xml.gz |") ) { 
			my $i = 0;
			while (!open ($xml_fh, "gzip -dc $pdb_atom_top_dir/$two/$pdb4.xml.gz |") && $i < $MAX_FILE_ATTEMPTS ) { 
				sleep(1);
				$i++;
			}
			if ($i == $MAX_FILE_ATTEMPTS) {
				die "ERROR $sub: Could not open $pdb_atom_top_dir/$two/$pdb4.xml.gz for reading:$!\n";
			}
		}
	}elsif(-f "$pdb_obs_top_dir/$two/$pdb4.xml.gz") { 
		if (!open ($xml_fh, "gzip -dc $pdb_obs_top_dir/$two/$pdb4.xml.gz |") ) { 
			my $i = 0;
			while (!open ($xml_fh, "gzip -dc $pdb_obs_top_dir/$two/$pdb4.xml.gz |") && $i < $MAX_FILE_ATTEMPTS ) { 
				sleep(1);
				$i++;
			}
			if ($i == $MAX_FILE_ATTEMPTS) {
				die "ERROR $sub: Could not open $pdb_obs_top_dir/$two/$pdb4.xml.gz for reading:$!\n";
			}
		}		
	}else{
		print "ERROR $sub: PDBml file not found for $pdb4 in $pdb_atom_top_dir, aborting...\n";
		return 0;
	}

	my $pdbml_pull	= new XML::LibXML::Reader(IO => $xml_fh) or
		die "ERROR $sub: Cannot read $pdb_atom_top_dir/$two/$pdb4.xml\n";

	my @atom_site;
	my $default_model_num = 1;

	my $i = 0;
	my $start = 0;
	$pdbml_pull->nextElement("PDBx:atom_site");	
	
	while (($i == 0 && $start == 0) || $pdbml_pull->nextSibling()) {

		my $depth 	= $pdbml_pull->depth;
		my $nodeType	= $pdbml_pull->nodeType;
		my $name	= $pdbml_pull->name;
		my $isempty	= $pdbml_pull->isEmptyElement;
		if ($nodeType == 14) { #NODE_TYPE 14 SIGNIFICANT WHITESPACE
			$start = 1;
			next;
		}
		my $text 	= 'NA';
		if ($pdbml_pull->hasValue) { 
			$text = $pdbml_pull->value;
		}

		if ($DEBUG > 2) { 
			print "DEBUG $sub: $depth $nodeType $name $isempty $text\n";
		}

		my $node = $pdbml_pull->copyCurrentNode(1);

		my $atom_site_id	= $node->findvalue('@id');

		if ($DEBUG > 2) { 
			print "DEBUG $sub: atom_site_id $atom_site_id\n";
		}

		my $B_iso_or_equiv	= $node->findvalue('PDBx:B_iso_or_equiv');

		#Coords
		my $Cartn_x	= $node->findvalue('PDBx:Cartn_x');
		my $Cartn_y	= $node->findvalue('PDBx:Cartn_y');
		my $Cartn_z	= $node->findvalue('PDBx:Cartn_z');

		#Residue indices
		my $label_asym_id	= $node->findvalue('PDBx:label_asym_id');
		if (! $asym_lookup{$label_asym_id} ) { $start = 1; next } 
		my $label_atom_id	= $node->findvalue('PDBx:label_atom_id');
		my $label_comp_id	= $node->findvalue('PDBx:label_comp_id');
		my $label_seq_id	= $node->findvalue('PDBx:label_seq_id');

		my $auth_seq_id		= $node->findvalue('PDBx:auth_seq_id'); # Not sure this is necessarily the same as pdb seq num :(

		#occupancy
		my $occupancy	= $node->findvalue('PDBx:occupancy');

		#B-Factor
		my $bfactor	= $node->findvalue('PDBx:B_iso_or_equiv');

		#insertion code
		my $pdbx_PDB_ins_code;
		if ($node->exists('PDBx:pdbx_PDB_ins_code')) { 
			$pdbx_PDB_ins_code = $node->findvalue('PDBx:pdbx_PDB_ins_code');
		}
		#alt_id 
		my $label_alt_id;
		if ($node->exists('PDBx:label_alt_id')) { 
			$label_alt_id	= $node->findvalue('PDBx:label_alt_id');
		}

		if ($DEBUG > 2) { 
			print "DEBUG $sub: label_asym_id $label_asym_id\n";
		}
		#model num
		my $pdbx_PDB_model_num = $node->findvalue('PDBx:pdbx_PDB_model_num');
		if ($pdbx_PDB_model_num != $default_model_num) { $start = 1; next } 

		$atom_site[$i]{atom_site_id}	= $atom_site_id;

		$atom_site[$i]{Cartn_x}	= $Cartn_x;
		$atom_site[$i]{Cartn_y}	= $Cartn_y;
		$atom_site[$i]{Cartn_z}	= $Cartn_z;

		$atom_site[$i]{occupancy}	= $occupancy;
		$atom_site[$i]{bfactor}		= $bfactor;

		$atom_site[$i]{label_asym_id}	= $label_asym_id;
		$atom_site[$i]{label_atom_id}	= $label_atom_id;
		$atom_site[$i]{label_comp_id}	= $label_comp_id;
		$atom_site[$i]{label_seq_id}	= $label_seq_id;

		$atom_site[$i]{auth_seq_id}	= $auth_seq_id;

		$atom_site[$i]{ins_code}	= $pdbx_PDB_ins_code;
		$atom_site[$i]{alt_id}		= $label_alt_id;
		$atom_site[$i]{model_num}	= $pdbx_PDB_model_num;
		$i++;
	}
	
	if ($DEBUG) { 
		print "DEBUG $sub: " . scalar(@atom_site) . "\n";
	}

	my @output_coords;
	my $j = 0;

	for (my $i = 0; $i < scalar(@atom_site); $i++) { 
		my $atom_site_id	= $atom_site[$i]{atom_site_id};
		if ($asym_lookup{ $atom_site[$i]{label_asym_id}}{ $atom_site[$i]{label_seq_id} }) { 
			if ($DEBUG > 1) { 
				print "DEBUG $sub: label_seq_id $atom_site[$i]{label_seq_id} label_asym_id $atom_site[$i]{label_asym_id} target_asym_id $atom_site[$i]{label_asym_id} ins_code $atom_site[$i]{ins_code}\n";
			}

			$output_coords[$j]{atom_site_id}	= $atom_site_id;
			$output_coords[$j]{label_atom_id}	= $atom_site[$i]{label_atom_id};
			$output_coords[$j]{label_comp_id}	= $atom_site[$i]{label_comp_id};
			$output_coords[$j]{label_asym_id}	= $atom_site[$i]{label_asym_id};
			$output_coords[$j]{label_seq_id}	= $atom_site[$i]{label_seq_id};
			$output_coords[$j]{auth_seq_id}		= $atom_site[$i]{auth_seq_id};
			$output_coords[$j]{ins_code}		= $atom_site[$i]{ins_code};
			$output_coords[$j]{alt_id}		= $atom_site[$i]{alt_id};

			$output_coords[$j]{bfactor}	= $atom_site[$i]{bfactor};
			$output_coords[$j]{occupancy}	= $atom_site[$i]{occupancy};
			$output_coords[$j]{Cartn_x}	= $atom_site[$i]{Cartn_x};
			$output_coords[$j]{Cartn_y}	= $atom_site[$i]{Cartn_y};
			$output_coords[$j]{Cartn_z}	= $atom_site[$i]{Cartn_z};

			$j++;
		}
	}
	return \@output_coords;
}

sub pdbml_coord_parse { 
	my $sub = 'pdbml_coord_parse';
	my ($pdb4, $seq_id_aref, $target_asym_id) = @_;
	if ($DEBUG) { 
		print "DEBUG: $sub pdb4:$pdb4 asym:$target_asym_id\n";
	}
	my $two;
	if ($pdb4 =~ /\w(\w{2})\w/) { 
		$two = $1;
	}else{ 
		die "WARNING! $sub: REGEXP fail on $pdb4\n";
	}

	my $xml_fh;
	if (-f "$pdb_atom_top_dir/$two/$pdb4.xml") { 
		open ($xml_fh, "$pdb_atom_top_dir/$two/$pdb4.xml") or 
			die "$sub: Could not open $pdb_atom_top_dir/$two/$pdb4.xml for reading:$!\n";
	}elsif(-f "$pdb_atom_top_dir/$two/$pdb4.xml.gz") { 
		open ($xml_fh, "gzip -dc $pdb_atom_top_dir/$two/$pdb4.xml.gz |") or
			die "$sub: Could not open gzip pipe $pdb_atom_top_dir/$two/$pdb4-atom.xml.gz for reading:$!\n";
	}elsif(-f  "$pdb_obs_top_dir/$two/$pdb4.xml.gz") { 
		if (!open ($xml_fh, "gzip -dc $pdb_obs_top_dir/$two/$pdb4.xml.gz |") ) { 
			my $i = 0;
			while (!open ($xml_fh, "gzip -dc $pdb_obs_top_dir/$two/$pdb4.xml.gz |") && $i < $MAX_FILE_ATTEMPTS ) { 
				sleep(1);
				$i++;
			}
			if ($i == $MAX_FILE_ATTEMPTS) {
				die "ERROR $sub: Could not open $pdb_obs_top_dir/$two/$pdb4.xml.gz for reading:$!\n";
			}
		}
	}else{

		print "PDBml file not found for $pdb4 in $pdb_atom_top_dir, aborting...\n";
		return 0;
	}

	my $pdbml = XML::LibXML->load_xml (
		IO => $xml_fh);
	close $xml_fh;

	#use model 1 by default
	my $default_model_num = '1';

	my $atom_siteCategoryXPath = "//PDBx:atom_siteCategory/PDBx:atom_site/PDBx:pdbx_PDB_model_num[text() ='$default_model_num']";	

	my $atom_site_nodes = $pdbml->findnodes($atom_siteCategoryXPath);
	my @atom_site;
	if ($DEBUG){ 
		print "DEBUG: $sub atom_site($pdb4) " . $atom_site_nodes->size() . "\n";
	}

	foreach my $node ($atom_site_nodes->get_nodelist() ) { 

		my $parent_node 	= $node->parentNode;
	
		my $atom_site_id	= $parent_node->findvalue('@id');
		
		my $B_iso_or_equiv	= $parent_node->findvalue('PDBx:B_iso_or_equiv');

		#Coords
		my $Cartn_x		= $parent_node->findvalue('PDBx:Cartn_x');
		my $Cartn_y		= $parent_node->findvalue('PDBx:Cartn_y');
		my $Cartn_z		= $parent_node->findvalue('PDBx:Cartn_z');

		#Residue indices
		my $label_asym_id	= $parent_node->findvalue('PDBx:label_asym_id');
		my $label_atom_id	= $parent_node->findvalue('PDBx:label_atom_id');
		my $label_comp_id	= $parent_node->findvalue('PDBx:label_comp_id');
		my $label_seq_id	= $parent_node->findvalue('PDBx:label_seq_id');

		#occupancy
		my $occupancy		= $parent_node->findvalue('PDBx:occupancy');

		#insertion code
		my $pdbx_PDB_ins_code = 'nil';
		if ($parent_node->findvalue('PDBx:pdbx_PDB_ins_code/@xsi:nil') ne 'true') { 
			$pdbx_PDB_ins_code = $parent_node->findvalue('PDBx:pdbx_PDB_ins_code');
		}

		#model number (always '1' in the case of xray structures)
		my $pdbx_PDB_model_num = $parent_node->findvalue('PDBx:pdbx_PDB_model_num');
	
		$atom_site[$atom_site_id]{atom_site_id} = $atom_site_id;

		#Cords
		$atom_site[$atom_site_id]{Cartn_x} 	= $Cartn_x;
		$atom_site[$atom_site_id]{Cartn_y} 	= $Cartn_y;
		$atom_site[$atom_site_id]{Cartn_z} 	= $Cartn_z;
	
		#label identifiers are required to be sequential, positive ascending integers
		$atom_site[$atom_site_id]{label_asym_id} = $label_asym_id;
		$atom_site[$atom_site_id]{label_atom_id} = $label_atom_id;#->CCD.chem_comp.atom_id
		$atom_site[$atom_site_id]{label_comp_id} = $label_comp_id;#->CCD.chem_comp.comp_id
		$atom_site[$atom_site_id]{label_seq_id}  = $label_seq_id; #->entity_poly_seq.num 
		$atom_site[$atom_site_id]{occupancy} 	= $occupancy;

		$atom_site[$atom_site_id]{ins_code} 	= $pdbx_PDB_ins_code;
		$atom_site[$atom_site_id]{model_num} 	= $pdbx_PDB_model_num;
	}	
	if ($DEBUG) { 
		print "DEBUG: $sub: " . scalar(@atom_site) . "\n";
	}
	my @output_coords;
	my $j = 0;
	for (my $i = 1; $i < scalar(@atom_site); $i++) { 
		my $atom_site_id = $atom_site[$i]{atom_site_id};
		if ($atom_site[$atom_site_id]{label_asym_id} eq $target_asym_id) { 
			if (grep {$_ eq $atom_site[$atom_site_id]{label_seq_id}} @$seq_id_aref) { 	
				if ($DEBUG > 2) { 
					print "DEBUG $sub: label_seq_id $atom_site[$atom_site_id]{label_seq_id} label_asym_id $atom_site[$atom_site_id]{label_asym_id} target_asym_id $target_asym_id\n";
				}
			
				$output_coords[$j]{atom_site_id}	= $atom_site_id; 
				$output_coords[$j]{label_atom_id}	= $atom_site[$atom_site_id]{label_atom_id};
				$output_coords[$j]{label_comp_id}	= $atom_site[$atom_site_id]{label_comp_id};
				$output_coords[$j]{label_asym_id}	= $atom_site[$atom_site_id]{label_asym_id};
				$output_coords[$j]{label_seq_id}	= $atom_site[$atom_site_id]{label_seq_id};

				$output_coords[$j]{Cartn_x}		= $atom_site[$atom_site_id]{Cartn_x};
				$output_coords[$j]{Cartn_y}		= $atom_site[$atom_site_id]{Cartn_y};
				$output_coords[$j]{Cartn_z}		= $atom_site[$atom_site_id]{Cartn_z};
				$j++;
			}
		}
	}

	return \@output_coords;
}
sub pdbml_seq_fetch  { 
	my $sub = 'pdbml_seq_fetch';

	my ( $pdb4, $asym_id, $chain, $seqid_aref ) = @_;


	$pdb4 = lc($pdb4);
	my $two;
	if ($pdb4 =~ /\w(\w{2})\w/) { 
		$two = $1;
	}else{
		print "ERROR $sub: $pdb4 does not regexp\n";
		return 0;
	}

	if ($DEBUG) { 
		printf "DEBUG $sub: $pdb4 $chain $asym_id %s\n", rangify(@$seqid_aref);
	}

	if (!$seqid_aref || scalar(@$seqid_aref) == 0) { 
		my $scalar = scalar(@$seqid_aref);
		die "ERROR! $sub: Passed empty array reference, $pdb4, $asym_id,$chain, $scalar\n";
	}

	my $pdbml = pdbml_load($pdb4);

	#Entity id?
	my $entity_poly_XPath = qq{//PDBx:entity_polyCategory/PDBx:entity_poly};
	my $entity_id;
	if ($pdbml->exists($entity_poly_XPath)) { 
		my $entity_poly_nodes	= $pdbml->findnodes($entity_poly_XPath);
		foreach my $node ($entity_poly_nodes->get_nodelist() ) { 
			my $chain_str = $node->findvalue('PDBx:pdbx_strand_id');
			my @chains = split(/\,/, $chain_str);
			foreach my $local_chain (@chains) { 
				if ($local_chain eq $chain) { 
					$entity_id = $node->findvalue('@entity_id');
					last;
				}
			}
		}
	}else{
		die "ERROR $sub: chain $chain entity id not found in $pdb4\n";
	}

	#ADD MODRES
	my @modres;
	my %mod_seq_lookup;
	my $mod_i = 0;
	my $struct_mod_res_XPath = '//PDBx:pdbx_struct_mod_residueCategory//PDBx:pdbx_struct_mod_residue';
	foreach my $smr_node ($pdbml->findnodes($struct_mod_res_XPath)->get_nodelist() ) { 

		my $smr_id		= $smr_node->findvalue('@id');
		my $details		= $smr_node->findvalue('PDBx:details');
		my $label_asym_id	= $smr_node->findvalue('PDBx:label_asym_id');
		my $label_comp_id	= $smr_node->findvalue('PDBx:label_comp_id');
		my $label_seq_id	= $smr_node->findvalue('PDBx:label_seq_id');
		my $parent_comp_id	= $smr_node->findvalue('PDBx:parent_comp_id');

		if ($asym_id eq $label_asym_id) { 
			$modres[$mod_i]{smr_id}		= $smr_id;
			$modres[$mod_i]{details}	= $details;
			$modres[$mod_i]{label_asym_id}	= $label_asym_id;
			$modres[$mod_i]{label_comp_id}	= $label_comp_id;
			$modres[$mod_i]{label_seq_id}	= $label_seq_id;
			$modres[$mod_i]{parent_comp_id}	= $parent_comp_id;
			$mod_seq_lookup{$label_seq_id} = $mod_i;
			if ($DEBUG) { 
				print "DEBUG $sub MOD $mod_i $label_seq_id $label_comp_id $parent_comp_id\n";
			}
			$mod_i++;
		}
	}
	#This is an absolutely terrible way to do this, use pdbx_poly_seq_scheme NOT entity_poly_seq
	#my $entity_poly_seq_XPath = qq{//PDBx:entity_poly_seqCategory/PDBx:entity_poly_seq[\@entity_id='$entity_id']};
	#my $entity_poly_seq_nodes	= $pdbml->findnodes($entity_poly_seq_XPath);
	#if (!$pdbml->exists($entity_poly_seq_XPath)) { 
	#	die "ERROR! $sub: entity_poly_seq not found for entity id $entity_id for $pdb4\n";
	#}

	my $pdbx_poly_seq_schemeCategoryXPath = qq{//PDBx:pdbx_poly_seq_schemeCategory/PDBx:pdbx_poly_seq_scheme[\@entity_id='$entity_id'][\@asym_id='$asym_id']};
	my $pdbx_poly_seq_scheme_nodes = $pdbml->findnodes($pdbx_poly_seq_schemeCategoryXPath);


	my %three_to_one = (
				ALA 	=> 'A',
				CYS 	=> 'C',
				ASP	=> 'D',
				GLU	=> 'E',
				PHE	=> 'F',
				GLY	=> 'G',
				HIS	=> 'H',
				ILE	=> 'I',
				LYS	=> 'K',
				LEU	=> 'L',
				MET	=> 'M',
				ASN	=> 'N',
				PRO	=> 'P',
				GLN	=> 'Q',
				ARG	=> 'R',
				SER	=> 'S',
				THR	=> 'T',
				VAL	=> 'V',
				TRP	=> 'W',
				TYR	=> 'Y');


	my %seq;
	foreach my $node ($pdbx_poly_seq_scheme_nodes->get_nodelist() ) { 
		my $seq_id 	= $node->findvalue('@seq_id');
		my $mon_id	= $node->findvalue('@mon_id');

		my $auth_seq_num	= $node->findvalue("PDBx:auth_seq_num");

		if (grep {$seq_id == $_} @$seqid_aref) { 
			my $one_letter_aa;
			if ($three_to_one{$mon_id}) { 
				$one_letter_aa	= $three_to_one{$mon_id};
			}elsif(exists $mod_seq_lookup{$seq_id}) { 
				$mon_id = $modres[$mod_seq_lookup{$seq_id}]{parent_comp_id};
			}else{
				$one_letter_aa	= 'X';
			}

			if ($DEBUG) { 
				print "DEBUG $sub: $seq_id $mon_id $one_letter_aa\n";
				printf "???%i %i\n", scalar keys %seq, scalar(keys %seq);
			}

			$seq{$seq_id} = $mon_id
					
		}
	}
	printf "%i\n", scalar(keys %seq);

	if (scalar(keys %seq) == 0) { 
		print "WARNING! $sub: No sequence found for $pdb4, $chain, $asym_id\n";
		return 0;
	}

	#my $seq_str = join('', @seq);
	#return $seq_str;
	return \%seq;
}
sub pdbml_obs_fasta_fetch { 
	my $sub = 'pdbml_obs_fasta_fetch';

	my ( $pdb4, $asym_id, $chain, $seqid_aref ) = @_;


	if ($DEBUG) { 
		print "DEBUG $sub: $pdb4 $asym_id $chain " . scalar(@$seqid_aref) . "\n";
	}
	$pdb4 = lc($pdb4);
	#$chain = uc($chain);
	$chain = $chain;
	my $two;
	if ($pdb4 =~ /\w(\w{2})\w/) { 
		$two = $1;
	}else{
		print "ERROR $sub: $pdb4 does not regexp\n";
		return 0;
	}

	if ($DEBUG) { 
		print "DEBUG $sub: $pdb4 $chain\n";
	}

	if (!$seqid_aref || scalar(@$seqid_aref) == 0) { 
		my $scalar = scalar(@$seqid_aref);
		die "ERROR! $sub: Passed empty array reference, $pdb4, $asym_id,$chain, $scalar\n";
	}
	
	my $xml_fh;
	if (-f "$pdb_obs_top_dir/$two/$pdb4.xml") { 
		if(!open ($xml_fh, "$pdb_obs_top_dir/$two/$pdb4.xml")) { 
			my $i = 0;
			while (!open($xml_fh, "$pdb_obs_top_dir/$two/$pdb4.xml") && $i < $MAX_FILE_ATTEMPTS) { 
				sleep(1);
				$i++
			}
			if ($i == $MAX_FILE_ATTEMPTS) { 
				die "ERROR $sub: Could not open $pdb_obs_top_dir/$two/$pdb4.xml for reading:$!\n";
			}
		}
	}elsif(-f  "$pdb_obs_top_dir/$two/$pdb4.xml.gz") { 
		if (!open ($xml_fh, "gzip -dc $pdb_obs_top_dir/$two/$pdb4.xml.gz |") ) { 
			my $i = 0;
			while (!open ($xml_fh, "gzip -dc $pdb_obs_top_dir/$two/$pdb4.xml.gz |") && $i < $MAX_FILE_ATTEMPTS ) { 
				sleep(1);
				$i++;
			}
			if ($i == $MAX_FILE_ATTEMPTS) {
				die "ERROR $sub: Could not open $pdb_obs_top_dir/$two/$pdb4.xml.gz for reading:$!\n";
			}
		}
	}else{
		print "ERROR $sub: PDBml file not found for $pdb4 in $pdb_obs_top_dir, aborting...\n";
		return 0;
	}
	my $pdbml = XML::LibXML->load_xml (
		IO => $xml_fh);
	close $xml_fh;

	#Entity id?
	my $entity_poly_XPath = qq{//PDBx:entity_polyCategory/PDBx:entity_poly};
	my $entity_id;
	if ($pdbml->exists($entity_poly_XPath)) { 
		my $entity_poly_nodes	= $pdbml->findnodes($entity_poly_XPath);
		foreach my $node ($entity_poly_nodes->get_nodelist() ) { 
			my $chain_str = $node->findvalue('PDBx:pdbx_strand_id');
			if ($chain_str =~ /$chain/) { 
				$entity_id = $node->findvalue('@entity_id');
				last;
			}
		}
	}else{
		die "ERROR $sub: chain $chain entity id not found in $pdb4\n";
	}
	#print "CHAIN $chain ENTITY $entity_id\n";


	#This is an absolutely terrible way to do this, use pdbx_poly_seq_scheme NOT entity_poly_seq
	
	#my $entity_poly_seq_XPath = qq{//PDBx:entity_poly_seqCategory/PDBx:entity_poly_seq[\@entity_id='$entity_id']};
	#my $entity_poly_seq_nodes	= $pdbml->findnodes($entity_poly_seq_XPath);
	#if (!$pdbml->exists($entity_poly_seq_XPath)) { 
	#	die "ERROR! $sub: entity_poly_seq not found for entity id $entity_id for $pdb4\n";
	#}

	#ADD MODRES
	my @modres;
	my %mod_seq_lookup;
	my $mod_i = 0;
	my $struct_mod_res_XPath = '//PDBx:pdbx_struct_mod_residueCategory//PDBx:pdbx_struct_mod_residue';
	foreach my $smr_node ($pdbml->findnodes($struct_mod_res_XPath)->get_nodelist() ) { 

		my $smr_id		= $smr_node->findvalue('@id');
		my $details		= $smr_node->findvalue('PDBx:details');
		my $label_asym_id	= $smr_node->findvalue('PDBx:label_asym_id');
		my $label_comp_id	= $smr_node->findvalue('PDBx:label_comp_id');
		my $label_seq_id	= $smr_node->findvalue('PDBx:label_seq_id');
		my $parent_comp_id	= $smr_node->findvalue('PDBx:parent_comp_id');

		if ($asym_id eq $label_asym_id) { 
			$modres[$mod_i]{smr_id}		= $smr_id;
			$modres[$mod_i]{details}	= $details;
			$modres[$mod_i]{label_asym_id}	= $label_asym_id;
			$modres[$mod_i]{label_comp_id}	= $label_comp_id;
			$modres[$mod_i]{label_seq_id}	= $label_seq_id;
			$modres[$mod_i]{parent_comp_id}	= $parent_comp_id;
			$mod_seq_lookup{$label_seq_id} = $mod_i;
			if ($DEBUG) { 
				print "DEBUG $sub MOD $mod_i $label_seq_id $label_comp_id $parent_comp_id\n";
			}
			$mod_i++;
		}
	}

	my $pdbx_poly_seq_schemeCategoryXPath = qq{//PDBx:pdbx_poly_seq_schemeCategory/PDBx:pdbx_poly_seq_scheme[\@entity_id='$entity_id'][\@asym_id='$asym_id']};
	my $pdbx_poly_seq_scheme_nodes = $pdbml->findnodes($pdbx_poly_seq_schemeCategoryXPath);

	my %three_to_one = (
				ALA 	=> 'A',
				CYS 	=> 'C',
				ASP	=> 'D',
				GLU	=> 'E',
				PHE	=> 'F',
				GLY	=> 'G',
				HIS	=> 'H',
				ILE	=> 'I',
				LYS	=> 'K',
				LEU	=> 'L',
				MET	=> 'M',
				ASN	=> 'N',
				PRO	=> 'P',
				GLN	=> 'Q',
				ARG	=> 'R',
				SER	=> 'S',
				THR	=> 'T',
				VAL	=> 'V',
				TRP	=> 'W',
				TYR	=> 'Y');


	my @seq;
	my %seq_ids;

	foreach my $node ($pdbx_poly_seq_scheme_nodes->get_nodelist() ) { 
		my $seq_id 	= $node->findvalue('@seq_id');
		my $mon_id	= $node->findvalue('@mon_id');

		if ($seq_ids{$seq_id}) { 
			next;
		}else{
			$seq_ids{$seq_id}++;
		}

		my $auth_seq_num	= $node->findvalue("PDBx:auth_seq_num");

		if (grep {$seq_id == $_} @$seqid_aref) { 
			my $one_letter_aa;
			if ($three_to_one{$mon_id}) { 
				$one_letter_aa	= $three_to_one{$mon_id};
			}elsif(exists $mod_seq_lookup{$seq_id} ) { 
				$one_letter_aa = $three_to_one{$modres[$mod_seq_lookup{$seq_id}]{parent_comp_id}};
			}else{
				$one_letter_aa	= 'X';
				#die "ERROR! $sub: No one letter abbreviation for mon_id $mon_id\n";
			}

			if ($DEBUG) { 
				print "DEBUG $sub: $seq_id $mon_id $one_letter_aa\n";
			}

			push (@seq, $one_letter_aa);
					
					
		}
	}

	if (scalar(@seq) == 0) { 
		print "WARNING! $sub: No sequence found for $pdb4, $chain\n";
		return 0;
	}

	my $seq_str = join('', @seq);
	return $seq_str;
}
sub pdbml_fasta_fetch { 
	my $sub = 'pdbml_fasta_fetch';

	my ( $pdb4, $asym_id, $chain, $seqid_aref ) = @_;

	if ($DEBUG) { 
		print "DEBUG $sub: $pdb4 $asym_id $chain " . scalar(@$seqid_aref) . "\n";
	}
	$pdb4 = lc($pdb4);
	#$chain = uc($chain);
	#$chain = $chain;
	if ($DEBUG) { 
		print "DEBUG $sub: $pdb4 $chain\n";
	}

	if (!$seqid_aref || scalar(@$seqid_aref) == 0) { 
		my $scalar = scalar(@$seqid_aref);
		die "ERROR! $sub: Passed empty array reference, $pdb4, $asym_id,$chain, $scalar\n";
	}
	
	my $pdbml = pdbml_load($pdb4);

	#Entity id?
	my $entity_poly_XPath = qq{//PDBx:entity_polyCategory/PDBx:entity_poly};
	my $entity_id;
	if ($pdbml->exists($entity_poly_XPath)) { 
		my $entity_poly_nodes	= $pdbml->findnodes($entity_poly_XPath);
		foreach my $node ($entity_poly_nodes->get_nodelist() ) { 
			my $chain_str = $node->findvalue('PDBx:pdbx_strand_id');
			my @chains = split(/,/, $chain_str);
			foreach my $i_chain (@chains) { 
				if ($i_chain eq $chain) { 
					$entity_id = $node->findvalue('@entity_id');
					last;
				}
			}
		}
	}else{
		die "ERROR $sub: chain $chain entity id not found in $pdb4\n";
	}

	#print "e: $entity_id\n";

	#This is an absolutely terrible way to do this, use pdbx_poly_seq_scheme NOT entity_poly_seq
	
	#my $entity_poly_seq_XPath = qq{//PDBx:entity_poly_seqCategory/PDBx:entity_poly_seq[\@entity_id='$entity_id']};
	#my $entity_poly_seq_nodes	= $pdbml->findnodes($entity_poly_seq_XPath);
	#if (!$pdbml->exists($entity_poly_seq_XPath)) { 
	#	die "ERROR! $sub: entity_poly_seq not found for entity id $entity_id for $pdb4\n";
	#}

	#ADD MODRES
	my @modres;
	my %mod_seq_lookup;
	my $mod_i = 0;
	my $struct_mod_res_XPath = '//PDBx:pdbx_struct_mod_residueCategory//PDBx:pdbx_struct_mod_residue';
	foreach my $smr_node ($pdbml->findnodes($struct_mod_res_XPath)->get_nodelist() ) { 

		my $smr_id		= $smr_node->findvalue('@id');
		my $details		= $smr_node->findvalue('PDBx:details');
		my $label_asym_id	= $smr_node->findvalue('PDBx:label_asym_id');
		my $label_comp_id	= $smr_node->findvalue('PDBx:label_comp_id');
		my $label_seq_id	= $smr_node->findvalue('PDBx:label_seq_id');
		my $parent_comp_id	= $smr_node->findvalue('PDBx:parent_comp_id');

		if ($asym_id eq $label_asym_id) { 
			$modres[$mod_i]{smr_id}		= $smr_id;
			$modres[$mod_i]{details}	= $details;
			$modres[$mod_i]{label_asym_id}	= $label_asym_id;
			$modres[$mod_i]{label_comp_id}	= $label_comp_id;
			$modres[$mod_i]{label_seq_id}	= $label_seq_id;
			$modres[$mod_i]{parent_comp_id}	= $parent_comp_id;
			$mod_seq_lookup{$label_seq_id} = $mod_i;
			if ($DEBUG) { 
				print "DEBUG $sub MOD $mod_i $label_seq_id $label_comp_id $parent_comp_id\n";
			}
			$mod_i++;
		}
	}

	my ($pdbx_poly_seq_schemeCategoryXPath, $pdbx_poly_seq_scheme_nodes);
	if ($entity_id) { 
		$pdbx_poly_seq_schemeCategoryXPath = qq{//PDBx:pdbx_poly_seq_schemeCategory/PDBx:pdbx_poly_seq_scheme[\@entity_id='$entity_id'][\@asym_id='$asym_id']};
		$pdbx_poly_seq_scheme_nodes = $pdbml->findnodes($pdbx_poly_seq_schemeCategoryXPath);
	}else{
		$pdbx_poly_seq_schemeCategoryXPath = qq{//PDBx:pdbx_poly_seq_schemeCategory/PDBx:pdbx_poly_seq_scheme[\@asym_id='$asym_id']};
		$pdbx_poly_seq_scheme_nodes = $pdbml->findnodes($pdbx_poly_seq_schemeCategoryXPath);
	}

	my %three_to_one = (
				ALA 	=> 'A',
				CYS 	=> 'C',
				ASP	=> 'D',
				GLU	=> 'E',
				PHE	=> 'F',
				GLY	=> 'G',
				HIS	=> 'H',
				ILE	=> 'I',
				LYS	=> 'K',
				LEU	=> 'L',
				MET	=> 'M',
				ASN	=> 'N',
				PRO	=> 'P',
				GLN	=> 'Q',
				ARG	=> 'R',
				SER	=> 'S',
				THR	=> 'T',
				VAL	=> 'V',
				TRP	=> 'W',
				TYR	=> 'Y');


	my @seq;
	my %seq_ids;

	my %search_seqids;
	foreach my $seq_id (@$seqid_aref) { 
		$search_seqids{$seq_id}++;
	}

	my ($c1, $c2, $c3);
	$c1 = 0; $c2 = 0; $c3 = 0;
	foreach my $node ($pdbx_poly_seq_scheme_nodes->get_nodelist() ) { 
		my $strand_id	= $node->findvalue('PDBx:pdb_strand_id');
		if ($strand_id ne $chain) { 
			$c1++; 
			next;
		} 
		my $seq_id 	= $node->findvalue('@seq_id');
		my $mon_id	= $node->findvalue('@mon_id');

		if ($seq_ids{$seq_id}) { 
			$c2++;
			next;
		}else{
			$seq_ids{$seq_id}++;
		}

		my $auth_seq_num	= $node->findvalue("PDBx:auth_seq_num");

		if ($search_seqids{$seq_id}) { 
			my $one_letter_aa;
			if ($three_to_one{$mon_id}) { 
				$one_letter_aa	= $three_to_one{$mon_id};
			}elsif(exists $mod_seq_lookup{$seq_id} ) { 
				$one_letter_aa = $three_to_one{$modres[$mod_seq_lookup{$seq_id}]{parent_comp_id}};
			}else{
				$one_letter_aa	= 'X';
				#die "ERROR! $sub: No one letter abbreviation for mon_id $mon_id\n";
				warn "ERROR! $sub: No abbrev for mon_id $mon_id\n";
			}

			if ($DEBUG) { 
				print "DEBUG $sub: $seq_id $mon_id $one_letter_aa\n";
			}

			push (@seq, $one_letter_aa);
		}else{
			$c3++;
		}
	}

	if (scalar(@seq) == 0) { 
		print "WARNING! $sub: No sequence found for $pdb4, $chain ($c1, $c2, $c3)\n";
		return 0;
	}

	my $seq_str = join('', @seq);
	return $seq_str;
}
sub pdbml_obs_seq_parse { 
	my $sub = 'pdbml_obs_seq_parse';
	my ( $pdb4, $chain ) = @_;

	$pdb4 = lc($pdb4);
	my $two;
	if ($pdb4 =~ /\w(\w{2})\w/) { 
		$two = $1;
	}else{
		print "ERROR $sub: $pdb4 does not regexp\n";
		return 0;
	}

	if ($DEBUG) { 
		print "DEBUG $sub: $pdb4 $chain\n";
	}
	

	my $pdbml = pdbml_open($pdb4);

	my $pdbx_poly_seq_schemeCategoryXPath = "//PDBx:pdbx_poly_seq_schemeCategory/PDBx:pdbx_poly_seq_scheme/PDBx:pdb_strand_id[text()=\'$chain\']";
	my $pdbx_poly_seq_scheme_nodes = $pdbml->findnodes($pdbx_poly_seq_schemeCategoryXPath);

	my @seq_ids;
	my %seq_ids;
	my @struct_seq_ids;
	my @pdb_seq_nums;

	if ($DEBUG) { 
		printf "DEBUG: $pdb4 $chain seqid nodes %i\n", scalar($pdbx_poly_seq_scheme_nodes->get_nodelist());
	}

	my %target_asym_id;
	my $target_asym_id = 'NA';
	foreach my $node ($pdbx_poly_seq_scheme_nodes->get_nodelist()) { 
		my $parent_node	= $node->parentNode;
		my $auth_seq_num	= $parent_node->findvalue("PDBx:auth_seq_num");
		my $seq_id		= $parent_node->findvalue('@seq_id');
		my $pdb_seq_num 	= $parent_node->findvalue("PDBx:pdb_seq_num");
		my $asym_id		= $parent_node->findvalue('@asym_id');
		my $hetero		= $parent_node->findvalue('PDBx:hetero');
		if ($hetero eq 'y') { 
			print "WARNING! $sub: Heterogeneity at $seq_id\n";
		}
		if ($seq_ids{$seq_id}) { 
			print "WARNING! $sub: duplicate seq id $seq_id in asym_id $asym_id, hetero?\n";
			next;
		}else{
			$seq_ids{$seq_id}++;
		}

		my $pdb_ins_code;
		if ($parent_node->findvalue('PDBx:pdb_ins_code/@xsi:nil') ne 'true') {
			$pdb_ins_code = $parent_node->findvalue('PDBx:pdb_ins_code');
		}
		if ($DEBUG > 1) { 
			print "DEBUG: $seq_id $pdb_seq_num $auth_seq_num\n";
		}
		if (defined($pdb_seq_num)) { 
			if ($pdb_ins_code) { 
				$pdb_seq_nums[$seq_id] = $pdb_seq_num . $pdb_ins_code;
			}else{
				$pdb_seq_nums[$seq_id] 	= $pdb_seq_num;
			}
		}

		if ($auth_seq_num =~ /\d+/) { 
			push (@struct_seq_ids, $seq_id);
		}
		push (@seq_ids, $seq_id);
		$target_asym_id{$asym_id}++;
	}
	if (scalar(keys %target_asym_id) > 1) { 
		die "$sub: More than one asym_id for chain $chain, aborting...\n" 
	}else{
		$target_asym_id = join(' ', keys %target_asym_id);
	}

	return (\@seq_ids, \@struct_seq_ids, \@pdb_seq_nums, $target_asym_id);
}
sub pdbml_seq_parse { 
	my $sub = 'pdbml_seq_parse';
	my ( $pdb4, $chain ) = @_;

	$pdb4 = lc($pdb4);
	my $two;
	if ($pdb4 =~ /\w(\w{2})\w/) { 
		$two = $1;
	}else{
		print "ERROR $sub: $pdb4 does not regexp\n";
		return 0;
	}

	if ($DEBUG) { 
		print "DEBUG $sub: $pdb4 $chain\n";
	}
	
	my $pdbml = pdbml_load($pdb4);

	my $pdbx_poly_seq_schemeCategoryXPath = "//PDBx:pdbx_poly_seq_schemeCategory/PDBx:pdbx_poly_seq_scheme/PDBx:pdb_strand_id[text()=\'$chain\']";
	my $pdbx_poly_seq_scheme_nodes = $pdbml->findnodes($pdbx_poly_seq_schemeCategoryXPath);

	my @seq_ids;
	my %seq_ids;
	my @struct_seq_ids;
	my @pdb_seq_nums;

	if ($DEBUG) { 
		printf "DEBUG: $pdb4 $chain seqid nodes %i\n", scalar($pdbx_poly_seq_scheme_nodes->get_nodelist());
	}

	my %target_asym_id;
	my $target_asym_id = 'NA';
	foreach my $node ($pdbx_poly_seq_scheme_nodes->get_nodelist()) { 
		my $parent_node	= $node->parentNode;
		my $auth_seq_num	= $parent_node->findvalue("PDBx:auth_seq_num");
		my $seq_id		= $parent_node->findvalue('@seq_id');
		my $pdb_seq_num 	= $parent_node->findvalue("PDBx:pdb_seq_num");
		my $asym_id		= $parent_node->findvalue('@asym_id');
		my $hetero		= $parent_node->findvalue('PDBx:hetero');
		if ($hetero eq 'y') { 
			print "WARNING! $sub: Heterogeneity at $seq_id\n";
		}
		if ($seq_ids{$seq_id}) { 
			print "WARNING! $sub: duplicate seq id $seq_id in asym_id $asym_id, hetero?\n";
			next;
		}else{
			$seq_ids{$seq_id}++;
		}

		my $pdb_ins_code;
		if ($parent_node->findvalue('PDBx:pdb_ins_code/@xsi:nil') ne 'true') {
			$pdb_ins_code = $parent_node->findvalue('PDBx:pdb_ins_code');
		}
		if ($DEBUG > 1) { 
			print "DEBUG: $seq_id $pdb_seq_num $auth_seq_num\n";
		}
		if (defined($pdb_seq_num)) { 
			if ($pdb_ins_code) { 
				$pdb_seq_nums[$seq_id] = $pdb_seq_num . $pdb_ins_code;
			}else{
				$pdb_seq_nums[$seq_id] 	= $pdb_seq_num;
			}
		}

		if ($auth_seq_num =~ /\d+/) { 
			push (@struct_seq_ids, $seq_id);
		}
		push (@seq_ids, $seq_id);
		$target_asym_id{$asym_id}++;
	}
	if (scalar(keys %target_asym_id) > 1) { 
		die "$sub: More than one asym_id for chain $chain, aborting...\n" 
	}else{
		$target_asym_id = join(' ', keys %target_asym_id);
	}

	return (\@seq_ids, \@struct_seq_ids, \@pdb_seq_nums, $target_asym_id);
}
sub pdbml_mc_seq_parse { 
	my $sub = 'pdbml_mc_seq_parse';
	my ( $pdb4, $chain_aref ) = @_;

	$pdb4 = lc($pdb4);

	my $pdbml	= pdbml_load($pdb4);

	my @seq_ids;
	my @struct_seq_ids;
	#my @pdb_seq_nums;
	my %pdb_seq_nums;

	my %asym_seq_nums;

	my @chain_ids;
	my @struct_chain_ids; #added 2/26/2013

	my %target_asym_id;
	my $target_asym_id = 'NA';

	foreach my $chain (@$chain_aref) { 
		my $pdbx_poly_seq_schemeCategoryXPath = "//PDBx:pdbx_poly_seq_schemeCategory/PDBx:pdbx_poly_seq_scheme/PDBx:pdb_strand_id[text()=\'$chain\']";
		my $pdbx_poly_seq_scheme_nodes = $pdbml->findnodes($pdbx_poly_seq_schemeCategoryXPath);
		if ($DEBUG) { 
			printf "DEBUG: $pdb4 $chain seqid nodes %i\n", scalar($pdbx_poly_seq_scheme_nodes->get_nodelist());
		}
	
		foreach my $node ($pdbx_poly_seq_scheme_nodes->get_nodelist()) { 
			my $chain_id 		= $node->findvalue('text()');

			my $parent_node		= $node->parentNode;
			my $auth_seq_num	= $parent_node->findvalue("PDBx:auth_seq_num");
			my $seq_id		= $parent_node->findvalue('@seq_id');
			my $pdb_seq_num 	= $parent_node->findvalue("PDBx:pdb_seq_num");
			my $asym_id		= $parent_node->findvalue('@asym_id');
			my $hetero		= $parent_node->findvalue('PDBx:hetero');
			my $pdb_ins_code;
			if ($parent_node->findvalue('PDBx:pdb_ins_code/@xsi:nil') ne 'true') {
				$pdb_ins_code = $parent_node->findvalue('PDBx:pdb_ins_code');
			}
			if ($DEBUG > 1) { 
				print "DEBUG: $seq_id $pdb_seq_num $auth_seq_num\n";
			}
			if (defined($pdb_seq_num)) { 
				if ($hetero eq 'y' && $pdb_seq_nums{$chain}{$seq_id}) { next } 
				if ($pdb_ins_code) { 
					$pdb_seq_nums{$chain}{$seq_id} 	= $pdb_seq_num . $pdb_ins_code;
				}else{
					$pdb_seq_nums{$chain}{$seq_id}	= $pdb_seq_num;
				}
			}

			$asym_seq_nums{$chain_id}{$seq_id} = $asym_id;

			if ($auth_seq_num =~ /\d+/) { 
				push (@struct_seq_ids, $seq_id);
				push (@struct_chain_ids, $chain_id);
			}
			push (@seq_ids, $seq_id);

			$target_asym_id{$asym_id}++;

			push (@chain_ids, $chain_id);
		}
	}
	$target_asym_id = join(',', keys %target_asym_id);

	return (\@seq_ids, \@struct_seq_ids, \%pdb_seq_nums, $target_asym_id, \@chain_ids, \@struct_chain_ids, \%asym_seq_nums);
}
sub pdbml_annot_parse { 
	my $sub = 'pdbml_annot_parse';
	my ( $pdb4, $chain ) = @_;

	$pdb4 = lc($pdb4);
	my $two;
	if ($pdb4 =~ /\w(\w{2})\w/) { 
		$two = $1;
	}else{
		print "ERROR $sub: $pdb4 does not regexp\n";
		return 0;
	}

	if ($DEBUG) { 
		print "DEBUG $sub: $pdb4 $chain\n";
	}

	my $pdbml = pdbml_load($pdb4);
	(!$pdbml && return 0);
	#struct
	my %struct;
	my $structCategoryXPath = '//PDBx:structCategory/PDBx:struct';
	my $struct_nodes = $pdbml->findnodes($structCategoryXPath);

	foreach my $node ($struct_nodes->get_nodelist()) { 
		my $entry_id = $node->findvalue('@entry_id');
		my $pdbx_CASP_flag;
		if ($node->findvalue('PDBx:pdbx_CASP_flag/@xsi:nil') eq 'true') {
			$pdbx_CASP_flag = 'FALSE';
		}else{
			$pdbx_CASP_flag = $node->findvalue('PDBx:pdbx_CASP_flag');
		}

		my $pdbx_descriptor = $node->findvalue('PDBx:pdbx_descriptor');
		my $title = $node->findvalue('PDBx:title');
		#print "$entry_id\t$pdbx_descriptor\t$title\n";
		$struct{$entry_id}{pdbx_descriptor} = $pdbx_descriptor;
		$struct{$entry_id}{title}	= $title;
	}

	#struct_ref;
	my %struct_ref;
	my $struct_refCategoryXPath ='//PDBx:struct_refCategory/PDBx:struct_ref';
	my $struct_ref_nodes = $pdbml->findnodes($struct_refCategoryXPath);

	foreach my $node ($struct_ref_nodes->get_nodelist() ) { 
		
		my $struct_ref_id = $node->findvalue('@id');
		my $struct_ref_entity_id = $node->findvalue('PDBx:entity_id');
		my $struct_ref_db_code	= $node->findvalue('PDBx:db_code');
		my $struct_ref_db_name	= $node->findvalue('PDBx:db_name');
		my $struct_ref_db_accession = $node->findvalue('PDBx:pdbx_db_accession');

		$struct_ref{$struct_ref_id}{entity_id}	= $struct_ref_entity_id;
		$struct_ref{$struct_ref_id}{db_code}	= $struct_ref_db_code;
		$struct_ref{$struct_ref_id}{db_name}	= $struct_ref_db_name;
		$struct_ref{$struct_ref_id}{db_accession}	= $struct_ref_db_accession;

	}

	#struct_ref_seq
	my %struct_ref_seq;
	my $struct_ref_seqCategoryXPath = '//PDBx:struct_ref_seqCategory/PDBx:struct_ref_seq';
	my $struct_ref_seq_nodes = $pdbml->findnodes($struct_ref_seqCategoryXPath);

	foreach my $node ($struct_ref_seq_nodes->get_nodelist() ) { 

		my $align_id 			= $node->findvalue('@align_id');

		my $db_align_beg		= $node->findvalue('PDBx:db_align_beg');
		my $db_align_end		= $node->findvalue('PDBx:db_align_end');

		my $pdbx_PDB_id_code		= $node->findvalue('PDBx:pdbx_PDB_id_code');

		my $pdbx_auth_seq_align_beg	= $node->findvalue('PDBx:pdbx_auth_seq_align_beg');
		my $pdbx_auth_seq_align_end	= $node->findvalue('PDBx:pdbx_auth_seq_align_end');

		my $pdbx_db_accession		= $node->findvalue('PDBx:pdbx_db_accession');
		my $pdbx_strand_id		= $node->findvalue('PDBx:pdbx_strand_id');

		my $pdbx_ref_id			= $node->findvalue('PDBx:ref_id');

		my $seq_align_beg		= $node->findvalue('PDBx:seq_align_beg');
		my $seq_align_end		= $node->findvalue('PDBx:seq_align_end');

		$struct_ref_seq{$align_id}{db_align_beg}	= $db_align_beg;
		$struct_ref_seq{$align_id}{db_align_end}	= $db_align_end;

		$struct_ref_seq{$align_id}{pdbx_PDB_id_code}	= $pdbx_PDB_id_code;

		$struct_ref_seq{$align_id}{pdbx_auth_seq_align_beg} = $pdbx_auth_seq_align_beg;
		$struct_ref_seq{$align_id}{pdbx_auth_seq_align_end} = $pdbx_auth_seq_align_end;

		$struct_ref_seq{$align_id}{pdbx_db_accession}	= $pdbx_db_accession;
		$struct_ref_seq{$align_id}{pdbx_strand_id}	= $pdbx_strand_id;

		$struct_ref_seq{$align_id}{pdbx_ref_id}		= $pdbx_ref_id;

		$struct_ref_seq{$align_id}{seq_align_beg}	= $seq_align_beg;
		$struct_ref_seq{$align_id}{seq_align_end}	= $seq_align_end;
	}
		

	#exptl;
	my %exptl;
	my $exptlCategoryXPath = '//PDBx:exptlCategory/PDBx:exptl';
	my $exptl_nodes = $pdbml->findnodes($exptlCategoryXPath);
	my ($entry_id, $method);
	foreach my $node ($exptl_nodes->get_nodelist()) { 
		$entry_id = $node->findvalue('@entry_id');
		$method = $node->findvalue('@method');
		push(@{$exptl{$entry_id}{method}}, $method);
	}

	#citation
	my %citation;
	my $citationCategoryXPath = '//PDBx:citationCategory/PDBx:citation';
	my $citation_nodes = $pdbml->findnodes($citationCategoryXPath);
	foreach my $node ($citation_nodes->get_nodelist() ) { 
		my $citation_id = $node->findvalue('@id');
		my $journal_abbrev = $node->findvalue('PDBx:journal_abbrev');
		my $journal_id_CSD = 'NA';
		if ($node->exists('PDBx:journal_id_CSD')) { 
			$journal_id_CSD	= $node->findvalue('PDBx:journal_id_CSD');
		}

		my $journal_id_ISSN = 'NA';
		if ($node->exists('PDBx:journal_id_ISSN')) { 
			$journal_id_ISSN = $node->findvalue('PDBx:journal_id_ISSN');
		}
			
		my $title 	= $node->findvalue('PDBx:title');

		my $DOI		= 'NA';
		if ($node->exists('PDBx:pdbx_database_id_DOI')) { 
			$DOI		= $node->findvalue('PDBx:pdbx_database_id_DOI');
		}
		my $PMID	= 'NA';
		if ($node->exists('PDBx:pdbx_database_id_PubMed')) { 
			$PMID	= $node->findvalue('PDBx:pdbx_database_id_PubMed');
		}

		$citation{$citation_id}{journal_abbrev}		= $journal_abbrev;
		$citation{$citation_id}{journal_id_CSD}		= $journal_id_CSD;
		$citation{$citation_id}{journal_id_ISSN}	= $journal_id_ISSN;
		$citation{$citation_id}{title}			= $title;
		$citation{$citation_id}{DOE}			= $DOI;
		$citation{$citation_id}{PMID}			= $PMID;
	}

	
	
	#citation_author
	my %citation_author;
	my $citation_auhorCategoryXPath = '//PDBx:citation_authorCategory/PDBx:citation_author';
	my $citation_author_nodes = $pdbml->findnodes($citation_auhorCategoryXPath);
	my ($citation_id, $name, $ordinal);
	foreach my $node ($citation_author_nodes->get_nodelist()) { 
		my $citation_id = $node->findvalue('@citation_id');
		my $name = $node->findvalue('@name');
		my $ordinal = $node->findvalue('@ordinal');
		$citation_author{$citation_id}{$ordinal} = $name;
	}



	#entity
	my %entity;
	my $target_entity_id;
	my $entityCategoryXPath = '//PDBx:entityCategory/PDBx:entity';
	my $entity_nodes = $pdbml->findnodes($entityCategoryXPath);
	my @polymer_entities;
	foreach my $node ($entity_nodes->get_nodelist()) { 
		my $entity_id = $node->findvalue('@id');
		my $pdbx_description = $node->findvalue("PDBx:pdbx_description");
		my $pdbx_mutation = $node->findvalue("PDBx:pdbx_mutation");
		my $pdbx_fragment = $node->findvalue("PDBx:pdbx_fragment");
		my $pdbx_number_of_molecules = $node->findvalue("PDBx:pdbx_number_of_molecules");
		my $src_method 	= $node->findvalue("PDBx:src_method");
		my $type 		= $node->findvalue("PDBx:type");

		$entity{$entity_id}{pdbx_description} = $pdbx_description;
		$entity{$entity_id}{pdbx_mutation}	= $pdbx_mutation;
		$entity{$entity_id}{pdbx_fragment}	= $pdbx_fragment;
		$entity{$entity_id}{src_method}	= $src_method;
		$entity{$entity_id}{type}	= $type;
	}

	#entity_src_gen
	my %entity_src_gen;
	my $entity_src_genCategoryXPath = '//PDBx:entity_src_genCategory/PDBx:entity_src_gen';
	my $entity_src_gen_nodes = $pdbml->findnodes($entity_src_genCategoryXPath);
	foreach my $node ($entity_src_gen_nodes->get_nodelist() ) { 

		my $entity_id	= $node->findvalue('@entity_id');

		my $pdbx_gene_src_gene			= $node->findvalue('PDBx:pdbx_gene_src_gene');
		my $pdbx_gene_src_ncbi_taxonomy_id 	= $node->findvalue('PDBx:pdbx_gene_src_ncbi_taxonomy_id');
		my $pdbx_gene_src_scientific_name	= $node->findvalue('PDBx:pdbx_gene_src_scientific_name');
		my $pdbx_host_org_ncbi_taxonomy_id	= $node->findvalue('PDBx:pdbx_host_org_ncbi_taxonomy_id');
		my $pdbx_host_org_scientific_name	= $node->findvalue('PDBx:pdbx_host_org_scientific_name');
		my $pdbx_host_org_strain		= $node->findvalue('PDBx:pdbx_host_org_strain');
		my $pdbx_host_org_vector_type		= $node->findvalue('PDBx:pdbx_host_org_vector_type');

		$entity_src_gen{$entity_id}{pdbx_gene_src_gene}			= $pdbx_gene_src_gene;
		$entity_src_gen{$entity_id}{pdbx_gene_src_ncbi_taxonomy_id} 	= $pdbx_gene_src_ncbi_taxonomy_id;
		$entity_src_gen{$entity_id}{pdbx_gene_src_scientific_name}	= $pdbx_gene_src_scientific_name;
		$entity_src_gen{$entity_id}{pdbx_host_org_ncbi_taxonomy_id}	= $pdbx_host_org_ncbi_taxonomy_id;
		$entity_src_gen{$entity_id}{pdbx_host_org_scientific_name}	= $pdbx_host_org_scientific_name;
		$entity_src_gen{$entity_id}{pdbx_host_org_strain}		= $pdbx_host_org_strain;
		$entity_src_gen{$entity_id}{pdbx_host_org_vector_type}		= $pdbx_host_org_vector_type;
	}

	#entity_src_nat
	my %entity_src_nat;
	my $entity_src_natCategoryXPath = '//PDBx:entity_src_natCategory/PDBx:entity_src_nat';
	my $entity_src_nat_nodes = $pdbml->findnodes($entity_src_natCategoryXPath);
	foreach my $node ($entity_src_nat_nodes->get_nodelist() ) { 

		my $entity_id	= $node->findvalue('@entity_id');
		my $pdbx_ncbi_taxonomy_id	= $node->findvalue('PDBx:pdbx_ncbi_taxonomy_id');
		my $pdbx_organism_scientific	= $node->findvalue('PDBx:pdbx_organism_scientific');

		$entity_src_nat{$entity_id}{pdbx_ncbi_taxonomy_id} = $pdbx_ncbi_taxonomy_id;
		$entity_src_nat{$entity_id}{pdbx_organism_scientific} = $pdbx_organism_scientific;
	}

	#entity_src_syn;
	my %entity_src_syn;
	my $entity_src_synCategoryXPath = '//PDBx:entity_src_synCategory/PDBx:entity_src_syn';
	my $entity_src_syn_nodes = $pdbml->findnodes($entity_src_synCategoryXPath);
	foreach my $node ($entity_src_syn_nodes->get_nodelist() ) { 

		my $entity_id = $node->findvalue('@entity_id');
		my $details	= $node->findvalue('PDBx:details');

		$entity_src_syn{$entity_id}{details} = $details;
	}

	
	my $entity_polyCategoryXPath = '//PDBx:entity_polyCategory/PDBx:entity_poly';
	my $entity_poly_nodes = $pdbml->findnodes($entity_polyCategoryXPath);

	foreach my $node ($entity_poly_nodes->get_nodelist()) { 
		my $entity_id = $node->findvalue('@entity_id');
		my $nstd_linkage = $node->findvalue('PDBx:nstd_linkage');
		my $nstd_monomer = $node->findvalue('PDBx:nstd_monomer');
		my $pdbx_seq_can = $node->findvalue('PDBx:pdbx_seq_one_letter_code_can');
		$pdbx_seq_can =~ s/\n//g;
		my $type = $node->findvalue('PDBx:type');
		my $pdbx_strand_id = $node->findvalue('PDBx:pdbx_strand_id');
		if ($nstd_linkage eq 'yes') { 
			#print "WARNING $sub: nstd_linkage in $entry_id, $entity_id\n";
		}
		if ($nstd_monomer eq 'yes') { 
			#print "WARNING $sub: nstd_monomer in $entry_id, $entity_id\n";
		}
		if ($type ne 'polypeptide(L)') { 
			#print "WARNING $sub: $entry_id, $entity_id is not polypeptide(L)\n";
		}else{
			#if (length($pdbx_seq_can) > 30) { 
			#	print OUT ">$entry_id,$entity_id,$pdbx_strand_id\n$pdbx_seq_can\n";
			#}
		}
		my $length = length($pdbx_seq_can);
		#print "\t$entity_id\t$type\t$pdbx_strand_id\t$length\t$chain\n";
		$entity{$entity_id}{type} = $type;
		$entity{$entity_id}{pdbx_strand_id} = $pdbx_strand_id;
		$pdbx_strand_id =~ s/\s+//;
		
		if ($pdbx_strand_id =~ /\,/) { 
			my @strands = split(/\,/, $pdbx_strand_id);
			foreach my $strand (@strands) { 
				if ($strand eq $chain) { 
					$target_entity_id = $entity_id;
				}
			}
		}else{
			if ($pdbx_strand_id eq $chain) { 
				$target_entity_id = $entity_id;
			}
		}
	}

	return (\%struct, \%exptl, \%entity, $target_entity_id, \%struct_ref, \%struct_ref_seq, \%citation, \%citation_author, \%entity_src_nat, \%entity_src_gen, \%entity_src_syn);

}


sub pdbml_asym_annotate_parse { 

	my ($pdb4, $target_asym_id) = @_;

	$pdb4 =~ /\w(\w{2})\w/;
	my $two = $1;

	my $sub = 'pdbml_asym_annotate_parse';
	my $pdbml = pdbml_load($pdb4);

	my $struct_confCategoryXPath = '//PDBx:struct_confCategory/PDBx:struct_conf';

	my $struct_conf_nodes = $pdbml->findnodes($struct_confCategoryXPath);
	
	if ($DEBUG) { 
		print "DEBUG $sub: struct_conf nodes: " . scalar($struct_conf_nodes->get_nodelist()) . "\n";
	}


	my @helices;
	my $i = 0;
	foreach my $node ($struct_conf_nodes->get_nodelist()) { 

		my $id		= $node->findvalue('@id');

		my $beg_label_asym_id	= $node->findvalue('PDBx:beg_label_asym_id');
		my $end_label_asym_id	= $node->findvalue('PDBx:end_label_asym_id');

		if ($beg_label_asym_id ne $target_asym_id || $end_label_asym_id ne $target_asym_id) { next } 

		my $beg_label_seq_id	= $node->findvalue('PDBx:beg_label_seq_id');
		my $end_label_seq_id	= $node->findvalue('PDBx:end_label_seq_id');

		my $conf_type		= $node->findvalue('PDBx:conf_type_id');

		if ($DEBUG) { 
			print "DEBUG: $sub $id $beg_label_seq_id $end_label_seq_id $conf_type\n";
		}

		$helices[$i]{id}		= $id;
		$helices[$i]{beg_seq_id} 	= $beg_label_seq_id;
		$helices[$i]{end_seq_id}	= $end_label_seq_id;
		$helices[$i]{conf_type}		= $conf_type;

		$i++;
	}
		
	my $struct_sheet_rangeCategoryXPath = '//PDBx:struct_sheet_rangeCategory/PDBx:struct_sheet_range';

	my $struct_sheet_range_nodes = $pdbml->findnodes($struct_sheet_rangeCategoryXPath);

	if ($DEBUG) { 
		print "DEBUG $sub: struct_sheet_range nodes: " . scalar($struct_sheet_range_nodes->get_nodelist()) . "\n";
	}

	my @sheets;
	$i = 0;
	foreach my $node ($struct_sheet_range_nodes->get_nodelist()) { 
		
		my $id		= $node->findvalue('@id');
		my $sheet_id	= $node->findvalue('@sheet_id');

		my $beg_label_asym_id	= $node->findvalue('PDBx:beg_label_asym_id');
		my $end_label_asym_id	= $node->findvalue('PDBx:end_label_asym_id');

		if ($beg_label_asym_id ne $target_asym_id || $end_label_asym_id ne $target_asym_id) { next } 

		my $beg_label_seq_id	= $node->findvalue('PDBx:beg_label_seq_id');
		my $end_label_seq_id	= $node->findvalue('PDBx:end_label_seq_id');

		if ($DEBUG) { 
			print "DEBUG $sub: $id $sheet_id $beg_label_seq_id $end_label_seq_id\n";
		}

		$sheets[$i]{id}		= $id;
		$sheets[$i]{sheet_id}	= $sheet_id;

		$sheets[$i]{beg_seq_id}	= $beg_label_seq_id;
		$sheets[$i]{end_seq_id}	= $end_label_seq_id;

		$i++;
	}

	my $struct_siteCategoryXPath = '//PDBx:struct_siteCategory/PDBx:struct_site';

	my $struct_site_nodes	= $pdbml->findnodes($struct_siteCategoryXPath);

	my %struct_site;
	foreach my $node ($struct_site_nodes->get_nodelist()) { 

		my $id		= $node->findvalue('@id');
		
		my $details	= $node->findvalue('PDBx:details');
		my $pdbx_evidence_code	= $node->findvalue('PDBx:pdbx_evidence_code');

		$struct_site{$id}{details} = $details;
		$struct_site{$id}{evidence}	= $pdbx_evidence_code;

	}


	

	my $struct_site_genCategoryXPath = '//PDBx:struct_site_genCategory/PDBx:struct_site_gen';

	my $struct_site_gen_nodes = $pdbml->findnodes($struct_site_genCategoryXPath);

	my @struct_site_gen;
	$i = 0;
	foreach my $node ($struct_site_gen_nodes->get_nodelist()) { 

		my $id		= $node->findvalue('@id');
		my $site_id	= $node->findvalue('@site_id');

		my $label_asym_id	= $node->findvalue('PDBx:label_asym_id');

		if ($label_asym_id ne $target_asym_id) { next } 

		my $label_seq_id	= $node->findvalue('PDBx:label_seq_id');

		my $label_comp_id	= $node->findvalue('PDBx:label_comp_id');
		my $symmetry		= $node->findvalue('PDBx:symmetry');

		if ($DEBUG) { 
			print "DEBUG $sub: $id $site_id $label_seq_id $label_comp_id $symmetry $struct_site{$site_id}{details} $struct_site{$site_id}{evidence}\n";
		}

		$struct_site_gen[$i]{id}		= $id;
		$struct_site_gen[$i]{site_id}		= $site_id;
		$struct_site_gen[$i]{seq_id}		= $label_seq_id;
		$struct_site_gen[$i]{comp_id}		= $label_comp_id;
		$struct_site_gen[$i]{symmetry}		= $symmetry;
		$struct_site_gen[$i]{details}		= $struct_site{$site_id}{details};
		$struct_site_gen[$i]{evidence}		= $struct_site{$site_id}{evidence};

		$i++;
	}

	my $struct_connCategoryXPath = '//PDBx:struct_connCategory/PDBx:struct_conn';
	
	my $struct_conn_nodes	= $pdbml->findnodes($struct_connCategoryXPath);

	$i = 0;
	my @disulf;

	foreach my $node ($struct_conn_nodes->get_nodelist()) { 

		my $id			= $node->findvalue('@id');

		my $ptnr1_label_asym_id	= $node->findvalue('PDBx:ptnr1_label_asym_id');
		my $ptnr2_label_asym_id	= $node->findvalue('PDBx:ptnr2_label_asym_id');

		if ($ptnr1_label_asym_id ne $target_asym_id && $ptnr2_label_asym_id ne $target_asym_id) { next } 

		my $conn_type_id	= $node->findvalue('PDBx:conn_type_id');

		my $ptnr1_label_seq_id	= $node->findvalue('PDBx:ptnr1_label_seq_id');
		my $ptnr2_label_seq_id	= $node->findvalue('PDBx:ptnr2_label_seq_id');

		my $ptnr1_label_comp_id	= $node->findvalue('PDBx:ptnr1_label_comp_id');
		my $ptnr2_label_comp_id	= $node->findvalue('PDBx:ptnr2_label_comp_id');

		if ($ptnr1_label_comp_id eq 'MSE' || $ptnr2_label_comp_id eq 'MSE') { next } 

		if ($DEBUG) { 
			print "DEBUG $sub: $id $conn_type_id $ptnr1_label_seq_id $ptnr1_label_comp_id $ptnr2_label_seq_id $ptnr2_label_comp_id\n";
		}

		$disulf[$i]{id}		= $id;
		$disulf[$i]{conn_type}	= $conn_type_id;

		$disulf[$i]{ptnr1_seq_id}	= $ptnr1_label_seq_id;
		$disulf[$i]{ptnr2_seq_id}	= $ptnr2_label_seq_id;

		$disulf[$i]{ptnr1_asym_id}	= $ptnr1_label_asym_id;
		$disulf[$i]{ptnr2_asym_id}	= $ptnr2_label_asym_id;

		$i++;
	}
	return(\@helices, \@sheets, \@disulf, \@struct_site_gen);
}

sub pdbml_date_method { 

	my $sub = 'pdbml_date_method';
	my ($pdb_id) = @_;
	$pdb_id = lc($pdb_id);
	my $pdbml = pdbml_load($pdb_id);

	#exptl
	my %exptl;
	my $exptlCategoryXPath = '//PDBx:exptlCategory/PDBx:exptl';
	my $exptl_nodes = $pdbml->findnodes($exptlCategoryXPath);
	my ($entry_id, $method);
	foreach my $node ($exptl_nodes->get_nodelist()) { 
		$entry_id = $node->findvalue('@entry_id');
		$method = $node->findvalue('@method');
		push (@{$exptl{$entry_id}{method}}, $method);
	}

	#deposition date (mod = 0)  
	my $date;
	my %database_PDB_rev;
	my $database_PDB_revCategory_XPath = '//PDBx:database_PDB_revCategory/PDBx:database_PDB_rev[@num="1"]';
	if ($pdbml->exists($database_PDB_revCategory_XPath)) { 

		my $rev_node = $pdbml->findnodes($database_PDB_revCategory_XPath)->get_node(1);
	#	if ($rev_node->exists('PDBx:date_original')) { 
	#		$date = $rev_node->findvalue('PDBx:date_original');
	#	}elsif($rev_node->exists('PDBx:date')) { 
			$date = $rev_node->findvalue('PDBx:date');
	#	}else{
	#		die "$pdb_id?";
	#	}
	}
	my $uc_pdb_id = uc($pdb_id);
	#print "$pdb_id $exptl{$uc_pdb_id}{method} $date\n";

	return ($date, $exptl{$uc_pdb_id}{method});

}
1;


