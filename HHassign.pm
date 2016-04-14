package Domains::HHassign;
require Exporter;

use warnings;
use strict;

use Domains::Partition;
use Domains::PDBml;
use Domains::Range;
use XML::Grishin;
#use ECOD::Reference

umask(0002);
my $UID = 1219;
my $GID = 902;

our @ISA    = qw(Exporter);
our @EXPORT = ( "&domain_assignment_by_hhsearch" );

my $DEBUG = 0;

my $HH_CUTOFF = 70; #Probability threshold for no-partition assignment

sub domain_assignment_by_hhsearch {
    my $sub = 'domain_assignment_by_hhsearch';

    my (
        $query_pdb,  $query_chain,       $pdb_chain_dir,    $input_mode,
        $reference,  $domain_xml_fn,     $ref_range_cache, $concise,
        $global_opt, $ref_chain_domains, $ref_domain_uid_lookup
    ) = @_;

    if ( -f $domain_xml_fn && !$$global_opt{force_overwrite} ) {
        print "WARNING! $sub: skipping $query_pdb, $query_chain, set FORCE_OVERWRITE bit to override...\n";
        return;
    }

    #Return the seqids and struct seqids of the query pdb_chain
    my ( $query_seqid_aref, $query_struct_seqid_aref, $query_pdbnum_aref, $query_asym_id ) =
      pdbml_seq_parse( $query_pdb, $query_chain );

    if ( !$query_seqid_aref ) {
        print "WARNING! $sub: $query_pdb not found! skipping...\n";
        return 0;
    }
    my $query_range_str        = rangify(@$query_seqid_aref);
    my $query_struct_range_str = rangify(@$query_struct_seqid_aref);

#Changed query sequence to be ungapped (i.e. may contain small unstructured regions, extremely important that GAP_TOL is the same for query generation as domain partition).
    my $query_ungapped_struct_range_str = ungap_range( $query_struct_range_str, $$global_opt{gap_tol} );
    my @gapped_query_struct_seqid =
      @$query_struct_seqid_aref;    #Need to retain this for dali_fillin (query pdbs are not ungapped)
    $query_struct_seqid_aref = range_expand($query_ungapped_struct_range_str);

    if ( !$query_struct_range_str ) {
        print "WARNING! $sub: $query_pdb $query_chain,  No structured range found =>  $query_struct_range_str\n";
    }

    #Return annotations for the query pdb_chain
    my (
        $struct_href,         $exptl_href,    $entity_href,          $target_entity_id,
        $struct_ref_href,     $citation_href, $citation_author_href, $entity_src_nat_href,
        $entity_src_gen_href, $entity_src_syn_href
    ) = pdbml_annot_parse( $query_pdb, $query_chain );

    my $pdb_chain = $query_pdb . "_" . $query_chain;

    my $src_method = $$entity_href{$target_entity_id}{src_method};

    my @unused_seq =
        $input_mode eq 'struct_seqid'
      ? @$query_struct_seqid_aref
      : @$query_seqid_aref;
    my @used_seq;

    if ($DEBUG) {
        print "DEBUG $sub: $query_pdb $query_chain $pdb_chain_dir $query_range_str $query_struct_range_str $domain_xml_fn\n";
    }

    my $domain_count = 1;    #Used for generating new domain ids

    my $domain_xml_doc  = XML::LibXML->createDocument;
    my $domain_doc_node = $domain_xml_doc->createElement('chain_domains_set_top');
    $domain_xml_doc->setDocumentElement($domain_doc_node);

    $domain_doc_node->setAttribute( 'pdb_id',        $query_pdb );
    $domain_doc_node->setAttribute( 'chain_id',      $query_chain );
    $domain_doc_node->setAttribute( 'mode',          'hh_assign' );
    $domain_doc_node->setAttribute( 'domain_prefix', $DOMAIN_PREFIX );

    my $lc_query_pdb = lc($query_pdb);

    if ( $src_method eq 'syn' ) {
        print "WARNING! $query_pdb, $query_chain is a known synthetic product\n";
        $domain_doc_node->setAttribute( 'known_synthetic', 'true' );
    }

    #Do some work for coil detection
    my $coils_fn = "$pdb_chain_dir/$pdb_chain.coils.xml";
    my @coils;
    if ( -f $coils_fn ) {
        print "Found coils definition for $pdb_chain\n";
        my $coils_xml_doc = xml_open($coils_fn);
        foreach my $coil_node ( $coils_xml_doc->findnodes('//coil') ) {
            my $coil_seqid_range = $coil_node->findvalue('coil_seqid_range');
            next unless $coil_seqid_range =~ /\d+/;
            my $method = $coil_node->findvalue('@method');
            my %coil_entry;
            $coil_entry{method} = $method;
            my ( $range_str, $chain_str ) = scop_range_split($coil_seqid_range);
            $coil_entry{seqid_range} = $range_str;
            push( @coils, \%coil_entry );
        }
    }

    #Build the stub of the domains xml file
    my $domain_list_node = $domain_xml_doc->createElement('domain_list');
    $domain_doc_node->appendChild($domain_list_node);

    #Open the sequence summary search file (i.e. the "blast_summ" file)
    #my $blast_summ_xml_doc = xml_open($blast_summ_fn);
    my $hhsearch_fn = "$pdb_chain_dir/$pdb_chain.$reference.rebuild.hh_summ.xml";
    if ( !-f $hhsearch_fn ) {
        die "ERROR! $sub: Could not find hhsearch file: $hhsearch_fn\n";
    }
    my $hhsearch_xml_doc = xml_open($hhsearch_fn);

    assign_hhsearch_domains(
        $hhsearch_xml_doc, $domain_xml_doc,          $ref_chain_domains, $ref_domain_uid_lookup,
        $ref_range_cache,  $query_struct_seqid_aref, $query_seqid_aref,  $query_pdbnum_aref,
         \@used_seq,               \@unused_seq,       $input_mode,
        $global_opt,       $reference,               
    );

    my $final_coverage = region_coverage( $query_struct_seqid_aref, \@used_seq );

    my $chain_domain_coverage_node = $domain_xml_doc->createElement('chain_domain_coverage');
    $chain_domain_coverage_node->appendTextNode($final_coverage);
    $chain_domain_coverage_node->setAttribute( 'used_res',   scalar @used_seq );
    $chain_domain_coverage_node->setAttribute( 'unused_res', scalar @unused_seq );
    $domain_doc_node->appendChild($chain_domain_coverage_node);

    xml_write( $domain_xml_doc, $domain_xml_fn );

    return 1;

}

#Assign entire unassigned range to top hhsearch hit regardless of coverage. Perilous. Use
#only as advised.
sub assign_hhsearch_domains {
    my $sub = 'assign_hhsearch_domains';

    my (
        $hhsearch_xml_doc, $domain_xml_doc,          $ref_chain_domains, $ref_domain_uid_lookup,
        $ref_range_cache,  $query_struct_seqid_aref, $query_seqid_aref,  $query_pdbnum_aref,
           $used_seq_aref,           $unused_seq_aref,   $input_mode,
        $global_opt,       $reference    ) = @_;

    my $domain_count =
        $domain_xml_doc->exists('//domain')
      ? $domain_xml_doc->findnodes('//domain')->size() + 1
      : 1;

    my $domain_list_node = $domain_xml_doc->findnodes('//domain_list')->get_node(1);

	my ($query_pdb, $query_chain) = get_pdb_chain($domain_xml_doc);

    foreach my $hit_node ( find_hh_hit_nodes($hhsearch_xml_doc)) {
        next if isObsolete($hit_node);

		my ($hit_uid, $hit_domain_id) = get_ids($hit_node);

        unless ( passesCoverageCheck($hit_node), $$global_opt{hit_cover} ) {
            print "WARNING! $hit_uid, $hit_domain_id has low hit coverage\n";
        }

        my ( $hh_prob, $hh_score ) = get_hh_scores($hit_node);

		my ($hit_query_reg, $hit_hit_seqid) = get_hh_hit_regions($hit_node);

        my $hit_query_struct_seqid_aref;
        if ( $input_mode eq 'seqid' ) {
            $hit_query_struct_seqid_aref = struct_region( range_expand($hit_query_reg), $query_seqid_aref );
        }
        elsif ( $input_mode eq 'struct_seqid' ) {
            $hit_query_struct_seqid_aref = struct_region( range_expand($hit_query_reg), $query_struct_seqid_aref );
        }
        else {
            die "ERROR! Unknown input_mode $input_mode\n";
        }

        my $query_struct_seqid = rangify(@$query_struct_seqid_aref);
        if ( !$query_struct_seqid_aref || !$query_struct_seqid ) {
            printf
"WARNING! No structured query range? $hit_domain_id $hit_query_reg $query_struct_seqid $$query_struct_seqid_aref[0] $$query_struct_seqid_aref[1] %i\n",
              scalar(@$query_struct_seqid_aref);
            next;
        }

        if (   !$hit_query_struct_seqid_aref
            || !rangify(@$hit_query_struct_seqid_aref) )
        {
            print "WARNING! No range for hit $hit_uid $hit_domain_id $query_pdb $query_chain\n";
            next;
        }

        my $hit_query_struct_seqid = rangify(@$hit_query_struct_seqid_aref);
        if ($DEBUG) {
            print
"DEBUG $sub: hhsearch $hit_uid $hit_domain_id socore:$hh_score posx:$hit_query_reg seqid:$hit_query_struct_seqid\n";
        }

        my $unused_coverage = region_coverage( $query_seqid_aref, $unused_seq_aref );

        if ( $unused_coverage < 0.05 || scalar(@$unused_seq_aref) < 10 ) {

            #print "query complete\n";
            last;
        }
        my $query_coverage = region_coverage( $hit_query_struct_seqid_aref, $unused_seq_aref );
        my $query_used_coverage = residue_coverage( $hit_query_struct_seqid_aref, $used_seq_aref );

        if ( $hh_prob > $HH_CUTOFF 
            && scalar(@$hit_query_struct_seqid_aref) > $$global_opt{gap_tol} )
        {

            print "DEFINE hh_assign: $sub: $hit_domain_id $hit_query_struct_seqid\n";

            my $domain_node = $domain_xml_doc->createElement('domain');
            my $domain_id   = "e" . lc($query_pdb) . "$query_chain$domain_count";
            $domain_node->setAttribute( 'domain_id', $domain_id );

            my $method      = 'hh_assign';
            my $method_node = $domain_xml_doc->createElement('method');
            $method_node->appendTextNode($method);
            $domain_node->appendChild($method_node);

            my $struct_seqid_range_node = $domain_xml_doc->createElement('struct_seqid_range');
            $struct_seqid_range_node->appendTextNode($query_struct_seqid);
            $domain_node->appendChild($struct_seqid_range_node);

            my $hit_struct_seqid_range_node = $domain_xml_doc->createElement('hit_struct_seqid_range');
            $hit_struct_seqid_range_node->appendTextNode($hit_hit_seqid);
            $domain_node->appendChild($hit_struct_seqid_range_node);

            my $ungapped_seqid_range      = ungap_range( $query_struct_seqid, $$global_opt{gap_tol} );
            my $ungapped_seqid_range_aref = range_expand($ungapped_seqid_range);
            my $ungapped_seqid_range_node = $domain_xml_doc->createElement('ungapped_seqid_range');
            $ungapped_seqid_range_node->setAttribute( 'gap_tolerance', $$global_opt{gap_tol} );
            $ungapped_seqid_range_node->appendTextNode($ungapped_seqid_range);
            $domain_node->appendChild($ungapped_seqid_range_node);

            my $pdb_range = pdb_rangify( $query_struct_seqid_aref, $query_pdbnum_aref );
            my $pdb_range_node = $domain_xml_doc->createElement('pdb_range');
            $pdb_range_node->appendTextNode($pdb_range);
            $domain_node->appendChild($pdb_range_node);

            my $ungapped_pdb_range = pdb_rangify( $ungapped_seqid_range_aref, $query_pdbnum_aref );
            my $ungapped_pdb_range_node = $domain_xml_doc->createElement('ungapped_pdb_range');
            $ungapped_pdb_range_node->setAttribute( 'gap_tolerance', $$global_opt{gap_tol} );
            $ungapped_pdb_range_node->appendTextNode($ungapped_pdb_range);
            $domain_node->appendChild($ungapped_pdb_range_node);

            my $hit_domain_id_node = $domain_xml_doc->createElement('hit_domain');
            $hit_domain_id_node->setAttribute( 'ecod_domain_id', $hit_domain_id );
            $hit_domain_id_node->setAttribute( 'reference',      $reference );
            $hit_domain_id_node->setAttribute( 'uid',            $$ref_domain_uid_lookup{$hit_domain_id} );

            $domain_node->appendChild($hit_domain_id_node);

            my $hh_score_node = $domain_xml_doc->createElement('hh_score');
            $hh_score_node->setAttribute( 'hh_prob',  $hh_prob );
            $hh_score_node->setAttribute( 'hh_score', $hh_score );
            $domain_node->appendChild($hh_score_node);

            $domain_list_node->appendChild($domain_node);

            range_exclude( $ungapped_seqid_range_aref, $unused_seq_aref );
            range_include( $ungapped_seqid_range_aref, $used_seq_aref );

            $domain_count++;

        }
        else {
            if ($DEBUG) {
                print
                  "query_coverage $query_coverage query_used_coverage $query_used_coverage hh_prob $hh_prob scalar "
                  . scalar(@$hit_query_struct_seqid_aref) . "\n";
            }

        }
    }
}

sub isObsolete {
    $_[0]->findvalue('@structure_obsolete') eq 'true';
}

sub get_hh_scores {
    ( $_[0]->findvalue('@hh_prob'), $_[0]->findvalue('@hh_score') );
}

sub passesCoverageCheck {
    my $hit_cover = $_[0]->findvalue('@hit_cover');
    return 1 if $hit_cover > $_[1];
    return 0;
}

sub find_hh_hit_nodes { 
	$_[0]->findnodes('//hh_hit');
}

sub get_hh_hit_regions { 
	($_[0]->findvalue('query_range'), $_[0]->findvalue('template_seqid_range'));
}

sub get_pdb_chain { 
	($_[0]->findvalue('//@pdb_id'), $_[0]->findvalue('//@chain_id'));
}

1;
