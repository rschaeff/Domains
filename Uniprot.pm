package Domains::Uniprot;
require Exporter;

use warnings;
use strict;

use Carp;
use XML::Grishin;
use Domains::Range;
use LWP::UserAgent;
use JSON::XS;
use XML::LibXML;

my $FORCE_OVERWRITE = 0;

our @ISA = 'Exporter';

our @EXPORT = (
				"&chainwise_chain_to_domain",
				"&pdb_unp_to_pdb_chain_unp",
				"&pdb_chain_to_chainwise_chain_unp",
				"&chainwise_domain_to_ecod_domain",
				"&assemble_chainwise_domain_fragments",
				"&sifts_pdb_update",
				);

sub sifts_pdb_update { 
	my $SIFTS_API_URL = "http://www.ebi.ac.uk/pdbe/api/mappings/uniprot/"; 

	my ($pdb_xml_doc) = @_;

	my $ua = new LWP::UserAgent;
	$ua->agent("auto_uniprot_to_ecod/1.0". $ua->_agent);
	$ua->from('ecod.database@gmail.com');
	$ua->proxy(['http', 'ftp'], "http://proxy.swmed.edu:3128");

	$ENV{http_proxy} = "proxy.swmed.edu:3128";
	my $i = 0;
	foreach my $pdb_node (find_pdb_nodes($pdb_xml_doc)) { 
		next unless has_any_chain($pdb_node);
		next if is_obsolete($pdb_node);
		if (has_unp_mapping($pdb_node) && !$FORCE_OVERWRITE) { next } 
		my $pdb_id = get_pdb_id($pdb_node);
		my $pdb_url = $SIFTS_API_URL . $pdb_id;
		if ($pdb_node->exists('unp_mapping')) { 
			next;
		}

		#print "url:$pdb_url\n";
		my $req = HTTP::Request->new(GET => $pdb_url);
		$req->content_type('application/x-www-form-urlencoded');

		my $response 	= $ua->request($req);
		my $json	= decode_json($response->content); 
		
		my $href = $$json{$pdb_id}{UniProt};
		if (keys %$href > 0) { 
			my @unp_ids = keys %$href;
			foreach my $unp_id (keys %$href) { 
				my $unp_node = $pdb_xml_doc->createElement('unp_mapping');
				$unp_node->setAttribute('unp_acc', $unp_id);
				$unp_node->setAttribute('id', $$href{$unp_id}{identifier});
				$unp_node->setAttribute('name', $$href{$unp_id}{name});

				my $aref = $$href{$unp_id}{mappings};
				my $unp_range 			= get_unp_range_from_sifts_json($aref);
				my $residue_range 		= get_res_range_from_sifts_json($aref); 
				my $pdb_residue_range  	= get_pdb_res_range_from_sifts_json($aref);
				$unp_node->appendTextChild('unp_range', $unp_range);
				$unp_node->appendTextChild('residue_range', $residue_range);
				$unp_node->appendTextChild('pdb_residue_range', $pdb_residue_range);
				$pdb_node->appendChild($unp_node);
			}
		}else{
			#print "No mappings!\n";
		}
	}
}

sub has_unp_mapping { 
	$_[0]->exists('unp_mapping');
}

sub is_obsolete { 
	$_[0]->findvalue('@structure_obsolete') eq 'true';
}

sub domain_range_report { 
	my ($domain_node) = @_;

	my ($domain_xml_doc, $domain_root_node) = xml_create('domain_residue_doc');

	my ($uid, $ecod_domain_id) = get_ids($domain_node);

	my $residue_list = $domain_xml_doc->createElement('residue_list');
	$domain_root_node->appendChild($residue_list);
	$domain_root_node->setAttribute('uid', $uid);
	$domain_root_node->setAttribute('ecod_domain_id', $ecod_domain_id);

	my $seqid_range = get_seqid_range($domain_node);

	my ($seqid_aref, $chain_aref) = multi_chain_range_expand($seqid_range);

	my $unp_href = get_unp_domain_ranges($domain_node);
	return 0 unless $unp_href;
	
	for (my $i = 0; $i < scalar(@$seqid_aref); $i++) { 
		my $seqid = $$seqid_aref[$i];
		my $chain = $$chain_aref[$i];

		my $residue_node = $domain_xml_doc->createElement('residue');
		$residue_node->setAttribute('seq_id', $seqid);
		$residue_node->setAttribute('chain_id', $chain);	
		$residue_list->appendChild($residue_node);
		#print "$seqid $chain\n";

		foreach my $unp_acc (keys %$unp_href) { 
			if (exists $$unp_href{$unp_acc}{$chain} && defined $$unp_href{$unp_acc}{$chain}[$seqid]) { 
				#print "$unp_acc $seqid\n";
				my $unp_residue_node = $domain_xml_doc->createElement('unp_residue');
				$unp_residue_node->setAttribute('unp_acc', $unp_acc);
				$unp_residue_node->setAttribute('unp_id', $$unp_href{$unp_acc}{$chain}[$seqid]);
				$residue_node->appendChild($unp_residue_node);
			}else{
				#warn "WARNING! Position $i in $uid has no unp map\n";
			}
		}
	}
	return $domain_xml_doc;
}
sub get_unp_range_from_sifts_json { 
	my ($aref) = @_;
	my @segs;
	foreach my $href (@$aref) { 
		my $r = $$href{unp_start} . "-" . $$href{unp_end};
		push (@segs, $r);
	}
	my $r = join (",", @segs);
	return $r;
}

sub get_res_range_from_sifts_json { 
	my $aref = shift;
	my @segs;
	foreach my $href (@$aref) { 
		my $chain_id = $$href{chain_id};
		my $s = $chain_id . ":" . $$href{start}{residue_number};
		my $e = $$href{end}{residue_number};
		push (@segs, "$s-$e");
	}
	my $r = join(",", @segs);
	return $r;
}
sub get_pdb_res_range_from_sifts_json { 
	my $aref = shift;
	my @segs;
	foreach my $href (@$aref) { 
		my $chain_id = $$href{chain_id};
		my $s = $chain_id . ":" . $$href{start}{residue_number} . $$href{start}{author_insertion_code} ;
		my $e = $$href{end}{residue_number} . $$href{start}{author_insertion_code} ;
		push (@segs, "$s-$e");
	}
	my $r = join(",", @segs);
	return $r;

}

sub chainwise_chain_to_domain {
    my ($chainwise_xml_doc) = @_;
    foreach my $chain_node ( find_chain_nodes($chainwise_xml_doc) ) {

        my $chain_id = get_chain_id($chain_node);

        if ( has_uniprot_mapping($chain_node) ) {
            foreach my $unp_node ( find_unp_mapping_nodes($chain_node) ) {
                my $unp_range       = get_unp_range($unp_node);
                my $unp_seqid_range = get_seq_range($unp_node);
                my $pdb_range       = get_unp_pdb_range($unp_node);

                my $unp_id   = get_unp_id($unp_node);
                my $unp_name = get_unp_name($unp_node);
                my $unp_acc  = get_unp_acc($unp_node);

                #print "1:$unp_range $unp_seqid_range $pdb_range\n";

                my ( $unp_seqid_aref, $unp_chain_aref ) = multi_chain_range_expand($unp_seqid_range);
                my ($unp_range_aref) = range_expand($unp_range);

                my $display_debug = 0;
                if ( scalar @$unp_range_aref != scalar @$unp_seqid_aref ) {
                    my $size1 = scalar(@$unp_range_aref);
                    my $size2 = scalar(@$unp_seqid_aref);

                    #warn "WARNING! $unp_acc aref size diff $size1 $size2 \n";
                    #$display_debug++;
                    $unp_node->setAttribute( 'mapping_discrepancy', 'true' );
                    next;
                }

                foreach my $domain_node ( find_domain_nodes($chain_node) ) {
                    my $domain_seqid_range = get_seqid_range($domain_node);

                    my ( $domain_seqid_aref, $domain_chain_aref ) = multi_chain_range_expand($domain_seqid_range);
                    my $coverage = multi_chain_region_coverage( $unp_seqid_aref, $unp_chain_aref,
                        $domain_seqid_aref, $domain_chain_aref );

                    my $domain_unp_node = $chainwise_xml_doc->createElement('uniprot');
                    $domain_unp_node->setAttribute( 'unp_acc', $unp_acc );
                    $domain_unp_node->setAttribute( 'id',      $unp_id );
                    $domain_unp_node->setAttribute( 'name',    $unp_name );
                    my $found = 0;

                    if ( $coverage == 1 ) {
                        $domain_unp_node->appendTextChild( 'unp_range',   $unp_range );
                        $domain_unp_node->appendTextChild( 'seqid_range', $unp_seqid_range );
                        $found = 1;
                    }
                    elsif ( $coverage > 0 ) {
                        my (@domain_unp_range);
                        my ( @seqid_unp_range, @seqid_unp_chain );
                        for ( my $i = 0 ; $i < scalar(@$domain_seqid_aref) ; $i++ ) {
                            for ( my $j = 0 ; $j < scalar(@$unp_seqid_aref) ; $j++ ) {
                                if (   $$domain_seqid_aref[$i] == $$unp_seqid_aref[$j]
                                    && $$domain_chain_aref[$i] eq $$unp_chain_aref[$j] )
                                {
                                    if (   !defined $$domain_seqid_aref[$i]
                                        || !defined $$domain_chain_aref[$i] )
                                    {
                                        confess
"!!!!!$i $$domain_seqid_aref[$i] $$domain_chain_aref[$i] $domain_seqid_range\n";
                                    }
                                    push( @domain_unp_range, $$unp_range_aref[$j] );
                                    push( @seqid_unp_range,  $$domain_seqid_aref[$i] );
                                    push( @seqid_unp_chain,  $$domain_chain_aref[$i] );

                                }
                            }
                        }
                        my $domain_unp_range = rangify(@domain_unp_range);
                        $domain_unp_node->appendTextChild( 'unp_range', $domain_unp_range );

                        my $seqid_unp_range = multi_chain_rangify( \@seqid_unp_range, \@seqid_unp_chain );
                        $domain_unp_node->appendTextChild( 'seqid_range', $seqid_unp_range );
						$found = 1;
                    }
                    else {
					     #Do nothing
                    }
                    $domain_node->appendChild($domain_unp_node) if $found;
                }
            }
        }
    }
}

#PDB pdb to PDB chain
sub pdb_unp_to_pdb_chain_unp {
    my ($pdb_xml_doc) = @_;
    foreach my $pdb_node ( find_pdb_nodes($pdb_xml_doc) ) {
        my $pdb_id = get_pdb_id($pdb_node);
        my %chains;
        my %unp_acc_ids;
        my %unp_acc_names;
        if ( has_uniprot_mapping($pdb_node) ) {
            foreach my $unp_node ( find_unp_mapping_nodes($pdb_node) ) {
                my $unp_acc = get_unp_acc($unp_node);

                my $unp_id = get_unp_id($unp_node);
                $unp_acc_ids{$unp_acc} = $unp_id;
                my $unp_name = get_unp_name($unp_node);
                $unp_acc_names{$unp_acc} = $unp_name;

                my $unp_range = get_unp_range($unp_node);
                my $seq_range = get_res_range($unp_node);
                my $pdb_range = get_unp_pdb_range($unp_node);

                my @unp_segs = split ',', $unp_range;
                my @seq_segs = split ',', $seq_range;
                my @pdb_segs = split ',', $pdb_range;

                if (   scalar @unp_segs != scalar @seq_segs
                    || scalar @seq_segs != scalar @pdb_segs )
                {
                    printf "%i %i %i\n", scalar @unp_segs, scalar @seq_segs, scalar @pdb_segs;
                    die "ERROR! $pdb_id has UNP segment issues\n";
                }
                for ( my $i = 0 ; $i < scalar(@unp_segs) ; $i++ ) {
                    if ( $seq_segs[$i] =~ /((\w+):\d+\-\d+)/ ) {
                        my $seg_chain     = $2;
                        my $seg_seq_range = $1;
                        push( @{ $chains{$seg_chain}{$unp_acc}{seq} }, $seg_seq_range );
                        push( @{ $chains{$seg_chain}{$unp_acc}{unp} }, $unp_segs[$i] );
                        push( @{ $chains{$seg_chain}{$unp_acc}{pdb} }, $pdb_segs[$i] );
                    }
                }
            }
        }
        foreach my $unp_chain ( keys %chains ) {
            if ( has_chain( $pdb_node, $unp_chain ) ) {
                my $chain_node = get_chain_node( $pdb_node, $unp_chain );
                my $unp_href = $chains{$unp_chain};
                foreach my $unp_acc ( keys %$unp_href ) {
                    my $unp_chain_node = $pdb_xml_doc->createElement('unp_mapping');
                    $unp_chain_node->setAttribute( 'unp_acc', $unp_acc );
                    my $unp_id = $unp_acc_ids{$unp_acc};
                    $unp_chain_node->setAttribute( 'id', $unp_id );
                    my $unp_name = $unp_acc_names{$unp_acc};
                    $unp_chain_node->setAttribute( 'name', $unp_acc );
                    my $seq_range = join ',', @{ $$unp_href{$unp_acc}{seq} };
                    my $unp_range = join ',', @{ $$unp_href{$unp_acc}{unp} };
                    my $pdb_range = join ',', @{ $$unp_href{$unp_acc}{pdb} };
                    $unp_chain_node->appendTextChild( 'seq_range', $seq_range );
                    $unp_chain_node->appendTextChild( 'unp_range', $unp_range );
                    $unp_chain_node->appendTextChild( 'pdb_range', $pdb_range );

                    $chain_node->appendChild($unp_chain_node);
                }
            }
        }
    }
}

#PDB to chainwise
sub pdb_chain_to_chainwise_chain_unp {
    my ( $pdb_xml_doc, $chainwise_xml_doc ) = @_;
    my %chains;
    foreach my $chain_node ( find_chain_nodes($pdb_xml_doc) ) {

        my $chain_id = get_chain_id($chain_node);
        my $pdb_id   = get_pdb_id( $chain_node->parentNode );

        if ( $chain_node->exists('unp_mapping') ) {
            $chains{$pdb_id}{$chain_id} = $chain_node;
        }
    }

    foreach my $chain_node ( find_chain_nodes($chainwise_xml_doc) ) {
        my $chain_id = get_chain_id($chain_node);
        my $pdb_id   = get_pdb_id( $chain_node->parentNode );
        if ( $chains{$pdb_id}{$chain_id} ) {
            foreach my $unp_node ( $chains{$pdb_id}{$chain_id}->findnodes('unp_mapping') ) {
                $chain_node->appendChild($unp_node);
            }
        }
    }
}

sub chainwise_domain_to_ecod_domain {

    #chainwise chain to ecod domain
    my ( $chainwise_xml_doc, $ecod_xml_doc, $strip ) = @_;
    my %domain_nodes;
    foreach my $domain_node ( find_domain_nodes($chainwise_xml_doc) ) {
        my ( $uid, $ecod_domain_id ) = get_ids($domain_node);
        if ( has_uniprot_range($domain_node) ) {
            $domain_nodes{$uid} = $domain_node;
        }
    }

    foreach my $domain_node ( find_domain_nodes($ecod_xml_doc) ) {
        my ( $uid, $ecod_domain_id ) = get_ids($domain_node);
        if ( has_uniprot_range($domain_node) ) {
            if ($strip) {
                drop_uniprot_nodes($domain_node);
            } else {
                next;
            }
        }
        if ( $domain_nodes{$uid} ) {
            foreach my $unp_head ( find_uniprot_nodes( $domain_nodes{$uid} ) ) {
                $domain_node->appendChild($unp_head);
            }
        }
    }
}

sub assemble_chainwise_domain_fragments {
    my ($chainwise_xml_doc) = @_;

    #Assemble chainwise domain fragments to mc_domains
    my %mc_domains;
    my %domain_nodes;
    foreach my $domain_node ( find_domain_nodes($chainwise_xml_doc) ) {
        my ( $uid, $ecod_domain_id ) = get_ids($domain_node);
        if ( is_domain_fragment($domain_node) ) {
            $uid =~ /(\d+)./;
            my $parent_uid = $1;
            if ( has_uniprot_range($domain_node) ) {
                push( @{ $mc_domains{$parent_uid} }, $domain_node );
            }
        }
        $domain_nodes{$uid} = $domain_node;
    }

	#Generate lookups
	my %sort;
	foreach my $domain_node ( find_domain_nodes($chainwise_xml_doc) ) { 
		my ( $uid, $ecod_domain_id ) = get_ids($domain_node);
		if ($mc_domains{$uid}) { 
			my ($seqid_range_aref, $chain_range_aref) = multi_chain_range_expand(get_seqid_range($domain_node));
			for (my $i = 0; $i < scalar(@$seqid_range_aref); $i++) { 
				my $chain = $$chain_range_aref[$i];
				my $seqid = $$seqid_range_aref[$i];
				$sort{$uid}{$chain}{$seqid} = $i;
			}
		}
	}

    my ( %unp_ids, %unp_names );
    foreach my $mc_uid ( keys %mc_domains ) {

        my %ranges;
        my %seqid_ranges;
        foreach my $domain_fragment ( @{ $mc_domains{$mc_uid} } ) {
            foreach my $unp_head ( find_uniprot_nodes($domain_fragment) ) {
                my $unp_acc     = get_unp_acc($unp_head);
                my $unp_id      = get_unp_id($unp_head);
                my $unp_name    = get_unp_name($unp_head);
                #my $unp_range   = $unp_head->findvalue('unp_range');
				my $unp_range	 = get_unp_range($unp_head);
                #my $seqid_range = $unp_head->findvalue('seqid_range');
				my $seqid_range	 = get_seqid_range($unp_head);

                $unp_ids{$unp_acc}    = $unp_id;
                $unp_names{$unp_name} = $unp_name;

                if ( !$ranges{$unp_acc} ) {
                    $ranges{$unp_acc}       = range_expand($unp_range);
                    ($seqid_ranges{$unp_acc}{seqid}, $seqid_ranges{$unp_acc}{chain}) = multi_chain_range_expand($seqid_range);
                }
                else {
                    range_include( range_expand($unp_range),   $ranges{$unp_acc} );
                    multi_chain_range_include( 
						multi_chain_range_expand($seqid_range), $seqid_ranges{$unp_acc}{seqid}, $seqid_ranges{$unp_acc}{chain}, $sort{$mc_uid}, 
						);
                }
            }
        }
        foreach my $unp_acc ( keys %ranges ) {
            my $total_range       = rangify( @{ $ranges{$unp_acc} } );
            my $total_seqid_range = multi_chain_rangify(  $seqid_ranges{$unp_acc}{seqid}, $seqid_ranges{$unp_acc}{chain} );

            my $unp_head_node = $chainwise_xml_doc->createElement('uniprot');

            $unp_head_node->setAttribute( 'unp_acc', $unp_acc );
            $unp_head_node->setAttribute( 'id',      $unp_ids{$unp_acc} );
            $unp_head_node->setAttribute( 'name',    $unp_names{$unp_acc} );

            $unp_head_node->appendTextChild( 'unp_range', $total_range );
            $unp_head_node->appendTextChild( 'seqid_range', $total_seqid_range );
            $domain_nodes{$mc_uid}->appendChild($unp_head_node);

        }

    }
}



sub has_any_chain { 
	$_[0]->exists('chain');
}
sub has_chain {
    $_[0]->exists(qq{chain[\@chain_id="$_[1]"]});
}

sub is_domain_fragment {
    $_[0]->findvalue('@domain_fragment') eq 'true';
}

sub get_chain_node {
    $_[0]->findnodes("chain[\@chain_id='$_[1]']")->get_node(1);
}

sub find_pdb_nodes {
    $_[0]->findnodes('.//pdb');
}

sub drop_uniprot_nodes {
    foreach ( $_[0]->findnodes('uniprot') ) {
        $_->unbindNode;
    }
}

sub has_uniprot_mapping {
    $_[0]->exists('.//unp_mapping');
}

sub has_uniprot_range {
    $_[0]->exists('.//unp_range');
}

sub find_unp_range_nodes {
    $_[0]->findnodes('unp_range');
}

sub find_unp_mapping_nodes {
    $_[0]->findnodes('.//unp_mapping');
}

sub get_unp_acc {
    $_[0]->findvalue('.//@unp_acc');
}

sub get_unp_id {
    $_[0]->findvalue('.//@id');
}

sub get_unp_name {
    $_[0]->findvalue('.//@name');
}

sub get_res_range {
    $_[0]->findvalue('residue_range');
}

sub get_unp_pdb_range {
    $_[0]->findvalue('pdb_residue_range');
}

sub get_seq_range {
    $_[0]->findvalue('seq_range');
}
1;
