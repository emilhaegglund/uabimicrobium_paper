#! /usr/bin/perl -w
use strict;


use Bio::SeqIO;
use Bio::Seq;
use Bio::Tools::Run::StandAloneBlast;
use Bio::SearchIO;
use Getopt::Long;

my $threads = 8;
GetOptions(
    't|threads:i' => \$threads,
    );


unless ($ARGV[0]) {
  die "\nYou must specify an input file/files\!\n\n\t
Runs blastp with the fastafile against the specified database and \n\t
parses the blast output\n\tor\n\t
Parsing is done directly on the specified blastfile\n\n
COG_class.pl 'fastafile' 'database name'\n
COG_class.pl 'blastfile'\n\n";
}

## ASSOCIATE COGs WITH THEIR CATEGORIES
my %cats = ();
print STDERR "Grabbing COG categories ... \n";
open COGNAMES, "/home/dtamarit/cog170403/cognames2003-2014.tab";
while (<COGNAMES>) {
    next if (/^#/);
    my @line = split /\t/;
    $cats{$line[0]} = $line[1];
}
close COGNAMES;


## IF PROVIDED FASTA AND DATABASE, DO BLAST AND READ REPORT. OTHERWISE, SIMPLY READ REPORT
my @seq_array = ();
my $blast_report = "";
if (defined $ARGV[1])  {
    print STDERR "Blasting $ARGV[0] against $ARGV[1] ...\n";
    my $fasta = Bio::SeqIO->new(-file => "$ARGV[0]",
				-format => 'Fasta');
    while (my $entry = $fasta->next_seq()) {
	push @seq_array, $entry;
    }
    system 'export BLASTDB=$ARGV[1]';
    my $database=$ARGV[1];
    my @params = ('program' => 'blastp',
		  'database' => $database,
		  'outfile' => 'ut.blastp',
		  'F' => 'F',
		  'a' => $threads,
		  '_READMETHOD' => 'Blast');
    my $factory = Bio::Tools::Run::StandAloneBlast->new(@params);
    $blast_report = $factory->blastall(\@seq_array);
} else {
    print STDERR "Parsing blast report $ARGV[0] ....\n";
    $blast_report =Bio::SearchIO->new(-file => "$ARGV[0]",
				      -format => 'blast');
}

my %class = ();
my %times_id = ();
my %descriptions = ();
my $counter = 0;
my $previous_desc=0;

print STDERR "Query number:\n";

QUERY: while( my $result = $blast_report->next_result ) {
    ## SHOW PROGRESS
    $counter++;
    if ($counter % 10 == 0) {
	print STDERR "\r$counter";
    }
    my $region = 0;
 #   my $j=0;
    my $name_q = $result->query_name();
    my $q_length = $result->query_length();

    ## AVOID DUPLICATING IDs EVEN IF THEY WERE PROVIDED
    if ($class{$name_q}) {
	$times_id{$name_q}++;
	print STDERR "Query name $name_q repeated. Using ",$name_q,$times_id{$name_q}, "now\n";
	$name_q = $name_q.$times_id{$name_q};
    }

    ## STORE DESCRIPTIONS FOR OUTPUT
    my $descr_q = $result->query_description();
    $descriptions{$name_q} = $descr_q;
    
    ## NO HITs
    my $num_hits=$result->num_hits;
    if ($num_hits == 0) {
	$class{$name_q}{'1'}{'status'} = "NO_HIT";
    }   

    HIT: while( my $hit = $result->next_hit ) {
#	my $i=0;
#	my $desc = $hit->name();
	my $hit_desc = $hit->description;
	my @hcog = split(",",$hit_desc);
	my $cog = $hcog[1];
	my $hit_length = $hit->length;


	HSP: while( my $hsp = $hit->next_hsp ) {
#	    $i++;

#	    my $length_hit= $hsp->length('hit');
	    my $length_query= $hsp->length('query');
#	    my $start_hit=$hsp->start('hit');
#	    my $end_hit=$hsp->end('hit');
	    my $start_query=$hsp->start('query');
	    my $end_query=$hsp->end('query');
#	    my $percentid = $hsp->percent_identity();
#	    my $short_percent = sprintf "%.0f", $percentid;

	    my $evalue= $hsp->evalue();
	    if ($evalue =~ /^e-\d*/) {
		$evalue =~ s/e/1e/;
	    }
	    
	    next QUERY if ($evalue > 1e-5);

	    
	    if (! defined($class{$name_q})) {
	  ### GENERATE FIRST REGION
		my $region = 1;
		$class{$name_q}{$region}{'start'} = $start_query;
		$class{$name_q}{$region}{'end'} = $end_query;
		$class{$name_q}{$region}{'cog'} = $cog;
		$class{$name_q}{$region}{'status'} = 1;

#		print STDERR "-- $start_query / $end_query / $cog\n";

	    } else {
	  ### ASSIGN TO DEFINED REGION
		## CHECK IF OVERLAP WITH ANY DEFINED REGION IS SIGNIFICANT ($assigned = 1)
		## CHECK IF NO OVERLAP EXISTS WITH ANY DEFINED REGION
		## IGNORE OTHERWISE

		my $assigned = 0;
		my $nooverlap = 0;
	      REGION: foreach my $r (sort {$a <=> $b} keys %{$class{$name_q}}) {
		  next REGION if ($assigned);
		  if ($class{$name_q}{$r}{'status'}) {
		      next REGION if ($class{$name_q}{$r}{'status'} eq "ambiguous");
		      next REGION if ($class{$name_q}{$r}{'status'} eq "classified");
		  }
		  my $start_r = $class{$name_q}{$r}{'start'};
		  my $end_r = $class{$name_q}{$r}{'end'};
		  my $length_r = $end_r-$start_r;
		  
		  if ($start_query >= $start_r && $end_query <= $end_r ) {
		      ## QUERY WITHIN REGION
		      if ($length_query > 0.5*$length_r) {
			  $assigned = 1;
		      }
		  } elsif ($start_query < $start_r && $end_query > $end_r) {
		      ## QUERY CONTAINS REGION
		      if ($length_r > 0.5*$length_query) {
			  $assigned = 1;
		      }
		  } elsif ($start_query < $start_r && $start_r < $end_query && $end_query <= $end_r) {
		      ## QUERY OVERLAPS WITH 3' OF REGION 
		      if ((($end_query-$start_r) > 0.5*$length_r) || 
			  (($end_query-$start_r) > 0.5*$length_query)) {
			  $assigned = 1;
		      }
		  } elsif ($start_query >= $start_r && $start_query < $end_r && $end_r < $end_query) {
		      ## OVERLAPS
		      if ((($end_r-$start_query) > 0.5*$length_r) || 
			  (($end_r-$start_query) > 0.5*$length_query)) {
			  $assigned = 1;
		      }
		  } 
		  else {
		      $nooverlap++;
		  }
		  
		  ## IMPROVE STATUS OF DEFINED REGION IF IT MATCHES, AND VALIDATE IF 3 HITS TO SAME COG
		  ## DECLARE AS AMBIGUOUS IF MATCHES AGAINST DIFFERENT COGS 
		  if ($assigned) {
		      if ($class{$name_q}{$r}{'cog'} eq $cog) {
			  $class{$name_q}{$r}{'status'}++;
			  $class{$name_q}{$r}{'status'} = "classified" if ($class{$name_q}{$r}{'status'} == 5);
		      } else {
			  $class{$name_q}{$r}{'status'} = "ambiguous";
		      }
		  }
	      }
		## IF NO OVERLAP EXISTS TO ANY DEFINED REGION, DECLARE A NEW ONE
		if ($nooverlap == scalar(keys %{$class{$name_q}})) {
		    $region = 1 + keys %{$class{$name_q}};
		    $class{$name_q}{$region}{'start'} = $start_query;
		    $class{$name_q}{$region}{'end'} = $end_query;
		    $class{$name_q}{$region}{'cog'} = $cog;
		    $class{$name_q}{$region}{'status'} = 1;
		    
		}
		
	    }
	}
    }
}


## PRINT RAW OUTPUT
foreach my $prot (sort keys %class) {
    print $prot, "\t", $descriptions{$prot}, "\t";
    foreach my $region (sort {$a <=> $b} keys %{$class{$prot}}) {
	if ($class{$prot}{$region}{'status'} eq "classified") {
	    print $class{$prot}{$region}{'cog'}, "[",$cats{$class{$prot}{$region}{'cog'}}, "],", $class{$prot}{$region}{'start'}, ",", $class{$prot}{$region}{'end'}, ";";
	} elsif ($class{$prot}{$region}{'status'} eq "ambiguous") {
	    print "NOT_ASSIGNED[?]", ",", $class{$prot}{$region}{'start'}, ",", $class{$prot}{$region}{'end'}, ";";
	} elsif ($class{$prot}{$region}{'status'} eq "NO_HIT") {
	    print "NO_HIT[-]", ";";
	} else {
	    print $class{$prot}{$region}{'status'}, "_HITs[+]", ",", $class{$prot}{$region}{'start'}, ",", $class{$prot}{$region}{'end'}, ";";
	}
    }
    print "\n";
}
