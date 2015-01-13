#!/usr/bin/perl

use warnings;
use strict;

use Getopt::Long;
use FileHandle;

=head1 VERSION

    Version 1.01

=head1 SYNOPSIS

    mergeSegments is an accessory script to VarScan 2 (http://varscan.sourceforge.net) that merges segments from the DNAcopy library and classifies them as focal/large-scale amplification/deletions based on user-specified thresholds.

=head1 USAGE    perl mergeSegments.pl [segments] OPTIONS

=head3  OPTIONS:

        --amp-threshold     Threshold above which a segment will be considered an amplification [0.25]

        --del-threshold     Threshold below which a segment will be considered a deletion [-0.25]

        --size-threshold    Fraction of a chromosome arm above which an event is considered large-scale [0.25]

        --ref-arm-sizes     Two column file of reference name and size in bp for calling by chromosome arm

        --output-basename   The basename for output files. Two will be created: basename.events.tsv and basename.summary.tsv

        --verbose           If set to 1, provide verbose output

=cut






our $VERSION = '1.01';

my $usage = qq{USAGE: mergeSegments.pl [segments] OPTIONS
        segments - A segments file with p-values from the DNAcopy library. This should be tab- or space-delimited
        with a header and the following columns: chrom, loc.start, loc.end, num.mark, seg.mean, bstat, pval, lcl, ucl.
        OPTIONS:
        --amp-threshold     Threshold above which a segment will be considered an amplification [0.25]
        --del-threshold     Threshold below which a segment will be considered a deletion [-0.25]
        --size-threshold    Fraction of a chromosome arm above which an event is considered large-scale [0.25]
        --ref-arm-sizes     Two column file of reference name and size in bp for calling by chromosome arm
        --output-basename   The basename for output files. Two will be created: basename.events.tsv and basename.summary.tsv
        --verbose           If set to 1, provide verbose output
};

die $usage if(!$ARGV[0]);
## Parse arguments ##

my $amp_threshold = 0.25;
my $del_threshold = -0.25;
my $size_threshold = 0.25;
my $ref_arm_sizes = "";
my $output_basename = "outfile";
my $verbose = 1;
my $ret = parse_arguments();
my %ref_sizes = ();

execute();

################################################################################

=head2	parse_arguments

    Parse user-defined options.

=cut
################################################################################

sub parse_arguments
{

    my $result = GetOptions (
                                "amp-threshold=s"   => \$amp_threshold,
                                "del-threshold=s"   => \$del_threshold,
                                "size-threshold=s"   => \$size_threshold,
                                "ref-arm-sizes=s"   => \$ref_arm_sizes,
                                "output-basename=s"   => \$output_basename,
                                "verbose=s"   => \$verbose,
    );    

}




################################################################################

=head2	execute 

    Main program execution.


=cut


sub execute
{
    my %stats = ();
    $stats{'num_variants'} = $stats{'num_merged_events'} = $stats{'num_fail_pos'} = $stats{'num_fail_strand'} = $stats{'num_fail_varcount'} = $stats{'num_fail_varfreq'} = $stats{'num_fail_mmqs'} = $stats{'num_fail_var_mmqs'} = $stats{'num_fail_mapqual'} = $stats{'num_fail_readlen'} = $stats{'num_fail_dist3'} = $stats{'num_pass_filter'} = 0;
    ## Load the ref sizes ##

    unless($ref_arm_sizes){
   	print "\nWARNING: Please give ref-arm-sizes file\n\n$usage\n";
	exit(1);
    }
    %ref_sizes = parse_ref_sizes($ref_arm_sizes);

	## Parse the segments file ##

	my $input = new FileHandle ($ARGV[0]);
	my $lineCounter = 0;

	my @merged_events = ();


	my $current_chrom = my $current_chr_start = my $current_chr_stop = my $current_event_type = "";
	my $current_list = "";
	my $current_segments = 0;
	
	while (<$input>)
	{
		chomp;
		my $line = $_;
		$lineCounter++;
		
                $line =~ s/\"//g;
		my ($chrom, $chr_start, $chr_stop, $num_mark, $seg_mean, $bstat, $p_value, $lcl, $ucl) = split(/\s+/, $line);

		if($chrom eq "chrom")
		{
			## Skip header ##
		}
		else
		{
			$chrom = "X" if($chrom eq "23");
			$chrom = "Y" if($chrom eq "24");
			$chrom =~ s/\"//g;		
			if (!exists($ref_sizes{$chrom})){
				print "ERROR: cannot find chromosome $chrom in $ref_arm_sizes, quitting\n";
				exit(1);
			}
			## Determine event type ##
			my $event_type = "neutral";		
	
			if($seg_mean >= $amp_threshold)
			{
				$event_type = "amplification";
				$stats{'num_amp_del_segments'}++;
			}
			elsif($seg_mean <= $del_threshold)
			{
				$event_type = "deletion";
				$stats{'num_amp_del_segments'}++;
			}
	
			## PROCEED ON JOINING DECISION ##
	
			## Case 1: No Events Yet Parsed ##
	
			if(!$current_chrom)
			{
				## Start a new event ##
				$current_chrom = $chrom;
				$current_chr_start = $chr_start;
				$current_chr_stop = $chr_stop;
				$current_event_type = $event_type;
				$current_segments = 1;
				$current_list = $line;
			}
			
			## Case 2: Chromosome or event type changes ##
			
			elsif ($chrom ne $current_chrom || $event_type eq "neutral" || $event_type ne $current_event_type)
			{
				if($current_segments)
				{
                                    $current_list = "" if(!$current_list);
					$merged_events[$stats{'num_merged_events'}] = process_event($current_chrom, $current_chr_start, $current_chr_stop, $current_event_type, $current_segments, $current_list);
					$stats{'num_merged_events'}++;		
				}
				
				if($event_type ne "neutral")
				{
					## Start new region ##
					
					$current_chrom = $chrom;
					$current_chr_start = $chr_start;
					$current_chr_stop = $chr_stop;
					$current_event_type = $event_type;
					$current_segments = 1;
					$current_list = $line;		
				}
				else
				{
					$current_chrom = $current_chr_start = $current_chr_stop = $current_event_type = $current_list = "";
					$current_segments = 0;
				}
			}
			
			## Case 3: Same type of event ##
			
			elsif($event_type eq $current_event_type)
			{
				$current_chr_stop = $chr_stop;
				$current_list .= "\n" . $line;
				$current_segments++;
			}
		}
	}
	
	close($input);
	
	if($current_segments)
	{
		$merged_events[$stats{'num_merged_events'}] = process_event($current_chrom, $current_chr_start, $current_chr_stop, $current_event_type, $current_segments, $current_list);
		$stats{'num_merged_events'}++;		
	}


#	print "$stats{'num_amp_del_segments'} copy-change segments\n";
#	print "$stats{'num_merged_events'} merged events\n";	

	my $large_scale_amp_arms = my $large_scale_del_arms = "";

	## Open outfile for amps and dels ##
	
	open(EVENTS, ">$output_basename.events.tsv") or die "Can't open outfile: $!\n";
	print EVENTS "chrom\tchr_start\tchr_stop\tseg_mean\tnum_segments\tnum_markers\tp_value\tevent_type\tevent_size\tsize_class\tchrom_arm\tarm_fraction\tchrom_fraction\n";
	## Process the merged segments ##

	foreach my $line (@merged_events)
	{
		my ($chrom, $chr_start, $chr_stop, $avg_seg_mean, $event_type, $num_segments, $num_mark, $p_value) = split(/\t/, $line);
		my $event_size = $chr_stop - $chr_start + 1;

		## Save event type ##

		$stats{$event_type}++;

		## Determine size category ##
		
		my $size_category = "focal";
		my $event_detail = "";	
		my $arm_name = my $chrom_arm_fraction = my $chrom_fraction = "";

		## Find most-affected chromosome arm ##

		my $p_arm_fraction = my $q_arm_fraction = 0;
                $p_arm_fraction = calculate_arm_fraction($ref_sizes{"$chrom\tp"}, $chr_start, $chr_stop) if($ref_sizes{"$chrom\tp"});
		$q_arm_fraction = calculate_arm_fraction($ref_sizes{"$chrom\tq"}, $chr_start, $chr_stop) if($ref_sizes{"$chrom\tq"});
		my ($q_start, $q_stop) = split(/\t/, $ref_sizes{"$chrom\tq"});
		my $chrom_size = $q_stop;

		if($q_arm_fraction > $p_arm_fraction)
		{
			$arm_name = $chrom . "q";
			$chrom_arm_fraction = $q_arm_fraction;
		}
		else
		{
			$arm_name = $chrom . "p";
			$chrom_arm_fraction = $p_arm_fraction;
		}

		## Determine fraction of chromosome affected ##

		$chrom_fraction = $event_size / $chrom_size;

		## If >50% of arm or >25% of chromosome affected, call it large-scale ##

		if(($chrom_arm_fraction && $chrom_arm_fraction >= $size_threshold) || ($chrom_fraction && $chrom_fraction >= ($size_threshold / 2)))
		{
			$size_category = "large-scale";
			$event_detail = join("\t", $chrom . $arm_name, sprintf("%.2f", $chrom_arm_fraction * 100) . "%", sprintf("%.2f", $chrom_fraction * 100) . "%");

		}
		
		## Save the event ##
		$stats{"$size_category $event_type"}++;

		## Save which arm ##
		$stats{"$size_category $event_type arms"} .= "," if($stats{"$size_category $event_type arms"});
		$stats{"$size_category $event_type arms"} .= $arm_name;

		$chrom_arm_fraction = sprintf("%.2f", $chrom_arm_fraction * 100) . '%';
		$chrom_fraction = sprintf("%.2f", $chrom_fraction * 100) . '%';
		print EVENTS join("\t", $chrom, $chr_start, $chr_stop, $avg_seg_mean, $num_segments, $num_mark, $p_value, $event_type, $event_size, $size_category, $arm_name, $chrom_arm_fraction, $chrom_fraction) . "\n";
	}

	close(EVENTS);



	## Print summary stats to file ##
	$stats{'large-scale amplification arms'} = "NA" if(!$stats{'large-scale amplification arms'});
	$stats{'large-scale deletion arms'} = "NA" if(!$stats{'large-scale deletion arms'});
	
            $stats{'amplification'} = 0 if(!$stats{'amplification'});
            $stats{'large-scale amplification'} = 0 if(!$stats{'large-scale amplification'});
            $stats{'focal amplification'} = 0 if(!$stats{'focal amplification'});
            $stats{'deletion'} = 0 if(!$stats{'deletion'});
            $stats{'large-scale deletion'} = 0 if(!$stats{'large-scale deletion'});
            $stats{'focal deletion'} = 0 if(!$stats{'focal deletion'});       
        
	open(SUMMARY, ">$output_basename.summary.txt") or die "Can't open outfile: $!\n";
	print SUMMARY "segments\tmerged_events\tamps\tlarge-scale\tfocal\tdels\tlarge-scale\tfocal\tamp_regions\tdel_regions\n";

	print SUMMARY join("\t", $stats{'num_amp_del_segments'}, $stats{'num_merged_events'}, $stats{'amplification'}, $stats{'large-scale amplification'}, $stats{'focal amplification'}, $stats{'deletion'}, $stats{'focal deletion'}, $stats{'large-scale deletion'}, $stats{'large-scale amplification arms'}, $stats{'large-scale deletion arms'}) . "\n";

	close(SUMMARY);

	if($verbose)
	{

            
		print $stats{'amplification'} . " classified as amplifications\n";
		print $stats{'large-scale amplification'} . " large-scale amplifications ";
		print "(" . $stats{'large-scale amplification arms'} . ")" if($stats{'large-scale amplification'});
		print "\n";
		print $stats{'focal amplification'} . " focal amplifications\n";
	
		print $stats{'deletion'} . " classified as deletions\n";
		print $stats{'large-scale deletion'} . " large-scale deletions ";
		print "(" . $stats{'large-scale deletion arms'} . ")" if($stats{'large-scale deletion'});
		print "\n";
		print $stats{'focal deletion'} . " focal deletions\n";				
	}


    return(0);
}




################################################################################

=head3	process_event

    Process a merged event


=cut

sub process_event
{
	my ($event_chrom, $event_chr_start, $event_chr_stop, $event_type, $num_segments, $segment_list) = @_;

	## Split merged events and calculate averages ##
	
	my @events = split(/\n/, $segment_list);
	
	my $total_num_mark = my $seg_mean_sum = my $p_val_sum = my $p_val_num = 0;
	
	foreach my $line (@events)
	{
		my ($id, $sample, $chrom, $chr_start, $chr_stop, $num_mark, $seg_mean, $bstat, $p_value, $lcl, $ucl) = split(/\s+/, $line);		
		$total_num_mark += $num_mark;
		$seg_mean_sum += $seg_mean;
		$p_val_sum += $p_value if($p_value ne "NA");
		$p_val_num++;
	}

	## Calculate averages ##
	
	my $avg_seg_mean = $seg_mean_sum / $num_segments;
	my $avg_p_value = "NA";
	$avg_p_value = $p_val_sum / $p_val_num if($p_val_num);

	## Calculate region_size ##
	
	my $event_size = $event_chr_stop - $event_chr_start + 1;

	return(join("\t", $event_chrom, $event_chr_start, $event_chr_stop, $avg_seg_mean, $event_type, $num_segments, $total_num_mark, $avg_p_value));
	
}



################################################################################

=head3	calculate_arm_fraction

    Compute the fraction of a chromosome (arm) encompassed by an event


=cut

sub calculate_arm_fraction
{
	my ($ref_sizes, $chr_start, $chr_stop) = @_;
	my ($arm_start, $arm_stop) = split(/\t/, $ref_sizes);
	my $arm_size = $arm_stop - $arm_start + 1;

	## Check to see if they overlap ##
	
	if($chr_start <= $arm_stop && $chr_stop >= $arm_start)
	{
		## Determine bounds of this event on this chromosome arm ##
		my $event_arm_start = $chr_start;
		my $event_arm_stop = $chr_stop;
		$event_arm_start = $arm_start if($event_arm_start < $arm_start);
		$event_arm_stop = $arm_stop if($event_arm_stop > $arm_stop);
		my $event_arm_size = $event_arm_stop - $event_arm_start + 1;				

		if($event_arm_size)
		{
			my $fraction_of_arm = $event_arm_size / $arm_size;
		}
		else
		{
			return(0);
		}
	}
	else
	{
		return(0);
	}
}




################################################################################

=head3	parse_ref_sizes

    Load the sizes of chromosome arms


=cut

sub parse_ref_sizes
{                               # replace with real execution logic.
	my $FileName = shift(@_);
	my %arms = ();

	my $input = new FileHandle ($FileName);
	my $lineCounter = 0;
	
	while (<$input>)
	{
		chomp;
		my $line = $_;
		$lineCounter++;
		
		my ($chrom, $chr_start, $chr_stop, $arm_name) = split(/\s+/, $line);
		$arms{$chrom . "\t" . $arm_name} = join("\t", $chr_start, $chr_stop);
		$arms{$chrom} = 1;
	}
	
	close($input);
	
	return(%arms);
}



=head1 AUTHOR

    Daniel C. Koboldt, << <dkoboldt at genome.wustl.edu> >>
    The Genome Institute at Washington University School of Medicine
    St. Louis, Missouri, USA

=head1 COPYRIGHT

    Copyright 2009-2012 Daniel C. Koboldt and Washington University
    All rights reserved.

=head1 LICENSE

    This program is free for non-commercial use. Please contact the author for commercial licensing information.

=cut
