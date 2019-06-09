#!/usr/bin/perl

use lib "$ENV{HOME}/bin/bioperl-live";
use strict;
use warnings;
use CGI;
use Bio::Graphics;
use Bio::SearchIO;
use Bio::SeqFeature::Generic;

my $cgi = CGI->new();
print $cgi->header(-type => 'image/png');


my %opts=(
	'chr'		=> "chr01",
	'gaps'		=> 0,
	'start' 	=> 1,
	'end'		=> 0,
	'width'		=> 2048,
	'length'	=> 0,
	'sname'		=> "",
	'format'	=> "png",
	'filter'	=> "",
	'maplen'	=> 0,
	'lenidx'	=> 5,
	'queidx'	=> 0,
	'subidx'	=> 2,
	'scoidx'	=> 4,
	'staidx'	=> 8,
	'endidx'	=> 9,
	'gname'		=> 1,
	'pname'		=> 1,
	'galgn'		=> 0,
);
my %param = map { $_ => scalar $cgi->param($_) } $cgi->param() ;
foreach my $k (keys %param) {
	$opts{$k} = $param{$k};
}

my @gaps=();
if (! $opts{'gaps'}) {
	open(my $fh, "<", "res/$opts{'chr'}.gaps");
	while(my $ln = <$fh>) {
		$ln =~ /^\((\d+)\,\s(\d+)\)/;
		push @gaps, Bio::SeqFeature::Generic->new(
			-start 			=> $1,
			-end 			=> $2,
			-strand 		=> 1,
			-id				=> $1,
			-display_name	=> $opts{"gname"}?$1:"",
		) if $ln =~ /^\((\d+)\,\s(\d+)\)/;
	}
	close $fh;
}

my %m = {};
if ($opts{'pname'}) {
	open(my $fh, "<", "res/seqmap.xls");
	while(my $ln = <$fh>) {
		chomp $ln;
		my @n = split "\t", $ln;
		$m{$n[2]}="$n[0] $n[1] $n[2]";
	}
	close $fh;
}

open(my $fh, "<", "res/$opts{'chr'}.xls") || die "Can not open BLAST file!\n";
my %res= {};
while(my $ln = <$fh>) {
	my @cls = split /\t/, $ln;
	next if $opts{'filter'} && $cls[$opts{'queidx'}]=~/$opts{'filter'}/ || $cls[$opts{'lenidx'}]<$opts{'maplen'};
	$opts{'sname'} = $cls[2] if !$opts{'sname'};
	$opts{'length'} = $cls[$opts{'staidx'}] if $opts{'length'} < $cls[$opts{'staidx'}];
	$opts{'length'} = $cls[$opts{'endidx'}] if $opts{'length'} < $cls[$opts{'endidx'}];
	$res{$cls[0]} = Bio::SeqFeature::Generic->new(
		-score        => $cls[2],
    	-display_name => $opts{'pname'}?$m{$cls[0]}:$cls[0],
	) if !$res{$cls[0]};
	my ($start, $end) = sort($cls[$opts{'staidx'}], $cls[$opts{'endidx'}]);
	my $sft = Bio::SeqFeature::Generic->new(
		-start => $start,
		-end => $end,
		-score => $cls[$opts{'scoidx'}],
	);
	if ($opts{'galgn'}) {
		my $flag = 0;
		foreach my $g (@gaps) {
			$flag = 1 if $g->start > $sft->start && $g->end < $sft->end;
		}
		$res{$cls[0]}->add_sub_SeqFeature($sft, 'EXPAND') if $flag;
	} else { 
		$res{$cls[0]}->add_sub_SeqFeature($sft, 'EXPAND');
	}
}
close $fh;
$opts{'end'} = $opts{'length'} if $opts{'end'} <= 0;


my $panel = Bio::Graphics::Panel->new(
									 -grid		=> 1,
									 -start		=> $opts{'start'},
                                     -end 	  	=> $opts{'end'},
                                     -width     => $opts{'width'},
                                     -pad_left  => 50,
                                     -pad_right => 50,
                                    );

my $full_length = Bio::SeqFeature::Generic->new(
                                    -start        => 1,
                                    -end          => $opts{'length'},
                                    -display_name => $opts{'sname'},
                                              );

$panel->add_track($full_length,
                 -glyph   => 'arrow',
                 -tick    => 2,
                 -fgcolor => 'black',
                 -double  => 1,
                 -label   => 1,
                );

my $gtrack1 = $panel->add_track(
		-glyph       => 'box',
		-label		 => 1,
		-bgcolor     => 'red',
		-font2color  => 'red',
);
$gtrack1->add_feature(@gaps);

my $ftrack = $panel->add_track(
            -glyph       => 'graded_segments',
            -label       => 1,
            -connector   => 'dashed',
            -bgcolor     => 'blue',
            -font2color  => 'red',
			-font		 => 'gdGiantFont',
            -description => sub {
                my $feature = shift;
                return unless $feature->has_tag('description');
                my ($description) = $feature->each_tag_value('description');
                my $score = $feature->score;
                "$description, score=$score";
            },
        );
foreach my $key (sort {$b cmp $a} keys %res) {
    $ftrack->add_feature($res{$key});
}
my $gtrack2 = $panel->add_track(
		-glyph       => 'box',
		-label		 => 1,
		-bgcolor     => 'red',
		-font2color  => 'red',
);
$gtrack2->add_feature(@gaps);

print $panel->png;
