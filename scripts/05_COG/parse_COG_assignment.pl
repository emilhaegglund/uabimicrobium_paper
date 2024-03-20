#!/usr/bin/perl

open COGNAMES, "/home/dtamarit/cog170403/cognames2003-2014.tab";
while (<COGNAMES>) {
    next if (/^#/);
    @line = split /\t/;
    $cats{$line[0]} = $line[1];
}
close COGNAMES;

open IN, $ARGV[0];
while (<IN>) {
    if (/^(\S+)\s+.*\s+(\S+)$/) {
        $id = $1;
        $cogline = $2;
        $cogline =~ s/;$//;
        @categories = ();
        
        @cogs = split /;/, $cogline;
        foreach $cog (@cogs) {
            @split = split /,/, $cog;
            if ($split[0] =~ /(COG\d+)/) {
                $cat = $cats{$1};
            } elsif ($split[0] =~ /NO_HIT/) {
                $cat = "-";
            } elsif ($split[0] =~ /NOT_ASSIGNED/) {
                $cat = "?";
            } elsif ($split[0] =~ /\d_HITs/) {
                $cat = "+";
            }
            push @categories, $cat;
        }
        @sorted_cat = sort @categories;
        $print_cat = join "", @sorted_cat;
        $print_cat =~ s/(\W+)(\w+)/$2$1/;

        print $id, "\t", $print_cat, "\n"; 
    }
}
