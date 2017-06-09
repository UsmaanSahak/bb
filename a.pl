#!/usr/bin/env perl

use LWP::Simple;
use JSON;
use strict;
use warnings;
use CGI;

my $o = new CGI;
print $o->header(-type => "application/json", -charset => "utf-8");

my $e = $o->param("entry");

my $url = 'http://test.borreliabase.org:8983/solr/testcore/suggest?suggest.q='.$e.'&wt=json&suggest=true&suggest.dictionary=suggLo';
my $content = get($url);
die "no url" unless defined $content;
my $result = decode_json($content);

my $documents = $result->{suggest}->{suggLo}->{$e}->{suggestions};
my $numDocs = @$documents;

print "<table position:relative>";
for (my $i = 1; $i < $numDocs; $i++) {
  print "<tr><th>$documents->[$i]->{term}</th></tr>";
}
print "</table>";
