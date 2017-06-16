#!/usr/bin/env perl

use LWP::Simple;
use JSON;
use strict;
use warnings;
use CGI;

my $o = new CGI;
print $o->header(-type => "application/json", -charset => "utf-8");

#my $e = "BB";
my $e = $o->param("entry");
#my $s = "loSuggestions";
my $s = $o->param("suggType");

my ($url,$content,$result,$documents,$numDocs);
if ($s eq "loSuggestions") {
 $url = 'http://test.borreliabase.org:8983/solr/testcore/suggest?suggest.q='.$e.'&wt=json&suggest=true&suggest.dictionary=suggLo';

 $content = get($url);
 die "no url" unless defined $content;
 $result = decode_json($content);
 $documents = $result->{suggest}->{suggLo}->{$e}->{suggestions};
 $numDocs = @$documents;
} 
if ($s eq "symSuggestions") {
 $url = 'http://test.borreliabase.org:8983/solr/testcore/suggest?suggest.q='.$e.'&wt=json&suggest=true&suggest.dictionary=suggSym';
 $content = get($url);
 die "no url" unless defined $content;
 $result = decode_json($content);
 $documents = $result->{suggest}->{suggSym}->{$e}->{suggestions};
 $numDocs = @$documents;
}


print "<table position:relative>";
print "<tr><th> $s is suggTypes </th></tr>";

for (my $i = 1; $i < $numDocs; $i++) {
  print "<tr><th>$documents->[$i]->{term}</th></tr>";
}
print "</table>";
