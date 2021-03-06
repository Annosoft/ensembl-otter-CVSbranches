#!/usr/bin/env perl

use strict;
use warnings;

use lib "${ENV{TEAM_TOOLS}}/t/tlib";
use CriticModule;

use Bio::Seq;

use Test::More tests => 2;

my $module;
BEGIN {
    $module = 'Bio::Vega::Utils::Align';
    use_ok($module); 
}

critic_module_ok($module);

1;

# Local Variables:
# mode: perl
# End:

# EOF
