#!/usr/bin/perl -w
use strict;

# Find the tests in this directory.
my @tests = split(/\s+/, `ls *.cc`);
my $goodTests = 0;
my $badTests = 0;

# Write the Makefile.
system("rm -f Makefile");
system("cat Makefile.head > Makefile");

open(MAKEFILE, ">>Makefile");
print MAKEFILE "\n";

foreach my $test (@tests) {
  chomp $test;
  $test =~ s/\.cc//;

  my $uc_test = uc($test);

  print MAKEFILE "$uc_test" . "_OBJS = $test.o\n";
  print MAKEFILE "$uc_test" . "_EXE = $test\n";
  print MAKEFILE "\n";
}

print MAKEFILE "EXECUTABLES = " . join(" ", @tests) . "\n\n";
print MAKEFILE "default : \$(EXECUTABLES)\n\n";

foreach my $test (@tests) {
  my $uc_test = uc($test);

  print MAKEFILE "$test : \$($uc_test" . "_OBJS)\n";
  print MAKEFILE "\t\$(CXX) -o \$($uc_test" . "_EXE) \$($uc_test" .
                 "_OBJS) \$(LDFLAGS)\n\n";
}

system("cat Makefile.tail >> Makefile");

# Now run the newly-created makefile.
system("make clean");
if (system("make -j 4") != 0) {
  system("rm -f Makefile");
  die "Compilation failed.\n";
}

# If build succeeded, run the tests.
print "\n";
print "===================\n";
foreach my $test (@tests) {
  chomp $test;

  system("./$test >& $test.result");

  open(RESULT_FILE, "$test.result");
  my @lines = <RESULT_FILE>;
  close(RESULT_FILE);

  if (@lines > 0) {
    print "FAILED: $test\n";
    $badTests++;
  } else {
    print "PASSED: $test\n";
    $goodTests++;
  }
}

print "$goodTests of " . @tests . " tests passed.\n";
print "===================\n\n";

# Clean up -- remove the outputs and the Makefile.
system("make clean");
system("rm -f Makefile");
